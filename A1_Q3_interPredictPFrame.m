function [predictedFrame, reconstructedFrame] = A1_Q3_interPredictPFrame(reference_frames, currentFrame, searchRange, blockSize, paddedHeight, paddedWidth, n, QP, MDiff_stream, QTC_stream, nRefFrames, MRFoverlay, frameIdx)
    Q_Matrix = A1_Q4_generateQMatrix(blockSize, QP);
    previous_mv = [0, 0];
    previous_ref_index = 0;

    colors = [
        255, 0, 0;    % Red for reference frame 1
        0, 255, 0;    % Green for reference frame 2
        0, 0, 255;    % Blue for reference frame 3
        255, 255, 0;  % Yellow for reference frame 4
    ];

    % Initialize overlay frame (3 channels for RGB overlay)
    overlayFrame = zeros(paddedHeight, paddedWidth, 3, 'double');

    % Loop over blocks
    for row = 1:blockSize:paddedHeight
        for col = 1:blockSize:paddedWidth
            % Boundary check
            if row + blockSize - 1 > paddedHeight || col + blockSize - 1 > paddedWidth
                continue;
            end

            % Current block
            currentBlock = currentFrame(row:row+blockSize-1, col:col+blockSize-1);

            % parameter for MFR
            best_mae = Inf;
            best_ref_index = 1;
            num_ava_refs = min(length(reference_frames), nRefFrames);

            bestMatch = [0, 0];
            
            % loop over reference frames for each current block
            for ref_idx = 1:num_ava_refs
                referenceFrame = reference_frames{ref_idx};

                % find the best match mv and predicted block
                [best_mv, predictedBlock, min_mae] = InterPredictBLK(referenceFrame, currentBlock, row, col, searchRange, blockSize, paddedWidth, paddedHeight);

                if min_mae < best_mae || ...
                   (min_mae == best_mae && abs(best_mv(1)) + abs(best_mv(2)) < abs(bestMatch(1)) + abs(bestMatch(2))) || ...
                   (min_mae == best_mae && abs(best_mv(1)) + abs(best_mv(2)) == abs(bestMatch(1)) + abs(bestMatch(2)) && ...
                   (best_mv(2) < bestMatch(2) || (best_mv(2) == bestMatch(2) && bestMatch(1) < bestMatch(1))))
                    best_mae = min_mae;
                    bestMatch = best_mv;
                    best_ref_index = ref_idx;
                    best_predictedBlk = predictedBlock;
                end
            end

            
            fprintf('best referencing frame idx: %d\n', best_ref_index);
            diff_ref_index = best_ref_index - previous_ref_index;
            previous_ref_index = best_ref_index;

            % Update predicted frame for next iteration
            predictedFrame(row:row+blockSize-1, col:col+blockSize-1) = best_predictedBlk;

            % update the motion vector differentiation
            diff_mv = bestMatch - previous_mv;
            previous_mv = bestMatch;

            % generate MDiff with MV
            fprintf(MDiff_stream, '%s %s %s %s\n', A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(diff_ref_index), A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)));
            
            % find Residual block
            approxResidualBlock = A1_Q3_calcResidual(best_predictedBlk, currentBlock, n);

            % encode the residual block
            encoded_residual_block = A1_Q4_quantizeBlockAfterDCT(approxResidualBlock, Q_Matrix);
            % generate QTC stream
            scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
            rle_encoded = A1_Q4_rleEncode(scanned_coeffs, blockSize);
            %QTC_stream = [QTC_stream, A1_Q4_expGolombEncode(rle_encoded)];
            encoded_value = A1_Q4_expGolombEncode(rle_encoded);
            fprintf(QTC_stream, '%s\n', encoded_value);   % each line is a block

            decoded_residual_block = A1_Q4_idctAfterDequantizeBlock(encoded_residual_block, Q_Matrix);

            % Reconstruct block and update the reconstructed frame
            reconstructedBlock = double(best_predictedBlk) + decoded_residual_block;
            %reconstructedBlock = max(min(reconstructedBlock, 255), 0);
            reconstructedFrame(row:row+blockSize-1, col:col+blockSize-1) = reconstructedBlock;

            if MRFoverlay == 1
                overlayColor = colors(best_ref_index, :);
                for channel = 1:3
                    overlayFrame(row:row+blockSize-1, col:col+blockSize-1, channel) = ...
                        (0.5 * reconstructedFrame(row:row+blockSize-1, col:col+blockSize-1)) + ...
                        (0.5 * overlayColor(channel)); % Alpha blending (50%)
                end
            end
            
        end
    end

    if MRFoverlay == 1
        % Clamp overlay frame values to valid range [0, 255]
        overlayFrame = uint8(max(min(overlayFrame, 255), 0));

        % Display the overlay frame
        % imshow(overlayFrame);

        outputFolder = 'output_frames';  % Directory to save the images
        % Create output folder if it doesn't exist
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end
        
        
        % Save the overlay image as a PNG file
        outputFileName = fullfile(outputFolder, sprintf('frame_overlay_%03d.png', frameIdx));
        imwrite(overlayFrame, outputFileName);
    
        fprintf('Saved overlay frame %d as %s\n', frameIdx, outputFileName);
    end
end

function [bestMatch, predictedBlock, bestMAE] = InterPredictBLK(referenceFrame, currentBlock, row, col, searchRange, blockSize, paddedWidth, paddedHeight)
    bestMAE = inf;
    bestMatch = [0, 0];

    % Full search within search range
    for yOffset = -searchRange:searchRange
        for xOffset = -searchRange:searchRange
            refRow = row + yOffset;
            refCol = col + xOffset;

            if refRow < 1 || refCol < 1 || refRow+blockSize-1 > paddedHeight || refCol+blockSize-1 > paddedWidth
                continue;
            end

            % Reference block
            referenceBlock = referenceFrame(refRow:refRow+blockSize-1, refCol:refCol+blockSize-1);
            cur_mae = mean(abs(double(currentBlock(:)) - double(referenceBlock(:))));

            % Update best match based on MAE and motion vector criteria
            if cur_mae < bestMAE || ...
               (cur_mae == bestMAE && abs(xOffset) + abs(yOffset) < abs(bestMatch(1)) + abs(bestMatch(2))) || ...
               (cur_mae == bestMAE && abs(xOffset) + abs(yOffset) == abs(bestMatch(1)) + abs(bestMatch(2)) && ...
               (yOffset < bestMatch(2) || (yOffset == bestMatch(2) && xOffset < bestMatch(1))))
                bestMAE = cur_mae;
                bestMatch = [xOffset, yOffset];
            end
        end
    end

    % Predicted block based on best match
    bestRow = row + bestMatch(2);
    bestCol = col + bestMatch(1);
    predictedBlock = referenceFrame(bestRow:bestRow+blockSize-1, bestCol:bestCol+blockSize-1);
end


% function yuvImageWithOverlay = addColorOverlayToYUV(yuvImage, alpha)
%     % Ensure input YUV image is double for calculations
%     yuvImage = double(yuvImage);
% 
%     colors = [
%         255, 0, 0;  % Red for reference frame 1
%         0, 255, 0;  % Green for reference frame 2
%         0, 0, 255;  % Blue for reference frame 3
%         255, 255, 0 % Yellow for reference frame 4 (if used)
%     ]; 
% 
%     % Get the color for this reference frame
%     color = colors(refFrameIndex, :);
% 
%     % Convert the overlay color from RGB to YUV
%     transformationMatrix = [
%         0.299,  0.587,  0.114;
%        -0.14713, -0.28886,  0.436;
%         0.615,  -0.51499, -0.10001
%     ];
%     offset = [0; 128; 128];
%     overlayColorYUV = transformationMatrix * color(:) + offset;
% 
%     % Separate Y, U, and V channels
%     Y = yuvImage(:, :, 1);
%     U = yuvImage(:, :, 2);
%     V = yuvImage(:, :, 3);
% 
%     % Blend the U and V channels with the overlay color
%     U_overlay = (1 - alpha) * U + alpha * overlayColorYUV(2);
%     V_overlay = (1 - alpha) * V + alpha * overlayColorYUV(3);
% 
%     % Combine the channels back into a YUV image
%     yuvImageWithOverlay = cat(3, Y, U_overlay, V_overlay);
% 
%     % Clamp the U and V channels to valid range [0, 255]
%     yuvImageWithOverlay(:, :, 2:3) = max(min(yuvImageWithOverlay(:, :, 2:3), 255), 0);
% 
%     % Convert back to uint8 for output
%     yuvImageWithOverlay = uint8(yuvImageWithOverlay);
% end