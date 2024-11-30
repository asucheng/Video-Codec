function [predictedFrame, reconstructedFrame] = A2_interPredictForPFrame(reference_frames, currentFrame, ...
    searchRange, blockSize, paddedHeight, paddedWidth, n, QP, ...
    MDiff_stream, MVPDiff_stream, QTC_stream, nRefFrames, frameIdx, ...
    VBSEnable, MRFoverlay, FMEEnable, FastME)

    predictedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
    reconstructedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
 
    Q_Matrix = A2_generateQMatrix(blockSize, QP);
    lambda = A2_getLambda(QP);

    previous_mv = [0, 0];
    mv = [0, 0];

    % determine which stream to store results
    if FastME
        Diff_stream = MVPDiff_stream;
    else
        Diff_stream = MDiff_stream;
    end

    % parameter for MRF
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
            block_height = min(blockSize, paddedHeight - row + 1);
            block_width = min(blockSize, paddedWidth - col + 1);

            % Boundary check
            if block_height < blockSize || block_width < blockSize
                continue;
            end

            % Current block
            currentBlock = currentFrame(row:row+block_height-1, col:col+block_width-1);
            [bestMatch_ns, predictedBlock_ns, best_ref_index_ns] = A2_interPredictForPBlock(reference_frames, ...
                    currentBlock, row, col, searchRange, blockSize, paddedWidth, paddedHeight, ...
                    nRefFrames, FMEEnable, FastME, mv);
            if VBSEnable
                % Calculate costs for non-split case
                % [bestMatch_ns, predictedBlock_ns, best_ref_index_ns] = A2_interPredictForPBlock(reference_frames, ...
                %     currentBlock, row, col, searchRange, blockSize, paddedWidth, paddedHeight, ...
                %     nRefFrames, FMEEnable, FastME, mv);

                residualBlock_ns = A1_Q3_calcResidual(predictedBlock_ns, currentBlock, n);
                [J_ns, encoded_block_ns] = A2_computeRD(residualBlock_ns, 0, QP, lambda);

                diff_ref_index_ns = best_ref_index_ns - previous_ref_index;

                % Add MV cost to non-split RD cost
                diff_mv_ns = bestMatch_ns - previous_mv;
                mv_bits_ns = strlength(A1_Q4_expGolombEncode(diff_mv_ns(1))) + strlength(A1_Q4_expGolombEncode(diff_mv_ns(2))) + strlength(A1_Q4_expGolombEncode(diff_ref_index_ns));
                J_ns = J_ns + lambda * mv_bits_ns;

                % Calculate split case (4 sub-blocks)
                split_size = blockSize/2;
                J_split = 0;
                relative_offsets = [0, 0; 0, split_size; split_size, 0; split_size, split_size];
                sub_mvs = zeros(4, 2);
                sub_ref_idxs = zeros(4, 1);   % might have a problem
                sub_predicted = zeros(blockSize);
                sub_reconstructed = zeros(blockSize);
                sub_encoded_blocks = cell(4, 1);
                temp_prev_mv = previous_mv;
                temp_prev_ref_idx = previous_ref_index;

                % Process sub-blocks in Z-order
                for idx = 1:4
                    relSubRow = relative_offsets(idx, 1) + 1;
                    relSubCol = relative_offsets(idx, 2) + 1;
                    absSubRow = row + relative_offsets(idx, 1);
                    absSubCol = col + relative_offsets(idx, 2);

                    sub_block = currentFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1);
                    [sub_mv, sub_pred, sub_ref_idx] = A2_interPredictForPBlock(reference_frames, ...
                        sub_block, absSubRow, absSubCol, searchRange, split_size, paddedWidth, paddedHeight, ...
                        nRefFrames, FMEEnable, FastME, mv);
                    sub_residual = A1_Q3_calcResidual(sub_pred, sub_block, n);
                    [sub_J, sub_encoded] = A2_computeRD(sub_residual, 0, QP, lambda);
 
                    % Add MV cost
                    diff_sub_ref_idx = sub_ref_idx - temp_prev_ref_idx;
                    diff_sub_mv = sub_mv - temp_prev_mv;
                    mv_bits = strlength(A1_Q4_expGolombEncode(diff_sub_mv(1))) + strlength(A1_Q4_expGolombEncode(diff_sub_mv(2))) + strlength(A1_Q4_expGolombEncode(diff_sub_ref_idx));
                    sub_J = sub_J + lambda * mv_bits;

                    J_split = J_split + sub_J;
                    sub_mvs(idx, :) = sub_mv;
                    sub_ref_idxs(idx, :) = sub_ref_idx;
                    
                    rowRange = relSubRow:relSubRow+split_size-1;
                    colRange = relSubCol:relSubCol+split_size-1;
                    sub_predicted(rowRange, colRange) = sub_pred;
                    sub_encoded_blocks{idx} = sub_encoded;

                    % Update temporary previous MV
                    temp_prev_mv = sub_mv;
                    temp_prev_ref_idx = sub_ref_idx;
                end

                % Add split flag cost
                J_split = J_split + lambda * strlength(A1_Q4_expGolombEncode(1));
                J_ns = J_ns + lambda * strlength(A1_Q4_expGolombEncode(0));
               
                % Choose between split and non-split based on RD cost
                if J_split < J_ns
                    % Use split blocks
                    fprintf(Diff_stream, '%s %s\n', A1_Q4_expGolombEncode(1), A1_Q4_expGolombEncode(0)); % Split=1, P-frame=0
                    
                    % Process each sub-block
                    for idx = 1:4
                        diff_mv = sub_mvs(idx, :) - previous_mv;
                        previous_mv = sub_mvs(idx, :);
                        relSubRow = relative_offsets(idx, 1) + 1;
                        relSubCol = relative_offsets(idx, 2) + 1;

                        diff_ref_index = sub_ref_idxs(idx, :) - previous_ref_index;
                        previous_ref_index = sub_ref_idxs(idx, :);

                        fprintf(Diff_stream, '%s %s %s\n', A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)), A1_Q4_expGolombEncode(diff_ref_index));
                        fprintf(QTC_stream, '%s\n', sub_encoded_blocks{idx});
                        
                        % Update reconstructed frame for this sub-block
                        absSubRow = row + relative_offsets(idx, 1);
                        absSubCol = col + relative_offsets(idx, 2);
                        
                        % Reconstruct sub-block
                        sub_decoded = A1_Q4_expGolombDecode(sub_encoded_blocks{idx});
                        sub_decoded = A1_Q4_rleDecode(sub_decoded, split_size);
                        sub_decoded = A1_Q4_inverseSScan(sub_decoded, split_size, split_size);
                        sub_residual = A1_Q4_idctAfterDequantizeBlock(sub_decoded, A1_Q4_generateQMatrix(split_size, QP));
                        
                        sub_recon = sub_predicted(relSubRow:relSubRow+split_size-1, relSubCol:relSubCol+split_size-1) + sub_residual;
                        sub_recon = max(min(sub_recon, 255), 0);
                        
                        reconstructedFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1) = sub_recon;
                        predictedFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1) = sub_predicted(relSubRow:relSubRow+split_size-1, relSubCol:relSubCol+split_size-1);
                    end
                else
                    % Use non-split block
                    % Update motion vector differentiation
                    diff_mv = bestMatch_ns - previous_mv;
                    previous_mv = bestMatch_ns;
                    
                    diff_ref_index = best_ref_index_ns - previous_ref_index;
                    previous_ref_index = best_ref_index_ns;
                    
                    % Write motion vector and reference index differences to stream
                    fprintf(Diff_stream, '%s %s %s %s %s\n', ...
                        A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(0), ...
                        A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)), ...
                        A1_Q4_expGolombEncode(diff_ref_index));
                    
                    % Compute residual block
                    residualBlock_ns = A1_Q3_calcResidual(predictedBlock_ns, currentBlock, n);
                    
                    % Encode residual block
                    [~, encoded_block_ns] = A2_computeRD(residualBlock_ns, 0, QP, lambda);
                    fprintf(QTC_stream, '%s\n', encoded_block_ns);
                    
                    % Decode residual block
                    decoded_block_ns = A1_Q4_expGolombDecode(encoded_block_ns);
                    decoded_block_ns = A1_Q4_rleDecode(decoded_block_ns, blockSize);
                    decoded_block_ns = A1_Q4_inverseSScan(decoded_block_ns, blockSize, blockSize);
                    decoded_residual_ns = A2_idctAfterDequantizeBlock(decoded_block_ns, Q_Matrix);
                    
                    % Reconstruct block
                    reconstructedBlock_ns = double(predictedBlock_ns) + decoded_residual_ns;
                    reconstructedBlock_ns = max(min(reconstructedBlock_ns, 255), 0);
                    
                    % Update reconstructed frame
                    reconstructedFrame(row:row+blockSize-1, col:col+blockSize-1) = reconstructedBlock_ns;
                    predictedFrame(row:row+blockSize-1, col:col+blockSize-1) = predictedBlock_ns;
                    
                    % Apply overlay if MRFoverlay is enabled
                    if MRFoverlay == 1
                        overlayColor = colors(best_ref_index_ns, :);
                        for channel = 1:3
                            overlayFrame(row:row+blockSize-1, col:col+blockSize-1, channel) = ...
                                (0.5 * reconstructedFrame(row:row+blockSize-1, col:col+blockSize-1)) + ...
                                (0.5 * overlayColor(channel)); % Alpha blending (50%)
                        end
                    end
                end
            else
                % Non-VBS path
                % find the best match mv and predicted block
                diff_ref_index = best_ref_index_ns - previous_ref_index;
                previous_ref_index = best_ref_index_ns;
                
                % Update motion vector differentiation
                diff_mv = bestMatch_ns - previous_mv;
                previous_mv = bestMatch_ns;
                
                % Write motion vector and reference index differences to the stream
                fprintf(Diff_stream, '%s %s %s %s %s\n', ...
                    A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(0), ...
                    A1_Q4_expGolombEncode(diff_mv(1)), ...
                    A1_Q4_expGolombEncode(diff_mv(2)), ...
                    A1_Q4_expGolombEncode(diff_ref_index));
                
                % Update predicted frame for the current block
                predictedFrame(row:row+blockSize-1, col:col+blockSize-1) = predictedBlock_ns;
                
                % Calculate residual block
                residualBlock = A1_Q3_calcResidual(predictedBlock_ns, currentBlock, n);
                
                % Encode residual block
                [~, encoded_block_ns] = A2_computeRD(residualBlock, 0, QP, lambda);
                fprintf(QTC_stream, '%s\n', encoded_block_ns);
                
                % Decode the residual block
                decoded_block_ns = A1_Q4_expGolombDecode(encoded_block_ns);
                decoded_block_ns = A1_Q4_rleDecode(decoded_block_ns, blockSize);
                decoded_block_ns = A1_Q4_inverseSScan(decoded_block_ns, blockSize, blockSize);
                decoded_residual_ns = A2_idctAfterDequantizeBlock(decoded_block_ns, Q_Matrix);
                
                % Reconstruct the block
                reconstructedBlock_ns = double(predictedBlock_ns) + decoded_residual_ns;
                reconstructedBlock_ns = max(min(reconstructedBlock_ns, 255), 0);
                
                % Update the reconstructed frame for the current block
                reconstructedFrame(row:row+blockSize-1, col:col+blockSize-1) = reconstructedBlock_ns;
                predictedFrame(row:row+blockSize-1, col:col+blockSize-1) = predictedBlock_ns;
                
                % Apply overlay if enabled
                if MRFoverlay == 1
                    overlayColor = colors(best_ref_index_ns, :);
                    for channel = 1:3
                        overlayFrame(row:row+blockSize-1, col:col+blockSize-1, channel) = ...
                            (0.5 * reconstructedFrame(row:row+blockSize-1, col:col+blockSize-1)) + ...
                            (0.5 * overlayColor(channel)); % Alpha blending (50%)
                    end
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
    
        %fprintf('Saved overlay frame %d as %s\n', frameIdx, outputFileName);
    end
end