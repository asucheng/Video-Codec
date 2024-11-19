function [predictedFrame, reconstructedFrame] = A2_Q34_interPredictForPFrame(referenceFrame, currentFrame, searchRange, ...
    blockSize, paddedHeight, paddedWidth, n, QP, ...
    MDiff_stream, MVPDiff_stream, QTC_stream, FMEEnable, FastME)

    Q_Matrix = A2_Q34_generateQMatrix(blockSize, QP);
    previous_mv = [0, 0];
    previous_mvp = [0, 0];
    mvp = [0, 0];

    predictedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
    reconstructedFrame = zeros(paddedHeight, paddedWidth, 'uint8');

    % Loop over blocks
    for row = 1:blockSize:paddedHeight
        for col = 1:blockSize:paddedWidth
            % Boundary check
            if row + blockSize - 1 > paddedHeight || col + blockSize - 1 > paddedWidth
                continue;
            end

            % Current block
            currentBlock = currentFrame(row:row+blockSize-1, col:col+blockSize-1);
            
            % find the best match mv and predicted block
            [best_mv, predictedBlock] = InterPredictBLK(referenceFrame, currentBlock, row, col, ...
                searchRange, blockSize, paddedWidth, paddedHeight, ...
                FMEEnable, FastME, mvp);

            % Update predicted frame for next iteration
            predictedFrame(row:row+blockSize-1, col:col+blockSize-1) = predictedBlock;

            % update the motion vector differentiation
            if FastME
                mvp = best_mv;
                mvp_diff = mvp - previous_mvp;
                previous_mvp = mvp;
                fprintf(MVPDiff_stream, '%s %s %s\n', A1_Q4_expGolombEncode(0), ...
                    A1_Q4_expGolombEncode(mvp_diff(1)), ...
                    A1_Q4_expGolombEncode(mvp_diff(2)));
            else
                diff_mv = best_mv - previous_mv;
                previous_mv = best_mv;
                % generate MDiff with MV
                fprintf(MDiff_stream, '%s %s %s\n', A1_Q4_expGolombEncode(0), ...
                    A1_Q4_expGolombEncode(diff_mv(1)), ...
                    A1_Q4_expGolombEncode(diff_mv(2)));
            end
            % Save the motion vector
            % fprintf(mvFile, 'Frame %d, Block (%d, %d): MV = (%d, %d)\n', frameIdx, row, col, best_mv(1), best_mv(2));
            
            % find Residual block
            residualBlock = uint16(currentBlock) - uint16(predictedBlock);

            % encode the residual block
            encoded_residual_block = A2_Q34_quantizeBlockAfterDCT(residualBlock, Q_Matrix);
            % generate QTC stream
            scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
            rle_encoded = A1_Q4_rleEncode(scanned_coeffs, blockSize);
            % QTC_stream = [QTC_stream, A1_Q4_expGolombEncode(rle_encoded)];
            encoded_value = A1_Q4_expGolombEncode(rle_encoded);
            fprintf(QTC_stream, '%s\n', encoded_value);   % each line is a block

            decoded_residual_block = A2_Q34_idctAfterDequantizeBlock(encoded_residual_block, Q_Matrix);

            % Reconstruct block and update the reconstructed frame
            reconstructedBlock = double(predictedBlock) + decoded_residual_block;
            reconstructedBlock = max(min(reconstructedBlock, 255), 0);
            reconstructedFrame(row:row+blockSize-1, col:col+blockSize-1) = reconstructedBlock;
        end
    end
end

function [bestMatch, bestPredictedBlock] = InterPredictBLK(referenceFrame, currentBlock, row, col, ...
    searchRange, blockSize, paddedWidth, paddedHeight, ...
    FMEEnable, FastME, mvp)

    bestMAE = inf;
    bestMatch = [0, 0];
    bestPredictedBlock = zeros(blockSize, blockSize, 'uint8');

    if FastME
        candidates = [
            mvp,
            mvp + [0, 1],
            mvp + [0, -1],
            mvp + [1, 0],
            mvp + [-1, 0]
        ];
    else
        [dy, dx] = meshgrid(-searchRange: searchRange, -searchRange: searchRange);
        candidates = [dy(:), dx(:)];
    end

    % Full search within search range
    for i = 1:size(candidates, 1)
        mv = candidates(i, :);
        refRow = row + mv(1);
        refCol = col + mv(2);

        if refRow < 1 || refCol < 1 || refRow+blockSize-1 > paddedHeight || refCol+blockSize-1 > paddedWidth
            continue;
        end

        % Reference block
        referenceBlock = referenceFrame(refRow:refRow+blockSize-1, refCol:refCol+blockSize-1);

        if FMEEnable
            referenceBlock = PerformInterpolation(referenceFrame, mv, refRow, refCol, blockSize);
        end
        mae = mean(abs(double(currentBlock(:)) - double(referenceBlock(:))));

        % Update best match based on MAE and motion vector criteria
        if mae < bestMAE || ...
           (mae == bestMAE && abs(mv(2)) + abs(mv(1)) < abs(bestMatch(1)) + abs(bestMatch(2))) || ...
           (mae == bestMAE && abs(mv(2)) + abs(mv(1)) == abs(bestMatch(1)) + abs(bestMatch(2)) && ...
           (mv(1) < bestMatch(2) || (mv(1) == bestMatch(2) && mv(2) < bestMatch(1))))
            bestMAE = mae;
            bestMatch = mv;
            bestPredictedBlock = referenceBlock;
        end
    end
end


function [interpolatedBlock] = PerformInterpolation(refFrame, mv, row, col, blockSize)
    int_dy = floor(mv(1) / 2); % Integer vertical motion
    int_dx = floor(mv(2) / 2); % Integer horizontal motion
    frac_dy = mod(mv(1), 2) / 2; % Fractional vertical motion (0 or 0.5)
    frac_dx = mod(mv(2), 2) / 2; % Fractional horizontal motion (0 or 0.5)

    % Perform interpolation for fractional motion estimation.
    interpolatedBlock = zeros(blockSize, blockSize, 'uint8');

    [ref_height, ref_width] = size(refFrame);

    row = row + int_dy;
    col = col + int_dx;

    for r = 1:blockSize
        for c = 1:blockSize
            % Get fractional pixel location
            y = row + r - 1;
            x = col + c - 1;

            top_left_y = max(1, min(y, ref_height));
            top_left_x = max(1, min(x, ref_width));
            top_right_x = max(1, min(x + 1, ref_width));
            bottom_left_y = max(1, min(y + 1, ref_height));
            bottom_right_y = max(1, min(y + 1, ref_height));
            bottom_right_x = max(1, min(x + 1, ref_width));
    
            top_left = refFrame(top_left_y, top_left_x);
            top_right = refFrame(top_left_y, top_right_x);
            bottom_left = refFrame(bottom_left_y, top_left_x);
            bottom_right = refFrame(bottom_right_y, bottom_right_x);

            % Direct reference pixel (MV = [0, 0])
            if frac_dx == 0 && frac_dy == 0
                interpolatedBlock(r, c) = top_left;

            % Horizontal interpolation (MV = [0, 1])
            elseif frac_dy == 0 && frac_dx ~= 0
                interpolatedBlock(r, c) = uint8(top_left * (1 - frac_dx) + top_right * frac_dx);

            % Vertical interpolation (MV = [1, 0])
            elseif frac_dy ~= 0 && frac_dx == 0
                interpolatedBlock(r, c) = uint8(top_left * (1 - frac_dy) + bottom_left * frac_dy);

            % Diagonal interpolation (MV = [1, 1])
            else
                top_interp = top_left * (1 - frac_dx) + top_right * frac_dx;
                bottom_interp = bottom_left * (1 - frac_dx) + bottom_right * frac_dx;
                interpolatedBlock(r, c) = uint8(top_interp * (1 - frac_dy) + bottom_interp * frac_dy);
            end
        end
    end

end
