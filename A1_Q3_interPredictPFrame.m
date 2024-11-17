function [predictedFrame, reconstructedFrame] = A1_Q3_interPredictPFrame(referenceFrame, currentFrame, searchRange, blockSize, paddedHeight, paddedWidth, n, QP, MDiff_stream, QTC_stream, VBSenable)
    Q_Matrix = A1_Q4_generateQMatrix(blockSize, QP);
    previous_mv = [0, 0];

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
            [bestMatch, predictedBlock] = InterPredictBLK(referenceFrame, currentBlock, row, col, searchRange, blockSize, paddedWidth, paddedHeight);

            % Update predicted frame for next iteration
            predictedFrame(row:row+blockSize-1, col:col+blockSize-1) = predictedBlock;

            % update the motion vector differentiation
            diff_mv = bestMatch - previous_mv;
            previous_mv = bestMatch;

            % generate MDiff with MV
            fprintf(MDiff_stream, '%s %s %s\n', A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)));
        
            % Save the motion vector
            % fprintf(mvFile, 'Frame %d, Block (%d, %d): MV = (%d, %d)\n', frameIdx, row, col, bestMatch(1), bestMatch(2));
            
            % find Residual block
            approxResidualBlock = A1_Q3_calcResidual(predictedBlock, currentBlock, n);

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
            reconstructedBlock = double(predictedBlock) + decoded_residual_block;
            %reconstructedBlock = max(min(reconstructedBlock, 255), 0);
            reconstructedFrame(row:row+blockSize-1, col:col+blockSize-1) = reconstructedBlock;
        end
    end
end

function [bestMatch, predictedBlock] = InterPredictBLK(referenceFrame, currentBlock, row, col, searchRange, blockSize, paddedWidth, paddedHeight)
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
            mae = mean(abs(double(currentBlock(:)) - double(referenceBlock(:))));

            % Update best match based on MAE and motion vector criteria
            if mae < bestMAE || ...
               (mae == bestMAE && abs(xOffset) + abs(yOffset) < abs(bestMatch(1)) + abs(bestMatch(2))) || ...
               (mae == bestMAE && abs(xOffset) + abs(yOffset) == abs(bestMatch(1)) + abs(bestMatch(2)) && ...
               (yOffset < bestMatch(2) || (yOffset == bestMatch(2) && xOffset < bestMatch(1))))
                bestMAE = mae;
                bestMatch = [xOffset, yOffset];
            end
        end
    end

    % Predicted block based on best match
    bestRow = row + bestMatch(2);
    bestCol = col + bestMatch(1);
    predictedBlock = referenceFrame(bestRow:bestRow+blockSize-1, bestCol:bestCol+blockSize-1);
end