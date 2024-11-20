function [predictedFrame, reconstructedFrame] = A2_Q2_interPredictForPFrame(referenceFrame, currentFrame, searchRange, blockSize, paddedHeight, paddedWidth, n, QP, MDiff_stream, QTC_stream, VBSenable)
    Q_Matrix = A2_Q2_generateQMatrix(blockSize, QP);
    previous_mv = [0, 0];
    predictedFrame = zeros(paddedHeight, paddedWidth);
    reconstructedFrame = zeros(paddedHeight, paddedWidth);
    lambda = A2_Q2_getLambda(QP);
    
    for i = 1:blockSize:paddedHeight
        for j = 1:blockSize:paddedWidth
            % Boundary check
            block_height = min(blockSize, paddedHeight - i + 1);
            block_width = min(blockSize, paddedWidth - j + 1);
            
            if block_height < blockSize || block_width < blockSize
                continue;
            end

            % Get current block
            currentBlock = currentFrame(i:i+block_height-1, j:j+block_width-1);

            if VBSenable
                % Calculate costs for non-split case
                [bestMatch_ns, predictedBlock_ns] = InterPredictBLK(referenceFrame, currentBlock, i, j, searchRange, blockSize, paddedWidth, paddedHeight);
                residualBlock_ns = A1_Q3_calcResidual(predictedBlock_ns, currentBlock, n);
                [J_ns, encoded_block_ns] = A2_Q2_computeRD(residualBlock_ns, 0, QP, lambda);
                
                % Add MV cost to non-split RD cost
                diff_mv_ns = bestMatch_ns - previous_mv;
                mv_bits_ns = strlength(A1_Q4_expGolombEncode(diff_mv_ns(1))) + strlength(A1_Q4_expGolombEncode(diff_mv_ns(2)));
                J_ns = J_ns + lambda * mv_bits_ns;

                % Calculate split case (4 sub-blocks)
                split_size = blockSize/2;
                J_split = 0;
                relative_offsets = [0, 0; 0, split_size; split_size, 0; split_size, split_size];
                sub_mvs = zeros(4, 2);
                sub_predicted = zeros(blockSize);
                sub_reconstructed = zeros(blockSize);
                sub_encoded_blocks = cell(4, 1);
                temp_prev_mv = previous_mv;

                % Process sub-blocks in Z-order
                for idx = 1:4
                    relSubRow = relative_offsets(idx, 1) + 1;
                    relSubCol = relative_offsets(idx, 2) + 1;
                    absSubRow = i + relative_offsets(idx, 1);
                    absSubCol = j + relative_offsets(idx, 2);

                    sub_block = currentFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1);
                    
                    [sub_mv, sub_pred] = InterPredictBLK(referenceFrame, sub_block, absSubRow, absSubCol, searchRange, split_size, paddedWidth, paddedHeight);
                    sub_residual = A1_Q3_calcResidual(sub_pred, sub_block, n);
                    [sub_J, sub_encoded] = A2_Q2_computeRD(sub_residual, 0, QP, lambda);

                    % Add MV cost
                    diff_sub_mv = sub_mv - temp_prev_mv;
                    mv_bits = strlength(A1_Q4_expGolombEncode(diff_sub_mv(1))) + strlength(A1_Q4_expGolombEncode(diff_sub_mv(2)));
                    sub_J = sub_J + lambda * mv_bits;

                    J_split = J_split + sub_J;
                    sub_mvs(idx, :) = sub_mv;
                    
                    rowRange = relSubRow:relSubRow+split_size-1;
                    colRange = relSubCol:relSubCol+split_size-1;
                    sub_predicted(rowRange, colRange) = sub_pred;
                    sub_encoded_blocks{idx} = sub_encoded;

                    % Update temporary previous MV
                    temp_prev_mv = sub_mv;
                end

                % Add split flag cost
                J_split = J_split + lambda * strlength(A1_Q4_expGolombEncode(1));
                J_ns = J_ns + lambda * strlength(A1_Q4_expGolombEncode(0));

                % Choose between split and non-split based on RD cost
                if J_split < J_ns
                    % Use split blocks
                    fprintf(MDiff_stream, '%s %s\n', A1_Q4_expGolombEncode(1), A1_Q4_expGolombEncode(0)); % Split=1, P-frame=0
                    
                    % Process each sub-block
                    for idx = 1:4
                        diff_mv = sub_mvs(idx, :) - previous_mv;
                        previous_mv = sub_mvs(idx, :);
                        relSubRow = relative_offsets(idx, 1) + 1;
                        relSubCol = relative_offsets(idx, 2) + 1;

                        fprintf(MDiff_stream, '%s %s\n', A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)));
                        fprintf(QTC_stream, '%s\n', sub_encoded_blocks{idx});
                        
                        % Update reconstructed frame for this sub-block
                        absSubRow = i + relative_offsets(idx, 1);
                        absSubCol = j + relative_offsets(idx, 2);
                        
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
                    % fprintf(MDiff_stream, '%s %s ', A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(0)); % Split=0, P-frame=0
                    
                    diff_mv = bestMatch_ns - previous_mv;
                    previous_mv = bestMatch_ns;
                    
                    fprintf(MDiff_stream, '%s %s %s %s\n', A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(0), ...
                    A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)));
                    fprintf(QTC_stream, '%s\n', encoded_block_ns);
                    
                    % Reconstruct block
                    decoded = A1_Q4_expGolombDecode(encoded_block_ns);
                    decoded = A1_Q4_rleDecode(decoded, blockSize);
                    decoded = A1_Q4_inverseSScan(decoded, blockSize, blockSize);
                    residual = A1_Q4_idctAfterDequantizeBlock(decoded, Q_Matrix);
                    
                    reconstructed = predictedBlock_ns + residual;
                    reconstructed = max(min(reconstructed, 255), 0);
                    
                    reconstructedFrame(i:i+block_height-1, j:j+block_width-1) = reconstructed;
                    predictedFrame(i:i+block_height-1, j:j+block_width-1) = predictedBlock_ns;
                end
            else
                % Non-VBS path
                [bestMatch, predictedBlock] = InterPredictBLK(referenceFrame, currentBlock, i, j, searchRange, blockSize, paddedWidth, paddedHeight);
                
                diff_mv = bestMatch - previous_mv;
                previous_mv = bestMatch;
                
                fprintf(MDiff_stream, '%s %s %s %s\n', A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(0), ...
                    A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)));
                
                residualBlock = A1_Q3_calcResidual(predictedBlock, currentBlock, n);
                [~, encoded_block] = A2_Q2_computeRD(residualBlock, 0, QP, lambda);
                fprintf(QTC_stream, '%s\n', encoded_block);
                
                % Reconstruct block
                decoded = A1_Q4_expGolombDecode(encoded_block);
                decoded = A1_Q4_rleDecode(decoded, blockSize);
                decoded = A1_Q4_inverseSScan(decoded, blockSize, blockSize);
                residual = A1_Q4_idctAfterDequantizeBlock(decoded, Q_Matrix);
                
                reconstructed = predictedBlock + residual;
                reconstructed = max(min(reconstructed, 255), 0);
                
                reconstructedFrame(i:i+block_height-1, j:j+block_width-1) = reconstructed;
                predictedFrame(i:i+block_height-1, j:j+block_width-1) = predictedBlock;
            end
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