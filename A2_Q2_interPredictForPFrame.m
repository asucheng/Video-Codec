function [predictedFrame, reconstructedFrame] = A2_Q2_interPredictForPFrame(referenceFrame, currentFrame, blockSize, searchRange, paddedHeight, paddedWidth, QP, MDiff_stream, QTC_stream, VBSenable)
    Q_Matrix = A2_Q2_generateQMatrix(blockSize, QP);
    previous_mv = [0, 0];  % Previous motion vector for differential encoding
    lambda = A2_Q2_getLambda(QP);  % RD cost parameter

    predictedFrame = zeros(paddedHeight, paddedWidth);
    reconstructedFrame = zeros(paddedHeight, paddedWidth);

    for row = 1:blockSize:paddedHeight
        for col = 1:blockSize:paddedWidth
            block_height = min(blockSize, paddedHeight - row + 1);
            block_width = min(blockSize, paddedWidth - col + 1);

            % Extract the current block
            currentBlock = currentFrame(row:row+block_height-1, col:col+block_width-1);

            % Perform motion estimation
            [bestMatch, predictedBlock] = A2_Q2_interPredictBlock(referenceFrame, currentBlock, row, col, searchRange, blockSize, paddedWidth, paddedHeight);

            % Compute residual for the non-split block
            residual_block = currentBlock - predictedBlock;
            [J_non_split, ~] = A2_Q2_computeRD(residual_block, 0, QP, lambda); % Mode = 0 for simplicity

            % Handle Variable Block Size (VBS)
            split_flag = false;
            if VBSenable && blockSize > 4
                % Divide the block into four sub-blocks
                split_size = floor(blockSize / 2);
                sub_blocks = {
                    currentBlock(1:split_size, 1:split_size),         % Top-left
                    currentBlock(1:split_size, split_size+1:end),    % Top-right
                    currentBlock(split_size+1:end, 1:split_size),    % Bottom-left
                    currentBlock(split_size+1:end, split_size+1:end) % Bottom-right
                };

                % Initialize RD cost for split
                J_split = 0;

                % Process each sub-block
                for k = 1:4
                    subRow = row + (k > 2) * split_size;
                    subCol = col + (mod(k - 1, 2) == 1) * split_size;

                    % Predict sub-block using motion estimation
                    subBlock = sub_blocks{k};
                    [bestSubMatch, predictedSubBlock] = A2_Q2_interPredictBlock(referenceFrame, subBlock, subRow, subCol, searchRange, split_size, paddedWidth, paddedHeight);

                    % Compute residuals and RD cost for each sub-block
                    residual_sub_block = subBlock - predictedSubBlock;
                    [J_sub, ~] = A2_Q2_computeRD(residual_sub_block, 0, QP - 1, lambda); % Mode = 0 for sub-blocks
                    J_split = J_split + J_sub;
                end

                % Compare RD costs
                if J_split < J_non_split
                    split_flag = true;
                end
            end

            % Encode and reconstruct based on split decision
            if split_flag
                % Encode split flag
                fprintf(MDiff_stream, '%s\n', A1_Q4_expGolombEncode(1)); % Split flag = 1

                % Process each sub-block
                reconstructedBlock = zeros(blockSize, blockSize);
                split_size = floor(blockSize / 2);
                offsets = [0, 0; 0, split_size; split_size, 0; split_size, split_size];
                for k = 1:4
                    subRow = row + offsets(k, 1);
                    subCol = col + offsets(k, 2);
                    subBlock = currentFrame(subRow:subRow+split_size-1, subCol:subCol+split_size-1);

                    % Predict and encode each sub-block
                    [subBestMatch, predictedSubBlock] = A2_Q2_interPredictBlock(referenceFrame, subBlock, subRow, subCol, searchRange, split_size, paddedWidth, paddedHeight);

                    % Encode motion vector for sub-block
                    diff_sub_mv = subBestMatch - previous_mv;
                    previous_mv = subBestMatch;
                    fprintf(MDiff_stream, '%s %s\n', A1_Q4_expGolombEncode(diff_sub_mv(1)), A1_Q4_expGolombEncode(diff_sub_mv(2)));

                    % Compute residual and encode
                    sub_residual = subBlock - predictedSubBlock;
                    [sub_encoded, sub_quantized] = A2_Q2_quantizeAndEncode(sub_residual, QP - 1);
                    fprintf(QTC_stream, '%s\n', sub_encoded);

                    % Decode and reconstruct sub-block
                    sub_residual_decoded = A2_Q2_idctAfterDequantizeBlock(sub_quantized, A2_Q2_generateQMatrix(split_size, QP - 1));
                    reconstructedSubBlock = predictedSubBlock + sub_residual_decoded;
                    reconstructedSubBlock = max(min(reconstructedSubBlock, 255), 0);

                    % Place sub-block into parent block
                    reconstructedBlock(offsets(k, 1)+1:offsets(k, 1)+split_size, offsets(k, 2)+1:offsets(k, 2)+split_size) = reconstructedSubBlock;
                end
            else
                % Encode non-split flag

                % Encode motion vector
                diff_mv = bestMatch - previous_mv;
                previous_mv = bestMatch;
                fprintf(MDiff_stream, '%s %s %s\n', A1_Q4_expGolombEncode(0), ...
                A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)));
                % Encode residuals for non-split block
                [encoded_block, quantized_block] = A2_Q2_quantizeAndEncode(residual_block, QP);
                fprintf(QTC_stream, '%s\n', encoded_block);

                % Decode and reconstruct non-split block
                decoded_residual_block = A2_Q2_idctAfterDequantizeBlock(quantized_block, Q_Matrix);
                reconstructedBlock = predictedBlock + decoded_residual_block;
                reconstructedBlock = max(min(reconstructedBlock, 255), 0);
            end

            % Update predicted and reconstructed frames
            predictedFrame(row:row+block_height-1, col:col+block_width-1) = predictedBlock;
            reconstructedFrame(row:row+block_height-1, col:col+block_width-1) = reconstructedBlock;
        end
    end
end