function [predictedFrame, reconstructedFrame] = A2_Q34_interPredictForPFrame(referenceFrame, currentFrame, searchRange, blockSize, paddedHeight, paddedWidth, n, QP, MDiff_stream, MVPDiff_stream, QTC_stream, VBSenable, FMEEnable, FastME)
    % Initialize parameters
    Q_Matrix = A1_Q4_generateQMatrix(blockSize, QP);
    Lambda = A2_Q2_getLambda(QP);
    previous_mv = [0, 0];
    previous_mvp = [0, 0];
    mvp = [0, 0];
    
    predictedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
    reconstructedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
    
    % Loop over blocks
    for i = 1:blockSize:paddedHeight
        for j = 1:blockSize:paddedWidth
            % Adjust block dimensions for edge cases
            block_height = min(blockSize, paddedHeight - i + 1);
            block_width = min(blockSize, paddedWidth - j + 1);
            if block_height < blockSize || block_width < blockSize
                continue;
            end
            currentBlock = currentFrame(i:i+block_height-1, j:j+block_width-1);
            
            if VBSenable && blockSize >= 8  % Minimum block size for splitting
                % **Variable Block Size Logic**
                
                % Compute RD cost for non-split case
                [best_mv_ns, predictedBlock_ns] = InterPredictBLK(referenceFrame, currentBlock, i, j, searchRange, blockSize, paddedWidth, paddedHeight, FMEEnable, FastME, mvp);
                residualBlock_ns = A1_Q3_calcResidual(predictedBlock_ns, currentBlock, n);
                mode = 0;  % Mode can be set to 0 for inter-prediction
                [J_ns, encoded_block_ns] = A2_Q2_computeRD(residualBlock_ns, mode, QP, Lambda);
                
                % Add MV cost to non-split RD cost
                diff_mv_ns = best_mv_ns - previous_mv;
                mv_bits_ns = A1_Q4_bitcountFromArray([diff_mv_ns(1), diff_mv_ns(2)]);
                J_ns = J_ns + Lambda * mv_bits_ns;
                
                % **Split Block Logic**
                split_size = blockSize / 2;
                J_split = 0;
                sub_mvs = zeros(4, 2);
                sub_encoded_blocks = cell(4, 1);
                temp_prev_mv = previous_mv;
                relative_offsets = [0, 0; 0, split_size; split_size, 0; split_size, split_size];
                sub_predicted = zeros(blockSize, blockSize, 'uint8');
                
                for idx = 1:4
                    relSubRow = relative_offsets(idx, 1) + 1;
                    relSubCol = relative_offsets(idx, 2) + 1;
                    absSubRow = i + relative_offsets(idx, 1);
                    absSubCol = j + relative_offsets(idx, 2);
                    
                    sub_block_height = min(split_size, paddedHeight - absSubRow + 1);
                    sub_block_width = min(split_size, paddedWidth - absSubCol + 1);
                    if sub_block_height < split_size || sub_block_width < split_size
                        continue;
                    end
                    
                    sub_block = currentFrame(absSubRow:absSubRow+sub_block_height-1, absSubCol:absSubCol+sub_block_width-1);
                    
                    % Motion estimation for sub-block
                    [sub_mv, sub_pred] = InterPredictBLK(referenceFrame, sub_block, absSubRow, absSubCol, searchRange, split_size, paddedWidth, paddedHeight, FMEEnable, FastME, mvp);
                    
                    sub_residual = A1_Q3_calcResidual(sub_pred, sub_block, n);
                    mode = 0;  % Mode for inter-prediction
                    [sub_J, sub_encoded] = A2_Q2_computeRD(sub_residual, mode, QP, Lambda);
                    
                    % Add MV cost
                    diff_sub_mv = sub_mv - temp_prev_mv;
                    mv_bits = A1_Q4_bitcountFromArray([diff_sub_mv(1), diff_sub_mv(2)]);
                    sub_J = sub_J + Lambda * mv_bits;
                    
                    J_split = J_split + sub_J;
                    sub_mvs(idx, :) = sub_mv;
                    temp_prev_mv = sub_mv;
                    
                    % Store predicted block and encoded data
                    rowRange = relSubRow:relSubRow+sub_block_height-1;
                    colRange = relSubCol:relSubCol+sub_block_width-1;
                    sub_predicted(rowRange, colRange) = sub_pred;
                    sub_encoded_blocks{idx} = sub_encoded;
                end
                
                % Add split flag cost
                J_split = J_split + Lambda * A1_Q4_bitcountFromArray(1);
                J_ns = J_ns + Lambda * A1_Q4_bitcountFromArray(0);
                
                % **Choose between split and non-split**
                if J_split < J_ns
                    % Use split blocks
                    fprintf(MDiff_stream, '%s %s \n', A1_Q4_expGolombEncode(1),A1_Q4_expGolombEncode(0));  % Split flag = 1
                    % Process each sub-block
                    for idx = 1:4
                        diff_mv = sub_mvs(idx, :) - previous_mv;
                        previous_mv = sub_mvs(idx, :);
                        fprintf(MDiff_stream, '%s %s\n', A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)));
                        fprintf(QTC_stream, '%s\n', sub_encoded_blocks{idx});
                        
                        % Reconstruct sub-block
                        relSubRow = relative_offsets(idx, 1) + 1;
                        relSubCol = relative_offsets(idx, 2) + 1;
                        absSubRow = i + relative_offsets(idx, 1);
                        absSubCol = j + relative_offsets(idx, 2);
                        sub_block_height = min(split_size, paddedHeight - absSubRow + 1);
                        sub_block_width = min(split_size, paddedWidth - absSubCol + 1);
                        rowRange = absSubRow:absSubRow+sub_block_height-1;
                        colRange = absSubCol:absSubCol+sub_block_width-1;
                        
                        % Decoding
                        decoded_block = A1_Q4_expGolombDecode(sub_encoded_blocks{idx});
                        decoded_coeffs = A1_Q4_rleDecode(decoded_block, split_size);
                        decoded_quantized_block = A1_Q4_inverseSScan(decoded_coeffs, split_size, split_size);
                        decoded_residual_block = A1_Q4_idctAfterDequantizeBlock(decoded_quantized_block, Q_Matrix(1:split_size, 1:split_size));
                        
                        reconstructedBlock = double(sub_predicted(relSubRow:relSubRow+sub_block_height-1, relSubCol:relSubCol+sub_block_width-1)) + decoded_residual_block;
                        reconstructedBlock = max(min(reconstructedBlock, 255), 0);
                        reconstructedFrame(rowRange, colRange) = uint8(reconstructedBlock);
                        predictedFrame(rowRange, colRange) = sub_predicted(relSubRow:relSubRow+sub_block_height-1, relSubCol:relSubCol+sub_block_width-1);
                    end
                else
                    % Use non-split block
                    diff_mv = best_mv_ns - previous_mv;
                    previous_mv = best_mv_ns;
                    fprintf(MDiff_stream, '%s %s %s\n', A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)));
                    fprintf(QTC_stream, '%s\n', encoded_block_ns);
                    
                    % Reconstruct block
                    decoded_block = A1_Q4_expGolombDecode(encoded_block_ns);
                    decoded_coeffs = A1_Q4_rleDecode(decoded_block, blockSize);
                    decoded_quantized_block = A1_Q4_inverseSScan(decoded_coeffs, blockSize, blockSize);
                    decoded_residual_block = A1_Q4_idctAfterDequantizeBlock(decoded_quantized_block, Q_Matrix);
                    
                    reconstructedBlock = double(predictedBlock_ns) + decoded_residual_block;
                    reconstructedBlock = max(min(reconstructedBlock, 255), 0);
                    reconstructedFrame(i:i+block_height-1, j:j+block_width-1) = uint8(reconstructedBlock);
                    predictedFrame(i:i+block_height-1, j:j+block_width-1) = predictedBlock_ns;
                end
            else
                % **Non-VBS Path**
                [best_mv, predictedBlock] = InterPredictBLK(referenceFrame, currentBlock, i, j, searchRange, blockSize, paddedWidth, paddedHeight, FMEEnable, FastME, mvp);
                
                % Update motion vectors and bitstreams
                if FastME
                    mvp = best_mv;
                    mvp_diff = mvp - previous_mvp;
                    previous_mvp = mvp;
                    % split flag, I-frame-flage
                    fprintf(MVPDiff_stream, '%s %s %s %s\n', A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(mvp_diff(1)), A1_Q4_expGolombEncode(mvp_diff(2)));
                else
                    diff_mv = best_mv - previous_mv;
                    previous_mv = best_mv;
                    fprintf(MDiff_stream, '%s %s %s %s\n',A1_Q4_expGolombEncode(0),A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(diff_mv(1)), A1_Q4_expGolombEncode(diff_mv(2)));
                end
                
                % Compute residual
                residualBlock = A1_Q3_calcResidual(predictedBlock, currentBlock, n);
                mode = 0;  % Mode for inter-prediction
                [~, encoded_block] = A2_Q2_computeRD(residualBlock, mode, QP, Lambda);
                fprintf(QTC_stream, '%s\n', encoded_block);
                
                % Reconstruct block
                decoded_block = A1_Q4_expGolombDecode(encoded_block);
                decoded_coeffs = A1_Q4_rleDecode(decoded_block, blockSize);
                decoded_quantized_block = A1_Q4_inverseSScan(decoded_coeffs, blockSize, blockSize);
                decoded_residual_block = A1_Q4_idctAfterDequantizeBlock(decoded_quantized_block, Q_Matrix);
                
                reconstructedBlock = double(predictedBlock) + decoded_residual_block;
                reconstructedBlock = max(min(reconstructedBlock, 255), 0);
                reconstructedFrame(i:i+block_height-1, j:j+block_width-1) = uint8(reconstructedBlock);
                predictedFrame(i:i+block_height-1, j:j+block_width-1) = predictedBlock;
            end
        end
    end
end

function [best_mv, predictedBlock] = InterPredictBLK(referenceFrame, currentBlock, row, col, searchRange, blockSize, paddedWidth, paddedHeight, FMEEnable, FastME, mvp)
    bestMAE = Inf;
    best_mv = [0, 0];
    block_height = size(currentBlock, 1);
    block_width = size(currentBlock, 2);
    predictedBlock = zeros(block_height, block_width, 'uint8');
    
    if FastME
        % Fast Motion Estimation using MVP
        candidates = [
            mvp;
            mvp + [0, 1];
            mvp + [0, -1];
            mvp + [1, 0];
            mvp + [-1, 0];
        ];
    else
        % Full search
        [dx, dy] = meshgrid(-searchRange:searchRange, -searchRange:searchRange);
        candidates = [dx(:), dy(:)];
    end
    
    for idx = 1:size(candidates, 1)
        mv = candidates(idx, :);
        refRow = row + mv(2);
        refCol = col + mv(1);
        
        % Boundary check
        if refRow < 1 || refCol < 1 || refRow + block_height - 1 > paddedHeight || refCol + block_width -1 > paddedWidth
            continue;
        end
        
        % Get reference block
        if FMEEnable
            referenceBlock = PerformInterpolation(referenceFrame, mv, refRow, refCol, block_height, block_width);
        else
            referenceBlock = referenceFrame(refRow:refRow+block_height-1, refCol:refCol+block_width-1);
        end
        
        % Compute MAE
        mae = mean(abs(double(currentBlock(:)) - double(referenceBlock(:))));
        
        % Update best match
        if mae < bestMAE || ...
           (mae == bestMAE && sum(abs(mv)) < sum(abs(best_mv))) || ...
           (mae == bestMAE && sum(abs(mv)) == sum(abs(best_mv)) && (mv(2) < best_mv(2) || (mv(2) == best_mv(2) && mv(1) < best_mv(1))))
            bestMAE = mae;
            best_mv = mv;
            predictedBlock = referenceBlock;
        end
    end
end

function interpolatedBlock = PerformInterpolation(refFrame, mv, refRow, refCol, block_height, block_width)
    int_mv = floor(mv);
    frac_mv = mv - int_mv;
    
    y = refRow;
    x = refCol;
    [ref_height, ref_width] = size(refFrame);
    
    interpolatedBlock = zeros(block_height, block_width, 'uint8');
    
    for r = 1:block_height
        for c = 1:block_width
            y_pos = y + r - 1 + frac_mv(2);
            x_pos = x + c - 1 + frac_mv(1);
            
            % Bilinear interpolation
            interpolatedBlock(r, c) = bilinearInterpolation(refFrame, y_pos, x_pos, ref_height, ref_width);
        end
    end
end

function value = bilinearInterpolation(refFrame, y, x, ref_height, ref_width)
    x1 = floor(x);
    x2 = ceil(x);
    y1 = floor(y);
    y2 = ceil(y);
    
    x1 = max(1, min(ref_width, x1));
    x2 = max(1, min(ref_width, x2));
    y1 = max(1, min(ref_height, y1));
    y2 = max(1, min(ref_height, y2));
    
    Q11 = double(refFrame(y1, x1));
    Q21 = double(refFrame(y1, x2));
    Q12 = double(refFrame(y2, x1));
    Q22 = double(refFrame(y2, x2));
    
    value = uint8(...
        Q11 * (x2 - x) * (y2 - y) + ...
        Q21 * (x - x1) * (y2 - y) + ...
        Q12 * (x2 - x) * (y - y1) + ...
        Q22 * (x - x1) * (y - y1));
end