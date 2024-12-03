function [predictedFrame, reconstructedFrame, bitConsumedThisFrame, totalRowCounts] = A3_intraPredictForIFrame( ...
    original_frame, block_size, QP, ...
    MDiff_stream, MVPDiff_stream, QTC_stream, ...
    VBSEnable, ~, FastME, RCflag, bitBudgetPerFrame, iFrameBudgetTable, ...
    blockSizeCollection, QPValueCollection, QPDiff_stream)

    Q_Matrix = A3_generateQMatrix(block_size, QP);
    [h, w] = size(original_frame);
    reconstructedFrame = zeros(h, w);
    encoded_frame = zeros(h, w);
    lambda = A3_getLambda(QP);

    previous_mode = 0;
    previous_qp = 0;

    if FastME
        Diff_stream = MVPDiff_stream;
    else
        Diff_stream = MDiff_stream;
    end

    bitConsumedThisFrame = 0;
    totalRowCounts = 0;

    remainingRowBlockCounts = ceil(h / block_size);
    remainingBitForThisFrame = bitBudgetPerFrame;

    if RCflag
        
        for i = 1:block_size:size(original_frame, 1)
            bitConsumedThisRowBlock = 0;
            totalRowCounts = totalRowCounts + 1;

            bitPerRemainingRowBlock = A3_getRowBitBudget(RCflag, remainingBitForThisFrame, remainingRowBlockCounts);
            QPValue = A3_findQPForBudget(bitPerRemainingRowBlock, block_size, blockSizeCollection, QPValueCollection, iFrameBudgetTable);
            Q_Matrix = A3_generateQMatrix(block_size, QP);

            qp_diff = QPValue - previous_qp;
            previous_qp = QPValue;
            fprintf(QPDiff_stream, '%s\n', A1_Q4_expGolombEncode(qp_diff));
            bitConsumedThisFrame = bitConsumedThisFrame + A3_calcBitCounts(A1_Q4_expGolombEncode(qp_diff));

            for j = 1:block_size:size(original_frame, 2)
                % get the real size of current block
                block_height = min(block_size, h - i + 1);
                block_width = min(block_size, w - j + 1);
                split_flag = false;
       
                % get the current block
                block = original_frame(i:i+block_height-1, j:j+block_width-1);
    
                % Perform intra prediction of current block
                [predicted_block, mode, split_flag] = A3_intraPredictForIBlock(original_frame, i, j, ...
                    block_height, block_width, VBSEnable, QPValue, lambda);
    
                if VBSEnable && split_flag
                     % Encode split flag, I_frame_flag
                     fprintf(Diff_stream, '%s %s\n', A1_Q4_expGolombEncode(1), A1_Q4_expGolombEncode(1)); % 1 indicates split
                     bitConsumedThisRowBlock = bitConsumedThisRowBlock + A3_calcBitCounts(A1_Q4_expGolombEncode(1));
                     bitConsumedThisRowBlock = bitConsumedThisRowBlock + A3_calcBitCounts(A1_Q4_expGolombEncode(1));
    
                     % Process each sub-block in Z-order
                     split_size = floor(block_size / 2);
                     relative_offsets = [0, 0; 0, split_size; split_size, 0; split_size, split_size];
                     sub_Q_Matrix = A3_generateQMatrix(split_size, QPValue);
    
                    for idx = 1:4
                        % Calculate relative and absolute positions for the sub-block
                        relSubRow = relative_offsets(idx, 1) + 1; % MATLAB indexing starts at 1
                        relSubCol = relative_offsets(idx, 2) + 1;
                        absSubRow = i + relative_offsets(idx, 1);
                        absSubCol = j + relative_offsets(idx, 2);
                
                        % Get the current sub-block
                        sub_block = original_frame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1);
                
                        % Perform intra prediction for sub-block
                        [sub_predicted, sub_mode, ~] = A3_intraPredictForIBlock(original_frame, absSubRow, absSubCol, ...
                            split_size, split_size, false, QPValue, lambda);
                
                        % Encode sub-block mode
                        diff_sub_mode = sub_mode - previous_mode;
                        previous_mode = sub_mode;
                        fprintf(Diff_stream, '%s\n', A1_Q4_expGolombEncode(diff_sub_mode));
                        bitConsumedThisRowBlock = bitConsumedThisRowBlock + A3_calcBitCounts(A1_Q4_expGolombEncode(diff_sub_mode));
    
                        % Calculate residual block
                        residual_block = sub_block - sub_predicted;
                
                        % Encode residual block
                        encoded_residual_block = A1_Q4_quantizeBlockAfterDCT(residual_block, sub_Q_Matrix);
                
                        % Perform scanning and RLE
                        scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
                        rle_encoded = A1_Q4_rleEncode(scanned_coeffs, split_size);
                
                        % Apply exp-Golomb encoding and write to QTC stream
                        encoded_value = A1_Q4_expGolombEncode(rle_encoded);
                        
                        [encoded_value, ~] = A3_quantizeAndEncode(residual_block, QPValue, false);
                        fprintf(QTC_stream, '%s\n', encoded_value);
                        bitConsumedThisRowBlock = bitConsumedThisRowBlock + A3_calcBitCounts(encoded_value);
    
                        % Reconstruct residual block (aligning with decoder)
                        decoded_rle = A1_Q4_expGolombDecode(encoded_value);  % Decode exp-Golomb
                        decoded_scanned = A1_Q4_rleDecode(decoded_rle, split_size);  % RLE decode
                        decoded_coeffs = A1_Q4_inverseSScan(decoded_scanned, split_size, split_size);  % Inverse scan
                        reconstructed_residual_block = A3_idctAfterDequantizeBlock(decoded_coeffs, sub_Q_Matrix);
                
                        % Reconstruct sub-block
                        sub_reconstructed = sub_predicted + reconstructed_residual_block;
                        sub_reconstructed = max(min(sub_reconstructed, 255), 0);
                
                        % Update reconstructed and predicted frames
                        reconstructedFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1) = sub_reconstructed;
                        predictedFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1) = sub_predicted;
                    end 
                else
                    predictedFrame(i:i+block_height-1, j:j+block_width-1) = predicted_block;
    
                    % encode the mode diff
                    differential_mode = mode - previous_mode;
                    previous_mode = mode;
                    encoded_diff = A1_Q4_expGolombEncode(differential_mode);
                    % split_flag, I_frame_flag, diff
                    fprintf(Diff_stream, '%s %s %s\n', A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(1), encoded_diff);
                    bitConsumedThisRowBlock = bitConsumedThisRowBlock + A3_calcBitCounts(A1_Q4_expGolombEncode(0));
                    bitConsumedThisRowBlock = bitConsumedThisRowBlock + A3_calcBitCounts(A1_Q4_expGolombEncode(1));
                    bitConsumedThisRowBlock = bitConsumedThisRowBlock + A3_calcBitCounts(encoded_diff);
             
                    % get the residuals of block
                    residual_block = block - predicted_block;
    
                    % encode the residual block
                    encoded_residual_block = A3_quantizeBlockAfterDCT(residual_block, Q_Matrix);
                    scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
                    rle_encoded = A1_Q4_rleEncode(scanned_coeffs, block_size);
                    encoded_value = A1_Q4_expGolombEncode(rle_encoded);
                    fprintf(QTC_stream, '%s\n', encoded_value);
                    bitConsumedThisRowBlock = bitConsumedThisRowBlock + A3_calcBitCounts(encoded_value);
    
                    % decode the encoded block and add it to prediction
                    decoded_residual_block = A3_idctAfterDequantizeBlock(encoded_residual_block, Q_Matrix);
                    reconstructed_block = predicted_block + decoded_residual_block;
                    reconstructed_block = max(min(reconstructed_block, 255), 0);
    
                    % store the encode residual block and reconstructed block
                    encoded_frame(i:i+block_height-1, j:j+block_width-1) = encoded_residual_block;
                    reconstructedFrame(i:i+block_height-1, j:j+block_width-1) = reconstructed_block;
                end
            end

            remainingRowBlockCounts = remainingRowBlockCounts - 1;
            remainingBitForThisFrame = remainingBitForThisFrame - bitConsumedThisRowBlock;
        end

        bitConsumedThisFrame = bitBudgetPerFrame - remainingBitForThisFrame;

    else
        % Simulate processing the frame in blocks
        for i = 1:block_size:size(original_frame, 1)
            totalRowCounts = totalRowCounts + 1;
 
            for j = 1:block_size:size(original_frame, 2)
                % get the real size of current block
                block_height = min(block_size, h - i + 1);
                block_width = min(block_size, w - j + 1);
                split_flag = false;
       
                % get the current block
                block = original_frame(i:i+block_height-1, j:j+block_width-1);
    
                % Perform intra prediction of current block
                [predicted_block, mode, split_flag] = A3_intraPredictForIBlock(original_frame, i, j, ...
                    block_height, block_width, VBSEnable, QP, lambda);
    
                if VBSEnable && split_flag
                     % Encode split flag, I_frame_flag
                     fprintf(Diff_stream, '%s %s\n', A1_Q4_expGolombEncode(1), A1_Q4_expGolombEncode(1)); % 1 indicates split
                     bitConsumedThisFrame = bitConsumedThisFrame + A3_calcBitCounts(A1_Q4_expGolombEncode(1));
                     bitConsumedThisFrame = bitConsumedThisFrame + A3_calcBitCounts(A1_Q4_expGolombEncode(1));
    
                     % Process each sub-block in Z-order
                     split_size = floor(block_size / 2);
                     relative_offsets = [0, 0; 0, split_size; split_size, 0; split_size, split_size];
                     sub_Q_Matrix = A3_generateQMatrix(split_size, QP);
    
                    for idx = 1:4
                        % Calculate relative and absolute positions for the sub-block
                        relSubRow = relative_offsets(idx, 1) + 1; % MATLAB indexing starts at 1
                        relSubCol = relative_offsets(idx, 2) + 1;
                        absSubRow = i + relative_offsets(idx, 1);
                        absSubCol = j + relative_offsets(idx, 2);
                
                        % Get the current sub-block
                        sub_block = original_frame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1);
                
                        % Perform intra prediction for sub-block
                        [sub_predicted, sub_mode, ~] = A3_intraPredictForIBlock(original_frame, absSubRow, absSubCol, ...
                            split_size, split_size, false, QP, lambda);
                
                        % Encode sub-block mode
                        diff_sub_mode = sub_mode - previous_mode;
                        previous_mode = sub_mode;
                        fprintf(Diff_stream, '%s\n', A1_Q4_expGolombEncode(diff_sub_mode));
                        bitConsumedThisFrame = bitConsumedThisFrame + A3_calcBitCounts(A1_Q4_expGolombEncode(diff_sub_mode));
    
                        % Calculate residual block
                        residual_block = sub_block - sub_predicted;
                
                        % Encode residual block
                        encoded_residual_block = A1_Q4_quantizeBlockAfterDCT(residual_block, sub_Q_Matrix);
                
                        % Perform scanning and RLE
                        scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
                        rle_encoded = A1_Q4_rleEncode(scanned_coeffs, split_size);
                
                        % Apply exp-Golomb encoding and write to QTC stream
                        encoded_value = A1_Q4_expGolombEncode(rle_encoded);
                        
                        [encoded_value, ~] = A3_quantizeAndEncode(residual_block, QP, false);
                        fprintf(QTC_stream, '%s\n', encoded_value);
                        bitConsumedThisFrame = bitConsumedThisFrame + A3_calcBitCounts(encoded_value);
    
                        % Reconstruct residual block (aligning with decoder)
                        decoded_rle = A1_Q4_expGolombDecode(encoded_value);  % Decode exp-Golomb
                        decoded_scanned = A1_Q4_rleDecode(decoded_rle, split_size);  % RLE decode
                        decoded_coeffs = A1_Q4_inverseSScan(decoded_scanned, split_size, split_size);  % Inverse scan
                        reconstructed_residual_block = A3_idctAfterDequantizeBlock(decoded_coeffs, sub_Q_Matrix);
                
                        % Reconstruct sub-block
                        sub_reconstructed = sub_predicted + reconstructed_residual_block;
                        sub_reconstructed = max(min(sub_reconstructed, 255), 0);
                
                        % Update reconstructed and predicted frames
                        reconstructedFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1) = sub_reconstructed;
                        predictedFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1) = sub_predicted;
                    end 
                else
                    predictedFrame(i:i+block_height-1, j:j+block_width-1) = predicted_block;
    
                    % encode the mode diff
                    differential_mode = mode - previous_mode;
                    previous_mode = mode;
                    encoded_diff = A1_Q4_expGolombEncode(differential_mode);
                    % split_flag, I_frame_flag, diff
                    fprintf(Diff_stream, '%s %s %s\n', A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(1), encoded_diff);
                    bitConsumedThisFrame = bitConsumedThisFrame + A3_calcBitCounts(A1_Q4_expGolombEncode(0));
                    bitConsumedThisFrame = bitConsumedThisFrame + A3_calcBitCounts(A1_Q4_expGolombEncode(1));
                    bitConsumedThisFrame = bitConsumedThisFrame + A3_calcBitCounts(encoded_diff);
             
                    % get the residuals of block
                    residual_block = block - predicted_block;
    
                    % encode the residual block
                    encoded_residual_block = A3_quantizeBlockAfterDCT(residual_block, Q_Matrix);
                    scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
                    rle_encoded = A1_Q4_rleEncode(scanned_coeffs, block_size);
                    encoded_value = A1_Q4_expGolombEncode(rle_encoded);
                    fprintf(QTC_stream, '%s\n', encoded_value);
                    bitConsumedThisFrame = bitConsumedThisFrame + A3_calcBitCounts(encoded_value);
    
                    % decode the encoded block and add it to prediction
                    decoded_residual_block = A3_idctAfterDequantizeBlock(encoded_residual_block, Q_Matrix);
                    reconstructed_block = predicted_block + decoded_residual_block;
                    reconstructed_block = max(min(reconstructed_block, 255), 0);
    
                    % store the encode residual block and reconstructed block
                    encoded_frame(i:i+block_height-1, j:j+block_width-1) = encoded_residual_block;
                    reconstructedFrame(i:i+block_height-1, j:j+block_width-1) = reconstructed_block;
                end
            end
        end
    end
end