function [predictedFrame, reconstructedFrame] = A2_Q34_intraPredictForIFrame(original_frame, block_size, QP, ...
    MDiff_stream, MVPDiff_stream, QTC_stream, VBSEnable, FastME)

    Q_Matrix = A1_Q4_generateQMatrix(block_size, QP);
    [h, w] = size(original_frame);
    %reconstructed_frame = zeros(h, w);
    encoded_frame = zeros(h, w);
    previous_mode = 0;

    % Simulate processing the frame in blocks
    for i = 1:block_size:size(original_frame, 1)
        for j = 1:block_size:size(original_frame, 2)

            % get the real size of current block
            block_height = min(block_size, h - i + 1);
            block_width = min(block_size, w - j + 1);
   
            % get the current block
            block = original_frame(i:i+block_height-1, j:j+block_width-1);

            % Perform intra prediction of current block
            [predicted_block, mode, split_flag] = A2_Q2_intraPredictForBlock(original_frame, i, j, block_height, block_width,VBSEnable, QP);
            if split_flag == 0
                predictedFrame(i:i+block_height-1, j:j+block_width-1) = predicted_block;
    
                differential_mode = mode - previous_mode;
                previous_mode = mode;
                %MDiff_stream = [MDiff_stream, A1_Q4_expGolombEncode(differential_mode)];
                encoded_diff = A1_Q4_expGolombEncode(differential_mode);
    
                if FastME
                    fprintf(MVPDiff_stream, '%s %s %s\n',A1_Q4_expGolombEncode(split_flag), A1_Q4_expGolombEncode(1), encoded_diff);
                else
                    fprintf(MDiff_stream, '%s %s %s\n',A1_Q4_expGolombEncode(split_flag), A1_Q4_expGolombEncode(1), encoded_diff);
                end
     
                % get the residuals of block
                residual_block = block - predicted_block;
    
                % encode the residual block
                encoded_residual_block = A1_Q4_quantizeBlockAfterDCT(residual_block, Q_Matrix);
    
                scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
                rle_encoded = A1_Q4_rleEncode(scanned_coeffs, block_size);
                %QTC_stream = [QTC_stream, A1_Q4_expGolombEncode(rle_encoded)];
                encoded_value = A1_Q4_expGolombEncode(rle_encoded);
                fprintf(QTC_stream, '%s\n', encoded_value);
    
                % decode the encoded block and add it to prediction
                decoded_residual_block = A1_Q4_idctAfterDequantizeBlock(encoded_residual_block, Q_Matrix);
                reconstructed_block = predicted_block + decoded_residual_block;
                reconstructed_block = max(min(reconstructed_block, 255), 0);
    
                % store the encode residual block and reconstructed block
                encoded_frame(i:i+block_height-1, j:j+block_width-1) = encoded_residual_block;
                reconstructedFrame(i:i+block_height-1, j:j+block_width-1) = reconstructed_block;
            else
                if FastME
                    fprintf(MVPDiff_stream, '%s %s\n',A1_Q4_expGolombEncode(split_flag), A1_Q4_expGolombEncode(1));
                else
                    fprintf(MDiff_stream, '%s %s\n',A1_Q4_expGolombEncode(split_flag), A1_Q4_expGolombEncode(1));
                end
                % Process each sub-block in Z-order
                split_size = floor(block_size / 2);
                relative_offsets = [0, 0; 0, split_size; split_size, 0; split_size, split_size];
                sub_Q_Matrix = A2_Q2_generateQMatrix(split_size, QP);
            
                for idx = 1:4
                    % Calculate relative and absolute positions for the sub-block
                    absSubRow = i + relative_offsets(idx, 1);
                    absSubCol = j + relative_offsets(idx, 2);
            
                    % Get the current sub-block
                    sub_block = original_frame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1);
            
                    % Perform intra prediction for sub-block
                    [sub_predicted, sub_mode, ~] = A2_Q2_intraPredictForBlock(original_frame, absSubRow, absSubCol, split_size, split_size, false, QP);
            
                    % Encode sub-block mode
                    diff_sub_mode = sub_mode - previous_mode;
                    previous_mode = sub_mode;
                    fprintf(MDiff_stream, '%s\n', A1_Q4_expGolombEncode(diff_sub_mode));
            
                    % Calculate residual block
                    residual_block = sub_block - sub_predicted;
            
                    % Encode residual block
                    encoded_residual_block = A1_Q4_quantizeBlockAfterDCT(residual_block, sub_Q_Matrix);
            
                    % Perform scanning and RLE
                    scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
                    rle_encoded = A1_Q4_rleEncode(scanned_coeffs, split_size);
            
                    % Apply exp-Golomb encoding and write to QTC stream
                    encoded_value = A1_Q4_expGolombEncode(rle_encoded);
                    fprintf(QTC_stream, '%s\n', encoded_value);
            
                    % Reconstruct residual block (aligning with decoder)
                    decoded_rle = A1_Q4_expGolombDecode(encoded_value);  % Decode exp-Golomb
                    decoded_scanned = A1_Q4_rleDecode(decoded_rle, split_size);  % RLE decode
                    decoded_coeffs = A1_Q4_inverseSScan(decoded_scanned, split_size, split_size);  % Inverse scan
                    reconstructed_residual_block = A1_Q4_idctAfterDequantizeBlock(decoded_coeffs, sub_Q_Matrix);
            
                    % Reconstruct sub-block
                    sub_reconstructed = sub_predicted + reconstructed_residual_block;
                    sub_reconstructed = max(min(sub_reconstructed, 255), 0);
            
                    % Update reconstructed and predicted frames
                    reconstructedFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1) = sub_reconstructed;
                    predictedFrame(absSubRow:absSubRow+split_size-1, absSubCol:absSubCol+split_size-1) = sub_predicted;
                end
            end
            % fprintf('Block (%d, %d) - Selected mode: %d\n', i, j, mode);
        end
    end
end