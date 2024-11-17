function [predictedFrame, reconstructedFrame] = A2_Q2_intraPredictForIFrame(original_frame, block_size, QP, MDiff_stream, QTC_stream, VBSenable)
    Q_Matrix = A2_Q2_generateQMatrix(block_size, QP);
    [h, w] = size(original_frame);
    reconstructedFrame = zeros(h, w);
    lambda = A2_Q2_getLambda(QP);
    previous_mode = 0;

    for i = 1:block_size:size(original_frame, 1)
        for j = 1:block_size:size(original_frame, 2)
            block_height = min(block_size, h - i + 1);
            block_width = min(block_size, w - j + 1);

            block = original_frame(i:i+block_height-1, j:j+block_width-1);
            [predicted_block, mode, split_flag] = A2_Q2_intraPredictForBlock(original_frame, i, j, block_height, block_width, VBSenable, QP, lambda);

            if VBSenable && split_flag 
                % Encode split flag
                fprintf(MDiff_stream, '%s\n', A1_Q4_expGolombEncode(1)); % 1 indicates split

                % Process each sub-block
                split_size = floor(block_size / 2);
                sub_blocks = {
                    [i, j, split_size, split_size], ...
                    [i, j+split_size, split_size, split_size], ...
                    [i+split_size, j, split_size, split_size], ...
                    [i+split_size, j+split_size, split_size, split_size]
                };

                for k = 1:4
                    si = sub_blocks{k}(1);
                    sj = sub_blocks{k}(2);
                    sh = sub_blocks{k}(3);
                    sw = sub_blocks{k}(4);

                    sub_block = original_frame(si:si+sh-1, sj:sj+sw-1);
                    [sub_predicted, sub_mode, ~] = A2_Q2_intraPredictForBlock(original_frame, si, sj, sh, sw, false, QP-1, lambda);

                    % Encode sub-block mode
                    diff_sub_mode = sub_mode - previous_mode;
                    previous_mode = sub_mode;
                    fprintf(MDiff_stream, '%s\n', A1_Q4_expGolombEncode(diff_sub_mode));

                    % Encode sub-block residual
                    sub_residual = sub_block - sub_predicted;
                    [sub_encoded, sub_quantized] = A2_Q2_quantizeAndEncode(sub_residual, QP-1);
                    fprintf(QTC_stream, '%s\n', sub_encoded);

                    % Reconstruct sub-block
                    sub_residual_decoded = A2_Q2_idctAfterDequantizeBlock(sub_quantized, A2_Q2_generateQMatrix(split_size, QP-1));
                    sub_reconstructed = sub_predicted + sub_residual_decoded;
                    sub_reconstructed = max(min(sub_reconstructed, 255), 0);

                    % Update predicted and reconstructed frames
                    predictedFrame(si:si+sh-1, sj:sj+sw-1) = sub_predicted;
                    reconstructedFrame(si:si+sh-1, sj:sj+sw-1) = sub_reconstructed;
                end
            else
                % Non-split logic: Encode split flag, frame type (I-frame), and mode difference
                block_height = min(block_size, h - i + 1);
                block_width = min(block_size, w - j + 1);
       
                % get the current block
                block = original_frame(i:i+block_height-1, j:j+block_width-1);
    
                % Perform intra prediction of current block
                [predicted_block, mode] = A2_Q2_intraPredictForBlock(original_frame, i, j, block_height, block_width, VBSenable);
    
                predictedFrame(i:i+block_height-1, j:j+block_width-1) = predicted_block;
    
                differential_mode = mode - previous_mode;
                previous_mode = mode;
                %MDiff_stream = [MDiff_stream, A1_Q4_expGolombEncode(differential_mode)];
                encoded_diff = A1_Q4_expGolombEncode(differential_mode);
                fprintf(MDiff_stream, '%s %s %s\n',A1_Q4_expGolombEncode(0), A1_Q4_expGolombEncode(1), encoded_diff);
     
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
            end
        end
    end
end