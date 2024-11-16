function [predictedFrame, reconstructedFrame] = A1_Q4_intraPredictForIFrame(original_frame, block_size, QP, MDiff_stream, QTC_stream)

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
            [predicted_block, mode] = A1_Q4_intraPredictForBlock(original_frame, i, j, block_height, block_width);

            predictedFrame(i:i+block_height-1, j:j+block_width-1) = predicted_block;

            differential_mode = mode - previous_mode;
            previous_mode = mode;
            %MDiff_stream = [MDiff_stream, A1_Q4_expGolombEncode(differential_mode)];
            encoded_diff = A1_Q4_expGolombEncode(differential_mode);
            fprintf(MDiff_stream, '%s %s\n', A1_Q4_expGolombEncode(1), encoded_diff);
 
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

            % fprintf('Block (%d, %d) - Selected mode: %d\n', i, j, mode);
        end
    end
end