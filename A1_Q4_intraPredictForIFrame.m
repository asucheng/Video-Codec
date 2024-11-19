function [predictedFrame, reconstructedFrame] = A1_Q4_intraPredictForIFrame(original_frame, block_size, QP, MDiff_stream, QTC_stream, VBSenable)
    Q_Matrix = A1_Q4_generateQMatrix(block_size, QP);

    [h, w] = size(original_frame);
    reconstructedFrame = zeros(h, w);
    predictedFrame = zeros(h, w);
    previous_mode = 0;
    Lambda = A2_Q2_getLambda(QP)
    % Process the frame in blocks
    for i = 1:block_size:h
        for j = 1:block_size:w
            % Determine the actual size of the block (boundary handling)
            block_height = min(block_size, h - i + 1);
            block_width = min(block_size, w - j + 1);
        
            % Get the current block
            block = original_frame(i:i+block_height-1, j:j+block_width-1);

            % Perform intra prediction with split/non-split based on VBSenable
            if VBSenable && block_height > 1 && block_width > 1
                % Use the advanced VBS-enabled block prediction
                [predicted_block, mode, split_flag] = A2_Q2_intraPredictForBlock(original_frame, i, j, block_height, block_width, VBSenable, QP, Lambda);
            else
                % Use non-split prediction only
                [predicted_block, mode] = A2_Q2_nonSplitPrediction(original_frame, i, j, block_height, block_width);
                split_flag = false;
            end

            % Update the predicted frame
            predictedFrame(i:i+block_height-1, j:j+block_width-1) = predicted_block;

            % Differential encode the mode
            differential_mode = mode - previous_mode;
            previous_mode = mode;
            encoded_diff = A1_Q4_expGolombEncode(differential_mode);
            fprintf(MDiff_stream, '%s %s\n', A1_Q4_expGolombEncode(1), encoded_diff);

            % Compute residuals
            residual_block = block - predicted_block;

            % Quantize and encode the residual block
            encoded_residual_block = A1_Q4_quantizeBlockAfterDCT(residual_block, Q_Matrix);
            scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
            rle_encoded = A1_Q4_rleEncode(scanned_coeffs, block_size);
            encoded_value = A1_Q4_expGolombEncode(rle_encoded);
            fprintf(QTC_stream, '%s\n', encoded_value);

            % Decode the residual block and reconstruct
            decoded_residual_block = A1_Q4_idctAfterDequantizeBlock(encoded_residual_block, Q_Matrix);
            reconstructed_block = predicted_block + decoded_residual_block;
            reconstructed_block = max(min(reconstructed_block, 255), 0);

            % Store the reconstructed block
            reconstructedFrame(i:i+block_height-1, j:j+block_width-1) = reconstructed_block;
        end
    end
end