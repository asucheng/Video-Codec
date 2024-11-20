function [J, encoded_block] = A2_Q2_computeRD(residual_block, mode, QP, Lambda)
    % Step 1: Perform DCT on the residual block
    dct_block = dct2(residual_block);

    % Step 2: Generate Q_Matrix and quantize the block
    block_size = size(dct_block, 1);
    Q_Matrix = A2_Q2_generateQMatrix(block_size, QP);
    quantized_block = round(dct_block ./ Q_Matrix);

    % Step 3: Encode the quantized block using the expanded logic
    % Quantize, scan, and encode the block
    encoded_residual_block = A1_Q4_quantizeBlockAfterDCT(residual_block, Q_Matrix);
    scanned_coeffs = A1_Q4_sScan(encoded_residual_block);
    rle_encoded = A1_Q4_rleEncode(scanned_coeffs, block_size);
    encoded_block = A1_Q4_expGolombEncode(rle_encoded);

    % Step 4: Rescale the quantized block to approximate original coefficients
    rescaled_block = quantized_block .* Q_Matrix; % Multiply quantized coefficients by Q matrix

    % Step 5: Perform inverse DCT to approximate the residual
    approximated_residual = idct2(rescaled_block);

    % Step 6: Compute distortion (SSD between original and reconstructed residuals)
    distortion = sum((residual_block(:) - approximated_residual(:)).^2);

    % Step 7: Compute rate as the length of the encoded string
    rate = strlength(encoded_block) + strlength(A1_Q4_expGolombEncode(mode)); % Add mode encoding cost

    % Step 8: Compute RD cost
    J = distortion + Lambda * rate;
end