function [J, encoded_block] = A2_Q2_computeRD(residual_block, mode, QP, Lambda)
    % Quantize and encode residuals
    [encoded_block, quantized_block] = A2_Q2_quantizeAndEncode(residual_block, QP);

    % Explicit rescaling to approximate original coefficients
    i = size(residual_block, 1); % Block size
    Q_Matrix = A2_Q2_generateQMatrix(i, QP); % Generate Q matrix
    rescaled_block = quantized_block .* Q_Matrix; % Multiply quantized coefficients by Q matrix
    
    % Inverse DCT to approximate the residual
    approximated_residual = idct2(rescaled_block);

    % Compute distortion (SSD between original and reconstructed residuals)
    distortion = sum((residual_block(:) - approximated_residual(:)).^2);

    % Compute rate as the length of the encoded string
    rate = strlength(encoded_block) + strlength(A1_Q4_expGolombEncode(mode)); % Add mode encoding cost

    % Compute RD cost
    J = distortion + Lambda * rate;
end