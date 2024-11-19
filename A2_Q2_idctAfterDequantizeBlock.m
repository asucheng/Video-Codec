function rescaled_block = A2_Q2_idctAfterDequantizeBlock(quantized_block, Q_Matrix)
    rescaled_coeffs = quantized_block .* Q_Matrix;
    rescaled_block = idct2(rescaled_coeffs);
end