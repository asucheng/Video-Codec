function quantized_block = A2_Q3_quantizeBlockAfterDCT(block, Q_matrix)
    dct_block = dct2(block);
    quantized_block = round(dct_block ./ Q_matrix);
end