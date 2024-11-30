function quantized_block = A2_quantizeBlockAfterDCT(block, Q_matrix)
    dct_block = dct2(block);
    quantized_block = round(dct_block ./ Q_matrix);
end