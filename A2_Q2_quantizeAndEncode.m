function [encoded_block, quantized_block] = A2_Q2_quantizeAndEncode(residual_block, QP, VBSenable)
    dct_block = dct2(residual_block);
    Q_Matrix = A2_Q2_generateQMatrix(size(dct_block, 1), QP);
    quantized_block = round(dct_block ./ Q_Matrix);

    scanned_coeffs = A1_Q4_sScan(quantized_block);
    rle_encoded = A1_Q4_rleEncode(scanned_coeffs, size(dct_block, 1));
    encoded_block = A1_Q4_expGolombEncode(rle_encoded);

end