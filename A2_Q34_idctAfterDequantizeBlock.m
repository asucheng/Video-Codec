function idct_block = A2_Q3_idctAfterDequantizeBlock(block, Q_matrix)
    deq_block = block .* Q_matrix;
    idct_block = idct2(deq_block);
end