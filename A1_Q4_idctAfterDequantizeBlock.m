function idct_block = A1_Q4_idctAfterDequantizeBlock(block, Q_matrix)
    deq_block = block .* Q_matrix;
    idct_block = idct2(deq_block);
end