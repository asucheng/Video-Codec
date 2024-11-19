function Q_Matrix = A2_Q2_generateQMatrix(block_size, QP)
    Q_Matrix = zeros(block_size, block_size);

    for x = 1:block_size
        for y = 1:block_size
            if x + y < block_size + 1
                Q_Matrix(x, y) = 2^QP;
            elseif x + y == block_size + 1
                Q_Matrix(x, y) = 2^(QP + 1);
            else
                Q_Matrix(x, y) = 2^(QP + 2);
            end
        end
    end
end