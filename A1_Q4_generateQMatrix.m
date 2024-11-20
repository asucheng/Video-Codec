function Q_Matrix = A1_Q4_generateQMatrix(i, QP)
    % Generate Q matrix based on QP
    Q_Matrix = zeros(i, i);

    for x = 1:i
        for y = 1:i
            if (x + y < i + 1)
                Q_Matrix(x, y) = 2^QP;
            elseif (x + y == i + 1)
                Q_Matrix(x, y) = 2^(QP + 1);
            else
                Q_Matrix(x, y) = 2^(QP + 2);
            end
        end
    end
end
