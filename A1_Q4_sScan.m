function reordered = A1_Q4_sScan(matrix)
    % Reorder a matrix in zigzag order (diagonal scan)
    [rows, cols] = size(matrix);
    reordered = zeros(1, rows * cols);
    index = 1;

    for s = 1:(rows + cols - 1)
        if s <= cols
            i = 1;
            j = s;
        else
            i = s - cols + 1;
            j = cols;
        end

        while i <= rows && j >= 1
            reordered(index) = matrix(i, j);
            index = index + 1;
            i = i + 1; 
            j = j - 1;
        end
    end
end
