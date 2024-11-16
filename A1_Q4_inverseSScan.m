function matrix = A1_Q4_inverseSScan(reordered, rows, cols)
    % Reconstruct a matrix from zigzag order (diagonal scan) back to 2D
    matrix = zeros(rows, cols);  % Initialize the matrix
    index = 1;

    for s = 1:(rows + cols - 1)
        if s <= cols
            i = 1;
            j = s;
        else
            i = s - cols + 1;
            j = cols;
        end

        % Assign each value from reordered back to matrix in zigzag order
        while i <= rows && j >= 1
            matrix(i, j) = reordered(index);
            index = index + 1;
            i = i + 1;
            j = j - 1;
        end
    end
end