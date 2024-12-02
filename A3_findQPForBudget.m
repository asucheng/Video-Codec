function [QP, avgBitCount] = A3_findQPForBudget(targetBudget, targetBlockSize, ...
    blockSizes, QP_values, budgetTable)
    
    % QP_Table is table-like:
    %   [ x x x x x x x x x]
    %   [ x x x x x x x x x]
    %   [ x x x x x x x x x]
    %   [ x x x x x x x x x]

    row_index = targetBlockSize == blockSizes;
    row_data = budgetTable(row_index, :);

    
    for idx = 1:length(row_data)
        if row_data(idx) <= targetBudget
            QP = QP_values(idx);
            avgBitCount = budgetTable(row_index, idx);
            return
        end
    end

    QP = -1; % represent no available found.
    avgBitCount = inf;
end