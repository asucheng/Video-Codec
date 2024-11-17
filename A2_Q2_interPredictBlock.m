function [bestMatch, predictedBlock] = A2_Q2_interPredictBlock(referenceFrame, currentBlock, row, col, searchRange, blockSize, paddedWidth, paddedHeight)
    % Initialize
    bestSAD = inf; % Initialize best SAD (Sum of Absolute Differences) to a large value
    bestMatch = [0, 0]; % Initialize best motion vector (row_offset, col_offset)

    blockHeight = size(currentBlock, 1);
    blockWidth = size(currentBlock, 2);

    % Full search within the search range
    for rowOffset = -searchRange:searchRange
        for colOffset = -searchRange:searchRange
            % Calculate the reference block position
            refRow = row + rowOffset;
            refCol = col + colOffset;

            % Check bounds
            if refRow < 1 || refCol < 1 || refRow + blockHeight - 1 > paddedHeight || refCol + blockWidth - 1 > paddedWidth
                continue; % Skip out-of-bound reference blocks
            end

            % Extract the reference block
            referenceBlock = referenceFrame(refRow:refRow+blockHeight-1, refCol:refCol+blockWidth-1);

            % Calculate SAD (Sum of Absolute Differences)
            SAD = sum(abs(double(currentBlock(:)) - double(referenceBlock(:))));

            % Update the best match if a lower SAD is found
            if SAD < bestSAD
                bestSAD = SAD;
                bestMatch = [rowOffset, colOffset];
            end
        end
    end

    % Use the best motion vector to predict the block
    bestRow = row + bestMatch(1);
    bestCol = col + bestMatch(2);
    predictedBlock = referenceFrame(bestRow:bestRow+blockHeight-1, bestCol:bestCol+blockWidth-1);
end