function [bestMatch_frame, best_predictedBlk, best_ref_index] = A2_interPredictForPBlock(reference_frames, currentBlock, ...
    row, col, searchRange, blockSize, paddedWidth, paddedHeight, ...
    nRefFrames, FMEEnable, FastME, mvp)

    % MRF parameter for frame
    bestMAE_frame = inf;
    best_ref_index = 1;
    num_ava_refs = min(length(reference_frames), nRefFrames);
    bestMatch_frame = [0, 0];

    % parameter for block
    bestMAE = inf;
    bestMatch = [0, 0];
    bestPredictedBlock = zeros(blockSize, blockSize, 'uint8');

    for ref_idx = 1:num_ava_refs
        referenceFrame = reference_frames{ref_idx};

        if FastME
            candidates = [
                mvp,
                mvp + [0, 1],
                mvp + [0, -1],
                mvp + [1, 0],
                mvp + [-1, 0]
            ];
        else
            [dy, dx] = meshgrid(-searchRange: searchRange, -searchRange: searchRange);
            candidates = [dy(:), dx(:)];
        end

        % Full search within search range
        for i = 1:size(candidates, 1)
            mv = candidates(i, :);
            refRow = row + mv(1);
            refCol = col + mv(2);

            if refRow < 1 || refCol < 1 || refRow+blockSize-1 > paddedHeight || refCol+blockSize-1 > paddedWidth
                continue;
            end

            % Reference block
            referenceBlock = referenceFrame(refRow:refRow+blockSize-1, refCol:refCol+blockSize-1);

            if FMEEnable
                referenceBlock = A2_performInterpolation(referenceFrame, mv, refRow, refCol, blockSize);
            end
            mae = mean(abs(double(currentBlock(:)) - double(referenceBlock(:))));

            % Update best match based on MAE and motion vector criteria
            if mae < bestMAE || ...
            (mae == bestMAE && abs(mv(2)) + abs(mv(1)) < abs(bestMatch(1)) + abs(bestMatch(2))) || ...
            (mae == bestMAE && abs(mv(2)) + abs(mv(1)) == abs(bestMatch(1)) + abs(bestMatch(2)) && ...
            (mv(1) < bestMatch(2) || (mv(1) == bestMatch(2) && mv(2) < bestMatch(1))))
                bestMAE = mae;
                bestMatch = mv;
                bestPredictedBlock = referenceBlock;
            end
        end

        if bestMAE < bestMAE_frame || ...
            (bestMAE == bestMAE_frame && abs(bestMatch(1)) + abs(bestMatch(2)) < abs(bestMatch_frame(1)) + abs(bestMatch_frame(2))) || ...
            (bestMAE == bestMAE_frame && abs(bestMatch(1)) + abs(bestMatch(2)) == abs(bestMatch_frame(1)) + abs(bestMatch_frame(2)) && ...
            (bestMatch(2) < bestMatch_frame(2) || (bestMatch(2) == bestMatch_frame(2) && bestMatch_frame(1) < bestMatch_frame(1))))
            bestMAE_frame = bestMAE;
            bestMatch_frame = bestMatch;
            best_ref_index = ref_idx;
            best_predictedBlk = bestPredictedBlock;
        end
    end
end