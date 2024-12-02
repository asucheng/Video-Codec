function [bitBudgetPerRowBlock] = A3_getRowBitBudget(RCflag, targetBR, rowBlockCounts, fps)
    if nargin < 4
        fps = 1;
    end

    if RCflag
        bitBudgetPerRowBlock = round(targetBR / fps / rowBlockCounts, 0);
    else
        bitBudgetPerRowBlock = inf; % means no limitation to bit count
    end
end
