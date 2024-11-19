function approxResidualBlock = A1_Q3_calcResidual(predictedBlock, currentBlock, n)
    % Residual block
    residualBlock = double(currentBlock) - double(predictedBlock);
    
    % Approximated residuals for n = 1, 2, 3
    approxResidualBlock = round(residualBlock / (2^n)) * (2^n);
    
    % Save the approximated residual
    % fprintf(residualFile, 'Frame %d, Block (%d, %d), n=%d:\n', frameIdx, row, col, n);
    % fprintf(residualFile, '%d ', approxResidualBlock(:));
    % fprintf(residualFile, '\n');
end