function residualBlock_de = A1_Q3_readResidualFile(residualFile, blockSize)
    % Read the residual block from the residualFile
    fgetl(residualFile); %skip one line
    residualLine = fgetl(residualFile);  % Read the residual block line
    residualValues = sscanf(residualLine, '%d');  % Extract residual values
    residualBlock_de = reshape(residualValues, [blockSize, blockSize]);
end