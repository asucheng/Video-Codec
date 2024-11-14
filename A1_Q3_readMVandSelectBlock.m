function [refRow, refCol] = A1_Q3_readMVandSelectBlock(mvFile, row, col)
    % Read the motion vector from the mvFile
    mvLine = fgetl(mvFile);  % Read a line from motion vector file
    % Use sscanf with a single output
    parsedValues = sscanf(mvLine, 'Frame %d, Block (%d, %d): MV = (%d, %d)');
    
    % Extract individual values from the parsed vector
    frameNum = parsedValues(1);
    blockRow = parsedValues(2);
    blockCol = parsedValues(3);
    xOffset = parsedValues(4);
    yOffset = parsedValues(5);
    
    % Use the motion vector to get the predictor block from the reference frame
    refRow = row + yOffset;
    refCol = col + xOffset;
end