function [refRow, refCol, previous_mv] = A1_Q3_readMVandSelectBlock(MDiff_file, row, col, previous_mv)
    % Read the motion vector from the mvFile
    % mvLine = fgetl(mvFile);  % Read a line from motion vector file
    % Use sscanf with a single output
    % parsedValues = sscanf(mvLine, 'Frame %d, Block (%d, %d): MV = (%d, %d)');
    
    % Extract individual values from the parsed vector
    % frameNum = parsedValues(1);
    % blockRow = parsedValues(2);
    % blockCol = parsedValues(3);
    % xOffset = parsedValues(4);
    % yOffset = parsedValues(5);
    
    % Use the motion vector to get the predictor block from the reference frame
    % refRow = row + bestMV(2); %row + yOffset; 
    % refCol = col + bestMV(1); %col + xOffset;

    % read the a line of MDiff file
    MDiff_line = fgetl(MDiff_file);
    diff_mv = A1_Q4_expGolombDecode(MDiff_line);
    
    % update the motion vector differentiation
    bestMatch = diff_mv + previous_mv;
    previous_mv = bestMatch;

    % Use the motion vector to get the predictor block from the reference frame
    refRow = row + bestMatch(2); %row + yOffset; 
    refCol = col + bestMatch(1); %col + xOffset;
end