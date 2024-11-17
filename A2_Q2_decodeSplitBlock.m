function reconstructedBlock = A2_Q2_decodeSplitBlock(referenceFrame, QTC_Coeff_file, MDiff_file, ...
                                                     row, col, blockSize, paddedHeight, paddedWidth, QP)
    % Calculate the size of the split sub-block
    splitSize = blockSize / 2;
    reconstructedBlock = zeros(blockSize, blockSize, 'double');  % Initialize the reconstructed block

    % Define relative offsets for Z-order traversal within the parent block
    relative_offsets = [0, 0; 0, splitSize; splitSize, 0; splitSize, splitSize];
    Q_Matrix = A2_Q2_generateQMatrix(splitSize, QP); % Pre-compute Q matrix for sub-blocks

    for idx = 1:4
        relSubRow = relative_offsets(idx, 1) + 1; % +1 for MATLAB indexing
        relSubCol = relative_offsets(idx, 2) + 1;

        % Calculate absolute positions in the padded frame
        absSubRow = row + relative_offsets(idx, 1);
        absSubCol = col + relative_offsets(idx, 2);

        % Check if sub-block is within bounds of the padded frame
        if absSubRow + splitSize - 1 > paddedHeight || absSubCol + splitSize - 1 > paddedWidth
            continue;
        end

        % Read the mode from MDiff file (no split flag needed for sub-blocks)
        MDiff_line = fgetl(MDiff_file);
        if isempty(MDiff_line)
            error('Unexpected end of MDiff file');
        end
        MDiff_line_array = A1_Q4_expGolombDecode(MDiff_line);
        mode = MDiff_line_array(1);  % Mode for sub-block

        % Generate the predictor block
        block_height = min(splitSize, paddedHeight - absSubRow + 1);
        block_width = min(splitSize, paddedWidth - absSubCol + 1);

        if mode == 0
            % Horizontal prediction
            if absSubRow == 1
                predictorBlock = 128 * ones(block_height, block_width);
            else
                predictorBlock = repmat(referenceFrame(absSubRow-1, absSubCol:absSubCol+block_width-1), block_height, 1);
            end
        else
            % Vertical prediction
            if absSubCol == 1
                predictorBlock = 128 * ones(block_height, block_width);
            else
                predictorBlock = repmat(referenceFrame(absSubRow:absSubRow+block_height-1, absSubCol-1), 1, block_width);
            end
        end

        % Read and decode residuals
        QTC_Line = fgetl(QTC_Coeff_file);
        % if isempty(QTC_Line)
        %     error('Unexpected end of QTC_Coeff file');
        % end
        encoded_rle = A1_Q4_expGolombDecode(QTC_Line);
        coeffs_scanned = A1_Q4_rleDecode(encoded_rle, splitSize);
        residual_block_encoded = A1_Q4_inverseSScan(coeffs_scanned, splitSize, splitSize);
        residualBlock_de = A2_Q2_idctAfterDequantizeBlock(residual_block_encoded, Q_Matrix);

        % Reconstruct the sub-block
        reconstructedSubBlock = double(predictorBlock) + double(residualBlock_de);
        reconstructedSubBlock = max(min(reconstructedSubBlock, 255), 0);

        % Place the sub-block in the relative location within the parent block
        reconstructedBlock(relSubRow:relSubRow+block_height-1, relSubCol:relSubCol+block_width-1) = reconstructedSubBlock;
    end
end