function psnrValues_verify = A3_Decoding(filename_prefix, nframes, blockSize, paddedWidth, paddedHeight, ...
    height, width, QP, nRefFrames, FMEEnable, FastME)

    % Decoder for Y-only frames using motion vectors and residuals
    psnrValues_verify = zeros(1, nframes);  % Store PSNR for each frame, decode vs reconstructed
    
    % Initialize the hypothetical reference frame for the first frame (all values = 128)
    referenceFrame_de = 128 * ones(paddedHeight, paddedWidth, 'uint8');
    
    % Open the necessary files
    splits = fopen(strcat(filename_prefix, 'splits.txt'), 'w');
    vid_decoded = fopen(strcat(filename_prefix, 'decoded_vid.yuv'), 'w');
    vid_reconstructed = fopen(strcat(filename_prefix, 'reconstructed_vid.yuv'), 'r');

    MDiff_file = fopen(strcat(filename_prefix, "MDiff.txt"), 'r');
    MVPDiff_file = fopen(strcat(filename_prefix, "MVPDiff.txt"), 'r');
    QTC_Coeff_file = fopen(strcat(filename_prefix, "QTC_Coeff.txt"), 'r');

    if FastME
        MDiff_file = MVPDiff_file;
    end
    
    % Loop through frames
    for frameIdx = 1:nframes
        % Initialize the decoded frame
        decodedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
        Q_Matrix = A3_generateQMatrix(blockSize, QP);

        previous_mv = [0, 0];
        previous_mode = 0;
        previous_ref_index = 0;

        % Loop over blocks in raster order
        for row = 1:blockSize:paddedHeight
            for col = 1:blockSize:paddedWidth
                % Check for boundary issues (skip if block exceeds frame size)
                if row + blockSize - 1 > paddedHeight || col + blockSize - 1 > paddedWidth
                    continue;
                end

                % read the a line of MDiff file----------------------------
                MDiff_line = fgetl(MDiff_file);
                if isempty(MDiff_line)
                    error('Unexpected end of MDiff file.');
                end
                if ~(ischar(MDiff_line) || isstring(MDiff_line))
                    disp('Error: MDiff_line is not a valid string.');
                    error('MDiff_line must be a character array or string.');
                end
                if ~strcmp(MDiff_line, strtrim(MDiff_line))
                    disp('Error: MDiff_line contains leading or trailing whitespace.');
                end

                string = strtrim(MDiff_line);
                MDiff_line_array = A1_Q4_expGolombDecode(MDiff_line);

                split_flag = MDiff_line_array(1);
                I_frame_flag = MDiff_line_array(2);
                fprintf(splits, '%d\n', split_flag);
 
                if split_flag == 1
                    % Handle split blocks inline (logic from A3_Q2_decodeSplitBlock)
                    splitSize = blockSize / 2;
                    % Define relative offsets for Z-order traversal within the parent block
                    relative_offsets = [0, 0; 0, splitSize; splitSize, 0; splitSize, splitSize];
                    sub_Q_Matrix = A3_generateQMatrix(splitSize, QP); % Pre-compute Q matrix for sub-blocks
                    for idx = 1:4
                        % Calculate absolute positions in the padded frame
                        absSubRow = row + relative_offsets(idx, 1);
                        absSubCol = col + relative_offsets(idx, 2);
                
                        % Check if sub-block is within bounds of the padded frame
                        if absSubRow + splitSize - 1 > paddedHeight || absSubCol + splitSize - 1 > paddedWidth
                            continue;
                        end
                
                        % Read the mode from MDiff file (no split flag for sub-blocks)
                        MDiff_line = fgetl(MDiff_file);
                        if isempty(MDiff_line)
                            error('Unexpected end of MDiff file');
                        end
                        
                        % Generate the predictor block
                        block_height = min(splitSize, paddedHeight - absSubRow + 1);
                        block_width = min(splitSize, paddedWidth - absSubCol + 1);
                
                        if I_frame_flag == 1
                            MDiff_line_array = A1_Q4_expGolombDecode(MDiff_line);
                            diff_mode = MDiff_line_array(1);  % Mode for sub-block
                            mode = previous_mode + diff_mode;
                            previous_mode = mode;

                            referenceFrames_de = [];
                            % I-frame prediction logic
                            if mode == 0
                                % Horizontal prediction
                                if absSubRow == 1
                                    predictorBlock = 128 * ones(block_height, block_width);
                                else
                                    predictorBlock = repmat(decodedFrame(absSubRow-1, absSubCol:absSubCol+block_width-1), block_height, 1);
                                end
                            else
                                % Vertical prediction
                                if absSubCol == 1
                                    predictorBlock = 128 * ones(block_height, block_width);
                                else
                                    predictorBlock = repmat(decodedFrame(absSubRow:absSubRow+block_height-1, absSubCol-1), 1, block_width);
                                end
                            end
                        else
                            % P-frame prediction logic
                            MDiff_line_array = A1_Q4_expGolombDecode(MDiff_line);
                            diff_ref_index = MDiff_line_array(3);
                            best_ref_index = diff_ref_index + previous_ref_index;
                            previous_ref_index = best_ref_index;
                            referenceFrame_de = referenceFrames_de{best_ref_index};

                            diff_mv = MDiff_line_array(1:2); % Extract motion vector difference
                            bestMatch = diff_mv + previous_mv; % Compute absolute motion vector
                            previous_mv = bestMatch;
    
                            if FMEEnable
                                int_dy = floor(bestMatch(2) / 2);
                                int_dx = floor(bestMatch(1) / 2);
                                frac_dy = mod(bestMatch(2), 2) / 2;
                                frac_dx = mod(bestMatch(1), 2) / 2;
        
                                refRow = absSubRow + int_dy;
                                refCol = absSubCol + int_dx;
                                refRow = max(1, min(refRow, paddedHeight - splitSize + 1));
                                refCol = max(1, min(refCol, paddedWidth - splitSize + 1));
        
                                predictorBlock = performReverseInterpolation(referenceFrame_de, refRow, refCol, frac_dy, frac_dx, splitSize);
                            else
                                % Use the motion vector to get the predictor block from the reference frame
                                refRow = absSubRow + bestMatch(2); % row + yOffset; 
                                refCol = absSubCol + bestMatch(1); % col + xOffset;
                                
                                % Check for boundary in reference frame
                                if refRow < 1 || refCol < 1 || refRow + splitSize - 1 > paddedHeight || refCol + splitSize - 1 > paddedWidth
                                    predictorBlock = 128 * ones(block_height, block_width); % Skip if the predictor block goes out of bounds
                                else
                                    predictorBlock = referenceFrame_de(refRow:refRow+block_height-1, refCol:refCol+block_width-1);
                                end
                            end
                        end
                
                        % Read and decode residuals
                        QTC_Line = fgetl(QTC_Coeff_file);
                        if isempty(QTC_Line)
                            error('Unexpected end of QTC_Coeff file');
                        end
                        encoded_rle = A1_Q4_expGolombDecode(QTC_Line);
                        coeffs_scanned = A1_Q4_rleDecode(encoded_rle, splitSize);
                        residual_block_encoded = A1_Q4_inverseSScan(coeffs_scanned, splitSize, splitSize);
                        residualBlock_de = A1_Q4_idctAfterDequantizeBlock(residual_block_encoded, sub_Q_Matrix);

                        % Reconstruct the sub-block
                        reconstructedSubBlock = double(predictorBlock) + double(residualBlock_de);
                        reconstructedSubBlock = max(min(reconstructedSubBlock, 255), 0);
                
                        % Update the decoded frame directly
                        decodedFrame(absSubRow:absSubRow+block_height-1, absSubCol:absSubCol+block_width-1) = reconstructedSubBlock;
                    end
                else 
                    block_height = min(blockSize, paddedHeight - row + 1);
                    block_width = min(blockSize, paddedWidth - col + 1);
                    % Decode non-split blocks (I- or P-frame)
                    if I_frame_flag == 0 % P-frame
                        % update the ref index differentiation
                        diff_ref_index = MDiff_line_array(5);
                        best_ref_index = diff_ref_index + previous_ref_index;
                        previous_ref_index = best_ref_index;
                        referenceFrame_de = referenceFrames_de{best_ref_index};

                        % P-frame: Motion vector-based prediction
                        diff_mv = MDiff_line_array(3:4); % Extract motion vector difference
                        bestMatch = diff_mv + previous_mv; % Compute absolute motion vector
                        previous_mv = bestMatch;

                        if FMEEnable
                            int_dy = floor(bestMatch(2) / 2);
                            int_dx = floor(bestMatch(1) / 2);
                            frac_dy = mod(bestMatch(2), 2) / 2;
                            frac_dx = mod(bestMatch(1), 2) / 2;
    
                            refRow = row + int_dy;
                            refCol = col + int_dx;
                            refRow = max(1, min(refRow, paddedHeight - blockSize + 1));
                            refCol = max(1, min(refCol, paddedWidth - blockSize + 1));
    
                            predictorBlock = performReverseInterpolation(referenceFrame_de, refRow, refCol, frac_dy, frac_dx, blockSize);
                        else
                            % Use the motion vector to get the predictor block from the reference frame
                            refRow = row + bestMatch(2); % row + yOffset; 
                            refCol = col + bestMatch(1); % col + xOffset;
                            
                            % Check for boundary in reference frame
                            if refRow < 1 || refCol < 1 || refRow + blockSize - 1 > paddedHeight || refCol + blockSize - 1 > paddedWidth
                                predictorBlock = 128 * ones(block_height, block_width);
                            else
                                predictorBlock = referenceFrame_de(refRow:refRow+block_height-1, refCol:refCol+block_width-1);
                            end
                        end
                    else
                        % I-frame: Intra-prediction mode
                        % Extract mode difference and compute the mode
                        diff_mode = MDiff_line_array(1);
                        mode = diff_mode + previous_mode;
                        previous_mode = mode;

                        referenceFrames_de = [];

                        % Perform intra-prediction based on mode
                        if mode == 0
                            % Horizontal prediction
                            if row == 1
                                horizontal_pred = 128 * ones(block_height, block_width);
                            else
                                horizontal_pred = repmat(decodedFrame(row-1, col:col+block_width-1), block_height, 1);
                            end
                            predictorBlock = horizontal_pred;
                        else
                            % Vertical prediction
                            if col == 1
                                vertical_pred = 128 * ones(block_height, block_width);
                            else
                                vertical_pred = repmat(decodedFrame(row:row+block_height-1, col-1), 1, block_width);
                            end
                            predictorBlock = vertical_pred;
                        end
                    end
                    % Handle residual decoding from QTC stream
                    QTC_Line = fgetl(QTC_Coeff_file); % Read encoded residual
                    encoded_rle = A1_Q4_expGolombDecode(QTC_Line);
                    coeffs_scanned = A1_Q4_rleDecode(encoded_rle, blockSize);
                    residual_block_encoded = A1_Q4_inverseSScan(coeffs_scanned, blockSize, blockSize);
                    residualBlock_de = A3_idctAfterDequantizeBlock(residual_block_encoded, Q_Matrix);

                    % Reconstruct the block
                    reconstructed_block = double(predictorBlock) + double(residualBlock_de);
                    reconstructed_block = max(min(reconstructed_block, 255), 0);

                    % Place the reconstructed block into the decoded frame
                    decodedFrame(row:row+blockSize-1, col:col+blockSize-1) = reconstructed_block;
                end
            end
        end

        % Update the reference frame
        referenceFrames_de = A3_updateFIFObuffer(referenceFrames_de, nRefFrames, decodedFrame);
        
        % Remove padding by cropping the frame
        unpaddedDecodeFrame = decodedFrame(1:height, 1:width);
        % Write the decoded frame to the output file
        fwrite(vid_decoded, unpaddedDecodeFrame', 'uint8');

        % Compute Mean Squared Error (MSE)
        reconstructedFrame_comp = fread(vid_reconstructed, [width, height], 'uint8')';
        mse = mean((double(reconstructedFrame_comp(:)) - double(unpaddedDecodeFrame(:))).^2);
        % Compute PSNR
        if mse == 0
            psnrValues_verify(frameIdx) = Inf;  % Perfect reconstruction
        else
            psnrValues_verify(frameIdx) = 10 * log10(255^2 / mse);
        end
    end
    
    % Close the files
    fclose(splits);
    fclose(vid_decoded);
    fclose(vid_reconstructed);
    fclose(MDiff_file);
    fclose(QTC_Coeff_file);
end


function predictedBlock = performReverseInterpolation(refFrame, ref_row, ref_col, ...
    frac_dy, frac_dx, blockSize)

    [ref_height, ref_width] = size(refFrame);
    predictedBlock = zeros(blockSize, blockSize, 'uint8');

    for r = 1:blockSize
        for c = 1:blockSize
            y = ref_row + r - 1;
            x = ref_col + c - 1;

            top_left = refFrame(y, x);
            top_right = refFrame(y, min(x + 1, ref_width));
            bottom_left = refFrame(min(y + 1, ref_height), x);
            bottom_right = refFrame(min(y + 1, ref_height), min(x + 1, ref_width));

            if frac_dy == 0 && frac_dx == 0
                predictedPixel = top_left;
            else
                top_interp = top_left * (1 - frac_dx) + top_right * frac_dx;
                bottom_interp = bottom_left * (1 - frac_dx) + bottom_right * frac_dx;
                predictedPixel = top_interp * (1 - frac_dy) + bottom_interp * frac_dy;
            end
            predictedBlock(r, c) = uint8(predictedPixel);
        end
    end
end