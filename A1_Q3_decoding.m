function psnrValues_verify = A1_Q3_decoding(nframes, blockSize, paddedWidth, paddedHeight, height, width, QP)
    % Decoder for Y-only frames using motion vectors and residuals
    psnrValues_verify = zeros(1, nframes);  % Store PSNR for each frame, decode vs reconstructed

    % Initialize the hypothetical reference frame for the first frame (all values = 128)
    referenceFrame_de = 128 * ones(paddedHeight, paddedWidth, 'uint8');

    % Open the necessary files
    splits = fopen('splits.txt', 'w');

    vid_decoded = fopen('decoded_vid.yuv', 'w'); % Output decoded video
    vid_reconstructed = fopen('reconstructed_vid.yuv', 'r'); % Input reconstructed video for PSNR computation
    MDiff_file = fopen("MDiff.txt", 'r'); % Motion vector and mode file
    QTC_Coeff_file = fopen("QTC_Coeff.txt", 'r'); % Residual coefficients file

    % Loop through each frame
    for frameIdx = 1:nframes
        % Initialize the decoded frame
        decodedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
        Q_Matrix = A2_Q2_generateQMatrix(blockSize, QP); % Generate quantization matrix

        % Initialize predictors for motion vectors and modes
        previous_mv = [0, 0];
        previous_mode = 0;

        % Loop over blocks in raster order
        for row = 1:blockSize:paddedHeight
            for col = 1:blockSize:paddedWidth
                % Skip boundary blocks exceeding the frame size
                if row + blockSize - 1 > paddedHeight || col + blockSize - 1 > paddedWidth
                    continue;
                end
                MDiff_line = fgetl(MDiff_file);

                if isempty(MDiff_line)
                    error('Unexpected end of MDiff file.');
                end
                if ~(ischar(MDiff_line) || isstring(MDiff_line))
                    disp('Error: MDiff_line is not a valid string.');
                    error('MDiff_line must be a character array or string.');
                end

                % Check if MDiff_line is already in strtrim format
                if ~strcmp(MDiff_line, strtrim(MDiff_line))
                    disp('Error: MDiff_line contains leading or trailing whitespace.');
                end

                % Read the MDiff line and decode it
                string = strtrim(MDiff_line);
                MDiff_line_array = A1_Q4_expGolombDecode(MDiff_line);

                % Extract the split flag
                split_flag = MDiff_line_array(1);
                I_frame_flag = MDiff_line_array(2); % 0 = P-frame, 1 = I-frame   
                fprintf(splits, '%d\n', split_flag);
                if split_flag == 1
                    % Handle split blocks inline (logic from A2_Q2_decodeSplitBlock)
                    splitSize = blockSize / 2;
                    % Define relative offsets for Z-order traversal within the parent block
                    relative_offsets = [0, 0; 0, splitSize; splitSize, 0; splitSize, splitSize];
                    sub_Q_Matrix = A2_Q2_generateQMatrix(splitSize, QP); % Pre-compute Q matrix for sub-blocks
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
                            diff_mv = MDiff_line_array(1:2); % Read motion vector differences
                            bestMatch = diff_mv + previous_mv; % Compute absolute motion vector
                            previous_mv = bestMatch;

                            % Reference block position
                            refRow = absSubRow + bestMatch(2);
                            refCol = absSubCol + bestMatch(1);
                
                            if refRow < 1 || refCol < 1 || refRow + splitSize - 1 > paddedHeight || refCol + splitSize - 1 > paddedWidth
                                predictorBlock = 128 * ones(block_height, block_width); % Fallback to neutral predictor
                            else
                                predictorBlock = referenceFrame_de(refRow:refRow+block_height-1, refCol:refCol+block_width-1);
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
                    % Decode non-split blocks (I- or P-frame)
                    
                    if I_frame_flag == 0
                        % P-frame: Motion vector-based prediction
                        diff_mv = MDiff_line_array(3:4); % Extract motion vector difference
                        bestMatch = diff_mv + previous_mv; % Compute absolute motion vector
                        previous_mv = bestMatch;

                        % Compute reference block position
                        refRow = row + bestMatch(2);
                        refCol = col + bestMatch(1);

                        % Check reference block boundaries
                        if refRow < 1 || refCol < 1 || refRow + blockSize - 1 > paddedHeight || refCol + blockSize - 1 > paddedWidth
                            continue; % Skip out-of-bound blocks
                        end

                        % Extract predictor block from reference frame
                        predictorBlock = referenceFrame_de(refRow:refRow+blockSize-1, refCol:refCol+blockSize-1);
                    else
                        % I-frame: Intra-prediction mode
                        block_height = min(blockSize, paddedHeight - row + 1);
                        block_width = min(blockSize, paddedWidth - col + 1);

                        % Extract mode difference and compute the mode
                        diff_mode = MDiff_line_array(3);
                        mode = diff_mode + previous_mode;
                        previous_mode = mode;

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

                    residualBlock_de = A2_Q2_idctAfterDequantizeBlock(residual_block_encoded, Q_Matrix);

                    % Reconstruct the block
                    reconstructed_block = double(predictorBlock) + double(residualBlock_de);
                    reconstructed_block = max(min(reconstructed_block, 255), 0);

                    % Place the reconstructed block into the decoded frame
                    decodedFrame(row:row+blockSize-1, col:col+blockSize-1) = reconstructed_block;
                end
            end
        end

        % Update the reference frame
        referenceFrame_de = decodedFrame;

        % Remove padding and write decoded frame to the output file
        unpaddedDecodeFrame = decodedFrame(1:height, 1:width);
        fwrite(vid_decoded, unpaddedDecodeFrame', 'uint8');

        % Compute PSNR for the frame
        reconstructedFrame_comp = fread(vid_reconstructed, [width, height], 'uint8')';
        mse = mean((double(reconstructedFrame_comp(:)) - double(unpaddedDecodeFrame(:))).^2);
        if mse == 0
            psnrValues_verify(frameIdx) = Inf; % Perfect reconstruction
        else
            psnrValues_verify(frameIdx) = 10 * log10(255^2 / mse);
        end
    end

    % Close all files
    fclose(splits);
    fclose(vid_decoded);
    fclose(vid_reconstructed);
    fclose(MDiff_file);
    fclose(QTC_Coeff_file);
end