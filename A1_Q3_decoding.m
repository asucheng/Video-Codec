function psnrValues_verify = A1_Q3_decoding(filename_prefix, nframes, blockSize, paddedWidth, paddedHeight, ...
    height, width, QP, FMEEnable, FastME)

    % Decoder for Y-only frames using motion vectors and residuals
    psnrValues_verify = zeros(1, nframes);  % Store PSNR for each frame, decode vs reconstructed
    
    % Initialize the hypothetical reference frame for the first frame (all values = 128)
    referenceFrame_de = 128 * ones(paddedHeight, paddedWidth, 'uint8');
    
    % Open the necessary files
    vid_decoded = fopen(strcat(filename_prefix, 'decoded_vid.yuv'), 'w');
    % mvFile = fopen('motion_vectors.txt', 'r');
    % residualFile = fopen('residuals.txt', 'r');
    vid_reconstructed = fopen(strcat(filename_prefix, 'reconstructed_vid.yuv'), 'r');

    MDiff_file = fopen(strcat(filename_prefix, "MDiff.txt"), 'r');
    MVPDiff_file = fopen(strcat(filename_prefix, "MVPDiff.txt"), 'r');
    QTC_Coeff_file = fopen(strcat(filename_prefix, "QTC_Coeff.txt"), 'r');

    if FastME
        Diff_line = MVPDiff_file;
    else
        Diff_line = MDiff_file;
    end
    
    % Loop through frames
    for frameIdx = 1:nframes
        % Initialize the decoded frame
        decodedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
        Q_Matrix = A1_Q4_generateQMatrix(blockSize, QP);

        previous_mv = [0, 0];
        previous_mode = 0;
        
        % Loop over blocks in raster order
        for row = 1:blockSize:paddedHeight
            for col = 1:blockSize:paddedWidth
                % Check for boundary issues (skip if block exceeds frame size)
                if row + blockSize - 1 > paddedHeight || col + blockSize - 1 > paddedWidth
                    continue;
                end

                % read the a line of MDiff file----------------------------
                MDiff_line = fgetl(Diff_line);
                MDiff_line_array = A1_Q4_expGolombDecode(MDiff_line);

                if MDiff_line_array(1) == 0
                    % decoder P frame
                    % update the motion vector differentiation
                    diff_mv = MDiff_line_array(end-1:end);
                    bestMatch = diff_mv + previous_mv;
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
                            continue;  % Skip if the predictor block goes out of bounds
                        end
                        
                        % Extract the predictor block from the reference frame
                        predictorBlock = referenceFrame_de(refRow:refRow+blockSize-1, refCol:refCol+blockSize-1);
                    end

                    % handle QTC file------------------------------------------
                    % residualBlock_de = A1_Q3_readResidualFile(QTC_Coeff_file, blockSize);
                    QTC_Line = fgetl(QTC_Coeff_file);  % Read the residual block line
                    encoded_rle = A1_Q4_expGolombDecode(QTC_Line);
                    coeffs_scanned = A1_Q4_rleDecode(encoded_rle, blockSize);
                    residual_block_encoded = A1_Q4_inverseSScan(coeffs_scanned, blockSize, blockSize);
                    residualBlock_de = A1_Q4_idctAfterDequantizeBlock(residual_block_encoded, Q_Matrix);
                    
                    % Add the residual block to the predictor block to get the decoded block
                    decodedBlock = double(predictorBlock) + double(residualBlock_de);
                    decodedBlock = max(min(decodedBlock, 255), 0);
                    
                    % Place the decoded block into the decoded frame
                    decodedFrame(row:row+blockSize-1, col:col+blockSize-1) = decodedBlock;
                else
                    % decode I frame
                    block_height = min(blockSize, paddedHeight - row + 1);
                    block_width = min(blockSize, paddedWidth - col + 1);
                    diff_mode = MDiff_line_array(3);

                    mode = diff_mode + previous_mode;
                    previous_mode = mode;

                    if mode == 0
                        if row == 1
                            horizontal_pred = 128 * ones(block_height, block_width);
                        else
                            horizontal_pred = repmat(decodedFrame(row-1, col:col+block_width-1), block_height, 1);
                        end
                        predictorBlock = horizontal_pred;
                    else
                        if col == 1
                            vertical_pred = 128 * ones(block_height, block_width);
                        else
                            vertical_pred = repmat(decodedFrame(row:row+block_height-1, col-1), 1, block_width);
                        end
                        predictorBlock = vertical_pred;
                    end

                    % handle QTC file------------------------------------------
                    QTC_Line = fgetl(QTC_Coeff_file);  % Read the residual block line
                    encoded_rle = A1_Q4_expGolombDecode(QTC_Line);
                    coeffs_scanned = A1_Q4_rleDecode(encoded_rle, blockSize);
                    residual_block_encoded = A1_Q4_inverseSScan(coeffs_scanned, blockSize, blockSize);
                    residualBlock_de = A1_Q4_idctAfterDequantizeBlock(residual_block_encoded, Q_Matrix);
                    
                    % Add the residual block to the predictor block to get the decoded block
                    reconstructed_block = double(predictorBlock) + double(residualBlock_de);
                    reconstructed_block = max(min(reconstructed_block, 255), 0);
                    
                    % Place the decoded block into the decoded frame
                    decodedFrame(row:row+blockSize-1, col:col+blockSize-1) = reconstructed_block;
                end
            end
        end

        % Update the reference frame for the next iteration
        referenceFrame_de = decodedFrame;
        
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
    fclose(vid_decoded);
    fclose(vid_reconstructed);
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