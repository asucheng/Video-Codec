function psnrValues_verify = A1_Q3_decoding(nframes, blockSize, paddedWidth, paddedHeight, height, width, QP)    
    % Decoder for Y-only frames using motion vectors and residuals
    psnrValues_verify = zeros(1, nframes);  % Store PSNR for each frame, decode vs reconstructed
    
    % Initialize the hypothetical reference frame for the first frame (all values = 128)
    referenceFrame_de = 128 * ones(paddedHeight, paddedWidth, 'uint8');
    
    % Open the necessary files
    vid_decoded = fopen('decoded_vid.yuv', 'w');
    vid_reconstructed = fopen('reconstructed_vid.yuv', 'r');

    MDiff_file = fopen("MDiff.txt", 'r');
    QTC_Coeff_file = fopen("QTC_Coeff.txt", 'r');
    
    % Loop through frames
    for frameIdx = 1:nframes
        % Initialize the decoded frame
        decodedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
        lambda = A2_Q2_getLambda(QP);
        Q_Matrix = A2_Q2_generateQMatrix(blockSize, QP);

        previous_mv = [0, 0];
        previous_mode = 0;
        
        % Loop over blocks in raster order
        for row = 1:blockSize:paddedHeight
            for col = 1:blockSize:paddedWidth
                % Check for boundary issues (skip if block exceeds frame size)
                if row + blockSize - 1 > paddedHeight || col + blockSize - 1 > paddedWidth
                    continue;
                end
                
                % Read the split flag from MDiff_stream
                MDiff_line = fgetl(MDiff_file);
                MDiff_line_array = A1_Q4_expGolombDecode(MDiff_line);

                split_flag = MDiff_line_array(1);
                
                if split_flag == 1
                    % Handle split blocks recursively
                    
                    sub_block = A2_Q2_decodeSplitBlock(referenceFrame_de, QTC_Coeff_file, MDiff_file, ...
                        row, col, blockSize, paddedHeight, paddedWidth, QP);
                    decodedFrame(row:row+blockSize-1, col:col+blockSize-1) = sub_block;
                else
                    % Decode non-split blocks (I or P frame)
                    if MDiff_line_array(2) == 0
                        % P-frame logic
                        diff_mv = MDiff_line_array(end-1:end);
                        bestMatch = diff_mv + previous_mv;
                        previous_mv = bestMatch;

                        % Get reference block
                        refRow = row + bestMatch(2);
                        refCol = col + bestMatch(1);

                        % Check for reference block boundary
                        if refRow < 1 || refCol < 1 || refRow + blockSize - 1 > paddedHeight || refCol + blockSize - 1 > paddedWidth
                            continue;  % Skip out-of-bound blocks
                        end
                        predictorBlock = referenceFrame_de(refRow:refRow+blockSize-1, refCol:refCol+blockSize-1);
                    else
                        % I-frame logic
                        block_height = min(blockSize, paddedHeight - row + 1);
                        block_width = min(blockSize, paddedWidth - col + 1);
                        diff_mode = MDiff_line_array(2);
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
                    end

                    % Handle residuals from QTC_stream
                    QTC_Line = fgetl(QTC_Coeff_file);
                    encoded_rle = A1_Q4_expGolombDecode(QTC_Line);
                    coeffs_scanned = A1_Q4_rleDecode(encoded_rle, blockSize);
                    residual_block_encoded = A1_Q4_inverseSScan(coeffs_scanned, blockSize, blockSize);
                    residualBlock_de = A2_Q2_idctAfterDequantizeBlock(residual_block_encoded, Q_Matrix);
                    
                    % Reconstruct block
                    reconstructed_block = double(predictorBlock) + double(residualBlock_de);
                    reconstructed_block = max(min(reconstructed_block, 255), 0);
                    
                    % Place the reconstructed block
                    decodedFrame(row:row+blockSize-1, col:col+blockSize-1) = reconstructed_block;
                end
            end
        end

        % Update reference frame
        referenceFrame_de = decodedFrame;
        
        % Write decoded frame to file
        unpaddedDecodeFrame = decodedFrame(1:height, 1:width);
        fwrite(vid_decoded, unpaddedDecodeFrame', 'uint8');

        % Compute PSNR
        reconstructedFrame_comp = fread(vid_reconstructed, [width, height], 'uint8')';
        mse = mean((double(reconstructedFrame_comp(:)) - double(unpaddedDecodeFrame(:))).^2);
        if mse == 0
            psnrValues_verify(frameIdx) = Inf;
        else
            psnrValues_verify(frameIdx) = 10 * log10(255^2 / mse);
        end
    end

    % Close files
    fclose(vid_decoded);
    fclose(vid_reconstructed);
end