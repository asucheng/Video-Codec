function psnrValues_verify = A1_Q3_decoding(nframes, blockSize, paddedWidth, paddedHeight, height, width)    
    % Decoder for Y-only frames using motion vectors and residuals
    psnrValues_verify = zeros(1, nframes);  % Store PSNR for each frame, decode vs reconstructed
    
    % Initialize the hypothetical reference frame for the first frame (all values = 128)
    referenceFrame_de = 128 * ones(paddedHeight, paddedWidth, 'uint8');
    
    % Open the necessary files
    vid_decoded = fopen('decoded_vid.yuv', 'w');
    mvFile = fopen('motion_vectors.txt', 'r');
    residualFile = fopen('residuals.txt', 'r');
    vid_reconstructed = fopen('reconstructed_vid.yuv', 'r');
    
    % Loop through frames
    for frameIdx = 1:nframes
        % Initialize the decoded frame
        decodedFrame = zeros(paddedHeight, paddedWidth, 'uint8');
        
        % Loop over blocks in raster order
        for row = 1:blockSize:paddedHeight
            for col = 1:blockSize:paddedWidth
                % Check for boundary issues (skip if block exceeds frame size)
                if row + blockSize - 1 > paddedHeight || col + blockSize - 1 > paddedWidth
                    continue;
                end
                
                [refRow, refCol] = A1_Q3_readMVandSelectBlock(mvFile, row, col);
                
                % Check for boundary in reference frame
                if refRow < 1 || refCol < 1 || refRow + blockSize - 1 > paddedHeight || refCol + blockSize - 1 > paddedWidth
                    continue;  % Skip if the predictor block goes out of bounds
                end
                
                % Extract the predictor block from the reference frame
                predictorBlock = referenceFrame_de(refRow:refRow+blockSize-1, refCol:refCol+blockSize-1);
                
                residualBlock_de = A1_Q3_readResidualFile(residualFile, blockSize);
                
                % Add the residual block to the predictor block to get the decoded block
                decodedBlock = double(predictorBlock) + double(residualBlock_de);
                
                % Place the decoded block into the decoded frame
                decodedFrame(row:row+blockSize-1, col:col+blockSize-1) = decodedBlock;
            end
        end
        
        % Remove padding by cropping the frame
        unpaddedDecodeFrame = decodedFrame(1:height, 1:width);
        % Write the decoded frame to the output file
        fwrite(vid_decoded, unpaddedDecodeFrame', 'uint8');
        
        % Update the reference frame for the next iteration
        referenceFrame_de = decodedFrame;
    
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
    fclose(mvFile);
    fclose(residualFile);
    fclose(vid_reconstructed);
end