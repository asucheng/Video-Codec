function [psnrValues] = A1_Q3_encoding(filename_prefix, nframes, paddedWidth, paddedHeight, ...
    blockSize, height, width, searchRange, n, ...
    I_period, QP_values, VBSenable, FMEEnable, FastME)

    % Open video files
    vid_out_Y_pad = fopen('y_only_padded.yuv', 'r');
    vid_out_Y = fopen('y_only.yuv', 'r');
    vid_reconstructed = fopen(strcat(filename_prefix, 'reconstructed_vid.yuv'), 'w');
    predicted_vid = fopen(strcat(filename_prefix, 'predicted.yuv'),'w');

    % mvFile = fopen('motion_vectors.txt', 'w');
    % residualFile = fopen('residuals.txt', 'w');

    MDiff_stream = fopen(strcat(filename_prefix, 'MDiff.txt'), 'w');
    QTC_stream = fopen(strcat(filename_prefix, 'QTC_Coeff.txt'), 'w');
    MVPDiff_stream = fopen(strcat(filename_prefix, 'MVPDiff.txt'), 'w');
    
    psnrValues = zeros(1, nframes); % Store PSNR for each frame
    referenceFrame = 128 * ones(paddedHeight, paddedWidth, 'uint8'); % Initial reference frame

    % Loop through frames
    for frameIdx = 1:nframes
        currentFrame = fread(vid_out_Y_pad, [paddedWidth, paddedHeight], 'uint8')';
        % currentIFrame = fread(vid_out_Y, [width, height], 'uint8')';

        if mod(frameIdx - 1, I_period) == 0
            % fprintf('Dealing with frame %d, I-Period\n', frameIdx);
            [predictedFrame, reconstructedFrame] = A2_Q34_intraPredictForIFrame(currentFrame, blockSize, QP_values, ...
                MDiff_stream, MVPDiff_stream, QTC_stream, FMEEnable, FastME); 
            type = 'I';
        else
            % fprintf('Dealing with frame %d, P-Period\n', frameIdx);
            % FMEEnable only works for P Frame
            [predictedFrame, reconstructedFrame] = A2_Q34_interPredictForPFrame(referenceFrame, currentFrame, searchRange, ...
                blockSize, paddedHeight, paddedWidth, n, QP_values, MDiff_stream, MVPDiff_stream, QTC_stream, ...
                FMEEnable, FastME);
            type = 'P';
        end

        referenceFrame = reconstructedFrame; % Update reference frame

        % Save predicted and reconstructed frames
        unpaddedPredictedFrame = predictedFrame(1:height, 1:width);
        fwrite(predicted_vid, unpaddedPredictedFrame', 'uint8');
        unpaddedReconstructedFrame = reconstructedFrame(1:height, 1:width);
        fwrite(vid_reconstructed, unpaddedReconstructedFrame', 'uint8');
        unpaddedCurrentFrame = currentFrame(1:height, 1:width);
        mse = mean((double(unpaddedCurrentFrame(:)) - double(unpaddedReconstructedFrame(:))).^2);
        
        % Compute PSNR
        psnrValues(frameIdx) = 10 * log10(255^2 / mse);
    end

    fclose(vid_out_Y_pad);
    fclose(vid_reconstructed);
    fclose(predicted_vid);
    fclose(MDiff_stream);
    fclose(QTC_stream);
end