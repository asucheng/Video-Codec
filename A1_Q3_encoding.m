function [psnrValues] = A1_Q3_encoding(nframes, paddedWidth, paddedHeight, blockSize, height, width, searchRange, n, I_period, QP, VBSenable)
    % Open video files
    vid_out_Y_pad = fopen('y_only_padded.yuv', 'r');
    vid_reconstructed = fopen('reconstructed_vid.yuv', 'w');
    predicted_vid = fopen('predicted.yuv', 'w');
    MDiff_stream = fopen('MDiff.txt', 'w');
    QTC_stream = fopen('QTC_Coeff.txt', 'w');
    psnrValues = zeros(1, nframes); % Store PSNR for each frame
    referenceFrame = 128 * ones(paddedHeight, paddedWidth, 'uint8'); % Initial reference frame

    % Loop through frames
    for frameIdx = 1:nframes
        currentFrame = fread(vid_out_Y_pad, [paddedWidth, paddedHeight], 'uint8')';
        
        if mod(frameIdx - 1, I_period) == 0
            fprintf('Processing I-frame %d\n', frameIdx);
            [predictedFrame, reconstructedFrame] = A2_Q2_intraPredictForIFrame(currentFrame, blockSize, QP, MDiff_stream, QTC_stream, VBSenable);
        else
            fprintf('Processing P-frame %d\n', frameIdx);
            [predictedFrame, reconstructedFrame] = A2_Q2_interPredictForPFrame(referenceFrame, currentFrame, searchRange, blockSize, paddedHeight, paddedWidth, n, QP, MDiff_stream, QTC_stream, VBSenable);
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