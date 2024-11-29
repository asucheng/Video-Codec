function [psnrValues] = Encoding(filename_prefix, nframes, paddedWidth, paddedHeight, ...
    blockSize, height, width, searchRange, n, I_period, QP_values, ...
    nRefFrames, VBSEnable, MRFoverlay, FMEEnable, FastME)

    % Open video files
    vid_out_Y_pad = fopen('y_only_padded.yuv', 'r');
    vid_reconstructed = fopen(strcat(filename_prefix, 'reconstructed_vid.yuv'), 'w');
    predicted_vid = fopen(strcat(filename_prefix, 'predicted.yuv'),'w');

    MDiff_stream = fopen(strcat(filename_prefix, 'MDiff.txt'), 'w');
    QTC_stream = fopen(strcat(filename_prefix, 'QTC_Coeff.txt'), 'w');
    MVPDiff_stream = fopen(strcat(filename_prefix, 'MVPDiff.txt'), 'w');
    
    psnrValues = zeros(1, nframes); % Store PSNR for each frame
    % Loop through frames
    for frameIdx = 1:nframes
        currentFrame = fread(vid_out_Y_pad, [paddedWidth, paddedHeight], 'uint8')';

        if mod(frameIdx - 1, I_period) == 0
            fprintf('Processing I-frame %d\n', frameIdx);
            [predictedFrame, reconstructedFrame] = A2_intraPredictForIFrame(currentFrame, blockSize, ...
                QP_values, MDiff_stream, MVPDiff_stream, QTC_stream, ...
                VBSEnable, FMEEnable, FastME); 
            % Clear the reference frames on I-frame
            reference_frames = [];
        else
            fprintf('Processing P-frame %d\n', frameIdx);
            % FMEEnable only works for P Frame
            [predictedFrame, reconstructedFrame] = A2_interPredictForPFrame(reference_frames, currentFrame, searchRange, ...
                blockSize, paddedHeight, paddedWidth, n, QP_values, MDiff_stream, MVPDiff_stream, QTC_stream, ...
                nRefFrames, frameIdx, VBSEnable, MRFoverlay, FMEEnable, FastME);
        end

        reference_frames = A2_updateFIFObuffer(reference_frames, nRefFrames, reconstructedFrame); % Update reference frame

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
    fclose(MVPDiff_stream);
end