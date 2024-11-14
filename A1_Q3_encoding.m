function [psnrValues, maes] = A1_Q3_encoding(nframes, paddedWidth, paddedHeight, blockSize, height, width, searchRange, n, I_period, QP_values)
    % Open video files
    vid_out_Y_pad = fopen('y_only_padded.yuv', 'r');
    vid_out_Y = fopen('y_only.yuv', 'r');
    vid_reconstructed = fopen('reconstructed_vid.yuv', 'w');
    predicted_vid = fopen('predicted.yuv','w');

    mvFile = fopen('motion_vectors.txt', 'w');
    residualFile = fopen('residuals.txt', 'w');
    
    % Store PSNR for each frame, reconstructed vs original
    psnrValues = zeros(1, nframes); 
    maes = zeros(1, nframes); 

    % initialize a reference frame for first frame
    referenceFrame = 128 * ones(paddedHeight, paddedWidth, 'uint8');

    % loop throught each frame
    for frameIdx = 1:nframes
        % Read current frame
        currentFrame = fread(vid_out_Y_pad, [paddedWidth, paddedHeight], 'uint8')';
        %currentIFrame = fread(vid_out_Y, [width, height], 'uint8')';

        if mod(frameIdx - 1, I_period) == 0
            fprintf('Dealing with frame %d, I-Period\n', frameIdx);
            [predictedFrame, reconstructedFrame, MDiff_stream, QTC_stream] = A1_Q4_intraPredictForIFrame(currentFrame, blockSize, QP_values); 
            type = 'I';
        else
            fprintf('Dealing with frame %d, F-Period\n', frameIdx);
            [predictedFrame, reconstructedFrame, MDiff_stream, QTC_stream] = A1_Q3_interPredictPFrame(referenceFrame, currentFrame, searchRange, blockSize, paddedHeight, paddedWidth, n, QP_values);
            type = 'P';
        end

        % Update the reference frame
        referenceFrame = reconstructedFrame;

        % Remove padding by cropping the frame
        % Write predicted frame
        unpaddedpredictedFrame = predictedFrame(1:height, 1:width);
        fwrite(predicted_vid, unpaddedpredictedFrame', 'uint8');

        % Remove padding by cropping the frame
        % Write reconstructed frame
        unpaddedReconstructedFrame = reconstructedFrame(1:height, 1:width);
        fwrite(vid_reconstructed, unpaddedReconstructedFrame', 'uint8');

        % Compute Mean Squared Error (MSE)
        unpaddedCurrentFrame = currentFrame(1:height, 1:width);
        mse = mean((double(unpaddedCurrentFrame(:)) - double(unpaddedReconstructedFrame(:))).^2);
        
        % Compute PSNR
        if mse == 0
            psnrValues(frameIdx) = Inf;  % Perfect reconstruction
        else
            psnrValues(frameIdx) = 10 * log10(255^2 / mse);
        end
    end

    % Close files
    fclose(vid_out_Y_pad);
    fclose(vid_out_Y);
    fclose(vid_reconstructed);
    fclose(predicted_vid);

    fclose(mvFile);
    fclose(residualFile);
end