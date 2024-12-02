function [psnrValues, frameBitUsage, rowBlockCount] = A3_Encoding(filename_prefix, nframes, paddedWidth, paddedHeight, ...
    blockSize, height, width, searchRange, n, I_period, QP_value, ...
    nRefFrames, VBSEnable, MRFoverlay, FMEEnable, FastME, ...
    RCflag, targetBR, iFrameAverageBitCount, pFrameAverageBitCount, ...
    blockSizeCollections, QPValueCollections)

    % Open video files
    vid_out_Y_pad = fopen('y_only_padded.yuv', 'r');
    vid_reconstructed = fopen(strcat(filename_prefix, 'reconstructed_vid.yuv'), 'w');
    predicted_vid = fopen(strcat(filename_prefix, 'predicted.yuv'),'w');

    MDiff_stream = fopen(strcat(filename_prefix, 'MDiff.txt'), 'w');
    QTC_stream = fopen(strcat(filename_prefix, 'QTC_Coeff.txt'), 'w');
    MVPDiff_stream = fopen(strcat(filename_prefix, 'MVPDiff.txt'), 'w');
    
    psnrValues = zeros(1, nframes); % Store PSNR for each frame
    frameBitUsage = zeros(1, nframes); % Track actual bit usage for each frame
    rowBlockCount = zeros(1, nframes);

    if RCflag
        % Initialize total bit budget
        totalBitBudget = targetBR; % Total budget for all frames
        remainingBitBudget = totalBitBudget; % Remaining bit budget
    end

    for frameIdx = 1:nframes
        currentFrame = fread(vid_out_Y_pad, [paddedWidth, paddedHeight], 'uint8')';

        if RCflag
            if frameIdx <= nframes - 1
                budgetPerFrame = remainingBitBudget / (nframes - frameIdx + 1);
            else
                budgetPerFrame = remainingBitBudget;
            end

            if mod(frameIdx - 1, I_period) == 0
                fprintf('Processing I-frame %d with frame budget: %.2f\n', frameIdx, budgetPerFrame);
                [predictedFrame, reconstructedFrame, bitCountConsumed, rowBlockCountThisFrame] = A3_intraPredictForIFrame( ...
                    currentFrame, blockSize, QP_value, MDiff_stream, MVPDiff_stream, QTC_stream, ...
                    VBSEnable, FMEEnable, FastME, ...
                    RCflag, budgetPerFrame, iFrameAverageBitCount, ...
                    blockSizeCollections, QPValueCollections);
                    
                % Clear the reference frames on I-frame
                reference_frames = [];
            else
                fprintf('Processing P-frame %d with frame budget: %.2f\n', frameIdx, budgetPerFrame);
                [predictedFrame, reconstructedFrame, bitCountConsumed, rowBlockCountThisFrame] = A3_interPredictForPFrame( ...
                    reference_frames, currentFrame, searchRange, ...
                    blockSize, paddedHeight, paddedWidth, n, QP_value, MDiff_stream, MVPDiff_stream, QTC_stream, ...
                    nRefFrames, frameIdx, VBSEnable, MRFoverlay, FMEEnable, FastME, ...
                    RCflag, budgetPerFrame, pFrameAverageBitCount, ...
                    blockSizeCollections, QPValueCollections);
            end
            remainingBitBudget = remainingBitBudget - bitCountConsumed;
        else
            if mod(frameIdx - 1, I_period) == 0
                % fprintf('Processing I-frame %d\n', frameIdx);
                [predictedFrame, reconstructedFrame, bitCountConsumed, rowBlockCountThisFrame] = A3_intraPredictForIFrame(currentFrame, blockSize, ...
                    QP_value, MDiff_stream, MVPDiff_stream, QTC_stream, ...
                    VBSEnable, FMEEnable, FastME, false, ...
                    inf, [], [], []);
    
                % Clear the reference frames on I-frame
                reference_frames = [];
            else
                % fprintf('Processing P-frame %d\n', frameIdx);
                % FMEEnable only works for P Frame
                [predictedFrame, reconstructedFrame, bitCountConsumed, rowBlockCountThisFrame] = A3_interPredictForPFrame(reference_frames, currentFrame, searchRange, ...
                    blockSize, paddedHeight, paddedWidth, n, QP_value, MDiff_stream, MVPDiff_stream, QTC_stream, ...
                    nRefFrames, frameIdx, VBSEnable, MRFoverlay, FMEEnable, FastME, ...
                    false, inf, [], [], []);
            end

        end

        frameBitUsage(frameIdx) = bitCountConsumed;
        rowBlockCount(frameIdx) = rowBlockCountThisFrame;

        reference_frames = A3_updateFIFObuffer(reference_frames, nRefFrames, reconstructedFrame); % Update reference frame

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