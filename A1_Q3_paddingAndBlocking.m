% Run the section, then double click "blocked_frame_cell_array" 
% in workspace to see the blocked frames
function [paddedWidth, paddedHeight] = A1_Q3_paddingAndBlocking(original_Width, original_Height, blockSize, block_array, nframes)
    %used to double check operation 
    vid_Y = fopen('y_only.yuv', 'r');
    vid_out_Y_pad = fopen('y_only_padded.yuv', 'w');
    
    % Iterate over all Y only frames
    for frame = 1:nframes
        % read frame into 352 x 288 mat, matlab read col by col, therefore
        % transpose needed, matlab mat dimension is in height x width
        Y_frame = fread(vid_Y, [original_Width, original_Height], 'uint8')';
    
        %padding each frame 
        paddedFrame = A1_Q3_paddingFrame(Y_frame, original_Height, original_Width, blockSize);
    
        % to save as video sequence to check result, transpose again
        frameToWrite = paddedFrame';
        fwrite(vid_out_Y_pad, frameToWrite, 'uint8'); 
    
        % split the frame into matrix of matrx store it as cell type, no need
        % to use "output mat transpose"
        blocked_frame = mat2cell(paddedFrame, repmat(blockSize, 1, size(paddedFrame, 1) / blockSize), ...
                          repmat(blockSize, 1, size(paddedFrame, 2) / blockSize));
    
        %save each blocked frame into a cell of 1x300
        block_array{1, frame} = blocked_frame;
    end
    
    % get new width and height for padded frame
    paddedWidth = size(paddedFrame, 2);
    paddedHeight = size(paddedFrame, 1);
    
    %used to double check operation 
    fclose(vid_Y);
    fclose(vid_out_Y_pad);
    % To view the padded video, need to find out the new resolution for video to play properly, no longer CIF.
    
    % Display a message when block splitting and padding are complete
    disp('Y frame is split to blocks and padded as needed!');
end