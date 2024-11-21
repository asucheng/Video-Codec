function display_y_only(yuvFile)
    % Function to read and display Y-only frames from a YUV file.
    %
    % INPUT:
    %   yuvFile - Path to the Y-only YUV file (string)

    % Set video parameters (you may adjust these as needed)
    width = 352;  % Frame width
    height = 288; % Frame height
    nframes = 30; % Total number of frames to read

    % Open the Y-only file for reading
    fileID = fopen(yuvFile, 'rb');
    if fileID == -1
        error('Error opening Y-only file.');
    end

    % Loop through each frame and display it
    for frameIdx = 1:nframes
        % Read one frame (Y-component only)
        Y = fread(fileID, width * height, 'uint8');
        
        if isempty(Y)
            disp('End of file reached.');
            break;
        end

        % Reshape the frame into 2D matrix (height x width)
        Y_frame = reshape(Y, [width, height])';

        % Display the frame
        figure(1); % Reuse the same figure window
        imshow(Y_frame, [0 255]); % Display as grayscale image
        title(sprintf('Y-Only Frame %d', frameIdx));
        colormap('gray');
        axis on;

        % Pause briefly to visualize frames (adjust as needed)
        pause(0.6);
    end

    % Close the file
    fclose(fileID);
    disp('Finished displaying all frames.');
end