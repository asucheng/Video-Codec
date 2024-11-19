function A1_Q3_extract_YOnly(originalVid, og_height, og_width, nframes)
    % output directory if save a frames
    % Y_only_dir = 'y_only';
    % if ~exist(Y_only_dir, 'dir')
    %     mkdir(Y_only_dir);
    % end
    
    vid = fopen(originalVid, 'r');
    vid_out_Y_only = fopen('y_only.yuv', 'w');
    
    % Iterate through the frames, carve out the U and V
    for frame = 1:nframes
        % read YUV frame from the video file
        Y = fread(vid, og_width * og_height, 'uint8');  
    
        % save as a video sequence
        fwrite(vid_out_Y_only, Y, 'uint8'); 
    
        % save as individual frames
        % outputFileName = fullfile(Y_only_dir, sprintf('frame%03d.yuv', frame));
        % outputFile = fopen(outputFileName, 'w');
        % fwrite(outputFile, Y, 'uint8');
        % fclose(outputFile);
    
        % skip UV for each frame
        fseek(vid, (og_width/2) * (og_height/2) * 2, 'cof');
    end
    
    fclose(vid);
    fclose(vid_out_Y_only);
    
    disp('Y-only files Y!');
end
