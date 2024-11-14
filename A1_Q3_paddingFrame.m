function paddedY = A1_Q3_paddingFrame(frame, height, width, blockSize)
    % calculate extra width and height that need to be padded
    padWidth = mod(width, blockSize);
    padHeight = mod(height, blockSize);

    % check padding needed
    % Apply padding if necessary
    if padWidth ~= 0
        padWidth = blockSize - padWidth;
    end
    if padHeight ~= 0
        padHeight = blockSize - padHeight;
    end

    paddedY = padarray(frame, [padHeight, padWidth], 128, 'post');
    
end 