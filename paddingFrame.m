function paddedY = paddingFrame(frame, height, width, blockSize)
    % calculate extra width and height that need to be padded
    padWidth = mod(width, blockSize);
    padHeight = mod(height, blockSize);

    % check padding needed
    % Apply padding if necessary
    if padWidth ~= 0
        padWidth_N = blockSize - padWidth;
    end
    if padHeight ~= 0
        padHeight_N = blockSize - padHeight;
    end

    paddedY = padarray(frame, [padHeight_N, padWidth_N], 128, 'post');
    
end 