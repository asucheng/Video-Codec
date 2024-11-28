function [interpolatedBlock] = A2_performInterpolation(refFrame, mv, row, col, blockSize)
    int_dy = floor(mv(1) / 2); % Integer vertical motion
    int_dx = floor(mv(2) / 2); % Integer horizontal motion
    frac_dy = mod(mv(1), 2) / 2; % Fractional vertical motion (0 or 0.5)
    frac_dx = mod(mv(2), 2) / 2; % Fractional horizontal motion (0 or 0.5)

    % Perform interpolation for fractional motion estimation.
    interpolatedBlock = zeros(blockSize, blockSize, 'uint8');

    [ref_height, ref_width] = size(refFrame);

    row = row + int_dy;
    col = col + int_dx;

    for r = 1:blockSize
        for c = 1:blockSize
            % Get fractional pixel location
            y = row + r - 1;
            x = col + c - 1;

            top_left_y = max(1, min(y, ref_height));
            top_left_x = max(1, min(x, ref_width));
            top_right_x = max(1, min(x + 1, ref_width));
            bottom_left_y = max(1, min(y + 1, ref_height));
            bottom_right_y = max(1, min(y + 1, ref_height));
            bottom_right_x = max(1, min(x + 1, ref_width));
    
            top_left = refFrame(top_left_y, top_left_x);
            top_right = refFrame(top_left_y, top_right_x);
            bottom_left = refFrame(bottom_left_y, top_left_x);
            bottom_right = refFrame(bottom_right_y, bottom_right_x);

            % Direct reference pixel (MV = [0, 0])
            if frac_dx == 0 && frac_dy == 0
                interpolatedBlock(r, c) = top_left;

            % Horizontal interpolation (MV = [0, 1])
            elseif frac_dy == 0 && frac_dx ~= 0
                interpolatedBlock(r, c) = uint8(top_left * (1 - frac_dx) + top_right * frac_dx);

            % Vertical interpolation (MV = [1, 0])
            elseif frac_dy ~= 0 && frac_dx == 0
                interpolatedBlock(r, c) = uint8(top_left * (1 - frac_dy) + bottom_left * frac_dy);

            % Diagonal interpolation (MV = [1, 1])
            else
                top_interp = top_left * (1 - frac_dx) + top_right * frac_dx;
                bottom_interp = bottom_left * (1 - frac_dx) + bottom_right * frac_dx;
                interpolatedBlock(r, c) = uint8(top_interp * (1 - frac_dy) + bottom_interp * frac_dy);
            end
        end
    end
end