function [predicted_block, mode, split_flag] = A2_Q2_intraPredictForBlock(frame, i, j, block_height, block_width, VBSenable, QP, Lambda)
    current_block = frame(i:i+block_height-1, j:j+block_width-1);

    % Generate vertical and horizontal predictors
    if j == 1
        vertical_pred = 128 * ones(block_height, block_width);
    else
        vertical_pred = repmat(frame(i:i+block_height-1, j-1), 1, block_width);
    end

    if i == 1
        horizontal_pred = 128 * ones(block_height, block_width);
    else
        horizontal_pred = repmat(frame(i-1, j:j+block_width-1), block_height, 1);
    end

    % Non-split block logic
    if ~VBSenable
        % Calculate MAE for horizontal and vertical predictions
        mae_horizontal = mean(abs(horizontal_pred(:) - current_block(:)), 'all');
        mae_vertical = mean(abs(vertical_pred(:) - current_block(:)), 'all');
    
        % Choose the mode with lower MAE
        if mae_horizontal <= mae_vertical
            predicted_block = horizontal_pred;
            mode = 0; % Horizontal prediction
        else
            predicted_block = vertical_pred;
            mode = 1; % Vertical prediction
        end
        split_flag = false; % No split for non-VBS
        return
    end

    % RD-based prediction for VBS
    % Compute residuals for horizontal and vertical predictions
    residual_horizontal = current_block - horizontal_pred;
    residual_vertical = current_block - vertical_pred;

    % Compute RD cost for horizontal and vertical modes
    [J_horizontal, ~] = A2_Q2_computeRD(residual_horizontal, 0, QP, Lambda);
    [J_vertical, ~] = A2_Q2_computeRD(residual_vertical, 1, QP, Lambda);

    % Choose the best prediction mode for non-split
    if J_horizontal <= J_vertical
        predicted_block = horizontal_pred;
        mode = 0; % Horizontal prediction
        J_nonsplit = J_horizontal;
    else
        predicted_block = vertical_pred;
        mode = 1; % Vertical prediction
        J_nonsplit = J_vertical;
    end

    % Split block logic if block size permits
    split_flag = false; % Default to no split
    if block_height > 2 && block_width > 2
        % Split the current block into four sub-blocks
        split_size = floor(block_height / 2);
        sub_blocks = {
            current_block(1:split_size, 1:split_size),         % Top-left
            current_block(1:split_size, split_size+1:end),    % Top-right
            current_block(split_size+1:end, 1:split_size),    % Bottom-left
            current_block(split_size+1:end, split_size+1:end) % Bottom-right
        };

        % Initialize RD cost for split
        J_split = 0;

        % Process each sub-block
        for k = 1:4
            % Compute sub-block coordinates
            sub_row = i + (k > 2) * split_size; % Top or bottom
            sub_col = j + (mod(k - 1, 2) == 1) * split_size; % Left or right
            sub_block = sub_blocks{k};

            % Generate predictors for sub-block
            if sub_col == j
                vertical_pred_sub = 128 * ones(split_size, split_size);
            else
                vertical_pred_sub = repmat(frame(sub_row:sub_row+split_size-1, sub_col-1), 1, split_size);
            end

            if sub_row == i
                horizontal_pred_sub = 128 * ones(split_size, split_size);
            else
                horizontal_pred_sub = repmat(frame(sub_row-1, sub_col:sub_col+split_size-1), split_size, 1);
            end

            % Compute RD cost for sub-block
            residual_horizontal_sub = sub_block - horizontal_pred_sub;
            residual_vertical_sub = sub_block - vertical_pred_sub;

            [J_horizontal_sub, ~] = A2_Q2_computeRD(residual_horizontal_sub, 0, QP, Lambda);
            [J_vertical_sub, ~] = A2_Q2_computeRD(residual_vertical_sub, 1, QP, Lambda);

            % Add the lower RD cost to the total split cost
            J_split = J_split + min(J_horizontal_sub, J_vertical_sub);
        end

        % Add overhead bit cost for split flag
        J_split = J_split + Lambda; % One bit overhead for split flag

        % Compare RD costs for non-split and split
        if J_split < J_nonsplit
            split_flag = true; % Decide to split
            predicted_block = NaN; % To indicate that block is split (handled at a higher level)
        end
    end
end