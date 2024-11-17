function [predicted_block, mode, split_flag] = A2_Q2_intraPredictForBlock(frame, i, j, block_height, block_width, VBSenable, QP, Lambda)
    current_block = frame(i:i+block_height-1, j:j+block_width-1);

    % Initialize predictors
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

    % Compute RD cost for non-split block
    residual_horizontal = current_block - horizontal_pred;
    residual_vertical = current_block - vertical_pred;

    % RD cost for horizontal prediction
    [J_horizontal, ~] = A2_Q2_computeRD(residual_horizontal, 0, QP, Lambda);

    % RD cost for vertical prediction
    [J_vertical, ~] = A2_Q2_computeRD(residual_vertical, 1, QP, Lambda);

    % Choose prediction mode with lower RD cost
    if J_horizontal <= J_vertical
        predicted_block = horizontal_pred;
        mode = 0; % Horizontal prediction
        J_non_split = J_horizontal;
    else
        predicted_block = vertical_pred;
        mode = 1; % Vertical prediction
        J_non_split = J_vertical;
    end

    % Handle Variable Block Size (VBS)
    split_flag = false;
    if VBSenable && block_height > 4 && block_width > 4
        % Divide the current block into four sub-blocks
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
            % Predict sub-block using the same logic
            [predicted_sub_block, sub_mode] = A2_Q2_nonSplitPrediction(frame, ...
                i + (k > 2) * split_size, j + (mod(k - 1, 2) == 1) * split_size, split_size, split_size);

            residual_sub_block = sub_blocks{k} - predicted_sub_block;

            % Compute RD cost for each sub-block
            [J_sub, ~] = A2_Q2_computeRD(residual_sub_block, sub_mode, QP, Lambda);
            J_split = J_split + J_sub;
        end

        % Compare RD costs
        if J_split < J_non_split
            split_flag = true;
        end
    end
end