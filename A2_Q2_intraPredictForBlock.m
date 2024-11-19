function [predicted_block, mode, split_flag] = A2_Q2_intraPredictForBlock(frame, i, j, block_height, block_width, VBSenable, QP, Lambda)
    current_block = frame(i:i+block_height-1, j:j+block_width-1);

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

    % Non-split block logic (matches old version)
    if ~VBSenable
        % Calculate MAE for horizontal and vertical predictions
        mae_horizontal = mean(abs(horizontal_pred(:) - current_block(:)), 'all');
        mae_vertical = mean(abs(vertical_pred(:) - current_block(:)), 'all');
    
        if mae_horizontal <= mae_vertical
            predicted_block = horizontal_pred;
            mode = 0;
        else
            predicted_block = vertical_pred;
            mode = 1;
        end
        split_flag= false;
        return
    end

    % RD-based prediction logic for VBS (if enabled)
    residual_horizontal = current_block - horizontal_pred;
    residual_vertical = current_block - vertical_pred;

    [J_horizontal, ~] = A2_Q2_computeRD(residual_horizontal, 0, QP, Lambda);
    [J_vertical, ~] = A2_Q2_computeRD(residual_vertical, 1, QP, Lambda);

    if J_horizontal <= J_vertical
        predicted_block = horizontal_pred;
        mode = 0;
        J_non_split = J_horizontal;
    else
        predicted_block = vertical_pred;
        mode = 1;
        J_non_split = J_vertical;
    end

    % Handle splitting if block size permits (currently optional)
    split_flag = false;
    if block_height > 4 && block_width > 4
        % Split logic omitted for simplicity in alignment with the old version
    end
end