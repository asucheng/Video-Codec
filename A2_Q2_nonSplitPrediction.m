function [predicted_block, mode] = A2_Q2_nonSplitPrediction(frame, i, j, block_height, block_width)
    % Extract the current block
    current_block = frame(i:i+block_height-1, j:j+block_width-1);

    % Vertical prediction
    if j == 1
        vertical_pred = 128 * ones(block_height, block_width);
    else
        vertical_pred = repmat(frame(i:i+block_height-1, j-1), 1, block_width);
    end

    % Horizontal prediction
    if i == 1
        horizontal_pred = 128 * ones(block_height, block_width);
    else
        horizontal_pred = repmat(frame(i-1, j:j+block_width-1), block_height, 1);
    end

    % Compare predictions using Mean Absolute Error (MAE)
    mae_horizontal = mean(abs(horizontal_pred(:) - current_block(:)), 'all');
    mae_vertical = mean(abs(vertical_pred(:) - current_block(:)), 'all');

    % Select the mode with the lower error
    if mae_horizontal <= mae_vertical
        predicted_block = horizontal_pred;
        mode = 0; % Horizontal mode
    else
        predicted_block = vertical_pred;
        mode = 1; % Vertical mode
    end
end