function [predicted_block, mode] = A1_Q4_intraPredictForBlock(frame, i, j, block_height, block_width, VBSenable)
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

    mae_horizontal = mean(abs(horizontal_pred(:) - current_block(:)), 'all');
    mae_vertical = mean(abs(vertical_pred(:) - current_block(:)), 'all');

    if mae_horizontal <= mae_vertical
        predicted_block = horizontal_pred;
        mode = 0;
    else
        predicted_block = vertical_pred;
        mode = 1;
    end
end