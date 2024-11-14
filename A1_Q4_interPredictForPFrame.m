function [encoded_frame, reconstructed_frame, bitstream] = Q4_interPredictForPFrame(frame, reference_frame, block_size, QP)
    Q_Matrix = Q4_generateQMatrix(block_size, QP);

    [h, w] = size(frame);
    reconstructed_frame = zeros(h, w);
    encoded_frame = zeros(h, w);
    bitstream = '';

    Q_Matrix = Q4_generateQMatrix(block_size, QP);
    previous_mv = [0, 0];
    search_range = 2;

    for i = 1:block_size:h
        for j = 1:block_size:w
            block_height = min(block_size, h - i + 1);
            block_width = min(block_size, w - j + 1);
            current_block = frame(i:i+block_height-1, j:j+block_width-1);

            [best_mv, predicted_block] = getMotionVector(reference_frame, current_block, i, j, search_range);

            % update the motion vector differentiation
            diff_mv = best_mv - previous_mv;
            previous_mv = best_mv;

            % add the bitstream of motion vector diff
            bitstream = [bitstream, Q4_expGolombEncode(diff_mv(1)), Q4_expGolombEncode(diff_mv(2))];

            % get the residuals of block
            residual_block = current_block - predicted_block;

            % encode the residual block
            encoded_residual_block = Q4_quantizeBlockAfterDCT(residual_block, Q_Matrix);

            % decode the encoded block and add it to prediction
            decoded_residual_block = Q4_idctAfterDequantizeBlock(encoded_residual_block, Q_Matrix);
            reconstructed_block = predicted_block + decoded_residual_block;
            reconstructed_block = max(min(reconstructed_block, 255), 0);
            reconstructed_frame(i:i+block_height-1, j:j+block_width-1) = reconstructed_block;
        end
    end
end


function [best_mv, predicted_block] = getMotionVector(reference_frame, current_block, i, j, search_range)
    [height, width] = size(reference_frame);
    block_height = size(current_block, 1);
    block_width = size(current_block, 2);

    min_mae = Inf;
    candidates = [];

    for dx = -search_range:search_range
        for dy = -search_range:search_range
            ref_row = i + dx;
            ref_col = j + dy;

            if ref_row > 0 && ref_col > 0 && ref_row + block_height - 1 <= height && ref_col + block_width - 1 <= width
                candidate_block = reference_frame(ref_row:ref_row+block_height-1, ref_col:ref_col+block_width-1);

                mae = mean(abs(double(current_block(:)) - double(candidate_block(:))));
                
                if mae < min_mae
                    min_mae = mae;
                    candidates = [dy, dx];
                elseif mae == min_mae
                    % Add this candidate if it has the same MAE
                    candidates = [candidates; [dy, dx]];
                end
            end
        end
    end

    % Apply the triple logic to resolve ties
    best_mv = resolveTies(candidates);

    % Calculate the coordinates of the predicted block
    best_row = i + best_mv(2);
    best_col = j + best_mv(1);

    % Extract the predicted block from the reference frame
    predicted_block = reference_frame(best_row:best_row+block_height-1, best_col:best_col+block_width-1);
end


% Helper function to resolve ties using L1 norm and coordinate checks
function bestMatch = resolveTies(candidates)
    % Step 1: Find candidates with the smallest L1 norm
    smallestL1 = min(sum(abs(candidates), 2));
    candidatesL1 = candidates(sum(abs(candidates), 2) == smallestL1, :);

    % Step 2: Find candidates with the smallest y-coordinate
    smallestY = min(candidatesL1(:, 2));
    candidatesY = candidatesL1(candidatesL1(:, 2) == smallestY, :);

    % Step 3: Find the candidate with the smallest x-coordinate
    smallestX = min(candidatesY(:, 1));
    bestMatch = candidatesY(candidatesY(:, 1) == smallestX, :);
end