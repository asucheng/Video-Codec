function decoded_coeffs = A1_Q4_rleDecode(rle_encoded, blockSize)
    % Initialize an empty array for decoded coefficients
    decoded_coeffs = [];
    i = 1;

    % Loop through the encoded data
    while i <= length(rle_encoded)
        value = rle_encoded(i);

        if value == 0
            % End of encoding
            break;
        elseif value < 0
            % Negative value indicates the length of a non-zero sequence
            count = -value;  % Get the number of non-zero coefficients
            non_zero_values = rle_encoded(i+1:i+count);  % Extract these coefficients
            decoded_coeffs = [decoded_coeffs, non_zero_values];  % Append to output
            i = i + count + 1;  % Move index past this sequence
        else
            % Positive value indicates the count of zeros
            decoded_coeffs = [decoded_coeffs, zeros(1, value)];  % Append zeros
            i = i + 1;  % Move to the next entry
        end
    end

    % Ensure the length matches blockSize^2 (e.g., 64 for 8x8 blocks)
    expectedLength = blockSize^2;
    if length(decoded_coeffs) < expectedLength
        %Pad with zeros if necessary
        decoded_coeffs = [decoded_coeffs, zeros(1, expectedLength - length(decoded_coeffs))];
    end

end