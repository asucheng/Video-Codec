function decoded_values = A1_Q4_expGolombDecode(encoded_bitstream)
    % Split the bitstream by spaces to get individual encoded values
    encoded_bitstream = strtrim(encoded_bitstream);
    encoded_values = strsplit(encoded_bitstream, ' ');
    
    % Initialize the decoded values array
    decoded_values = zeros(1, numel(encoded_values));
    
    % Decode each encoded string
    for i = 1:numel(encoded_values)
        if isempty(encoded_values{i})
            continue;  % Skip if it's an empty entry (in case of trailing spaces)
        end
        
        decoded_values(i) = exp_golomb_decode_single_value(encoded_values{i});
    end
end

function value = exp_golomb_decode_single_value(encoded_str)
    % Count the leading zeros to determine the prefix length
    num_leading_zeros = find(encoded_str ~= '0', 1) - 1;
    
    % Extract the binary part (mapped_value + 1)
    binary_str = encoded_str(num_leading_zeros + 1:end);
    mapped_value_plus_one = bin2dec(binary_str);
    
    % Calculate the mapped value
    mapped_value = mapped_value_plus_one - 1;
    
    % Reverse the mapping to get the original value
    if mod(mapped_value, 2) == 0
        value = mapped_value / 2;
    else
        value = -(mapped_value + 1) / 2;
    end
end