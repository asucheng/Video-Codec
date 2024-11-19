function encoded_bitstream = A1_Q4_expGolombEncode(value)
    if numel(value) > 1
        encoded_bitstream = '';
        for i = 1:length(value)
            encoded_str = exp_golomb_encode_single_value(value(i));
            encoded_bitstream = [encoded_bitstream, encoded_str, ' '];  % 用空格分隔
        end
    else
        encoded_bitstream = exp_golomb_encode_single_value(value);
    end
end


function encode_str = exp_golomb_encode_single_value(value)
    if value >= 0
        mapped_value = 2 * value;
    else
        mapped_value = -2 * value - 1;
    end
    binary_str = dec2bin(mapped_value + 1);
    prefix = repmat('0', 1, length(binary_str) - 1);
    encode_str = [prefix, binary_str];
end