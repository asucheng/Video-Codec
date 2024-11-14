function bitcount = Q4_bitcountFromArray(data)
    binary_str = '';

    for i = 1:length(data)
        bin_value = dec2bin(typecast(int8(data(i)), 'uint8'), 8); 
        binary_str = [binary_str, bin_value];
    end

    bitcount = length(binary_str);
end
