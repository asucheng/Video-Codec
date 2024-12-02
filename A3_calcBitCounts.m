function bitCounts = A3_calcBitCounts(encoded_bitstream)
    encoded_bitstream_no_space = strrep(encoded_bitstream, ' ', '');
    bitCounts = length(encoded_bitstream_no_space);
end