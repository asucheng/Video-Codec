function rle_encoded = A1_Q4_rleEncode(coeffs)
    % Run-Length Encoding (RLE) for quantized coefficients
    rle_encoded = [];
    n = length(coeffs);
    i = 1;

    while i <= n
        if coeffs(i) ~= 0
            j = i;
            while j <= n && coeffs(j) ~= 0
                j = j + 1;
            end
            rle_encoded = [rle_encoded, -(j - i), coeffs(i:j-1)];
            i = j;
        else
            j = i;
            while j <= n && coeffs(j) == 0
                j = j + 1;
            end

            if j > n
                break;
            end
            
            rle_encoded = [rle_encoded, j - i];
            i = j;
        end
    end

    rle_encoded = [rle_encoded, 0];
end