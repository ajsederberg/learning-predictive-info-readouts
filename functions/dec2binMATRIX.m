function binmat = dec2binMATRIX(dec_nums, num_bits)
% dec_nums is an Nx1 list of base-10 integers. num_bits is the number of
% binary bits to write numbers with. Converts to base-2, with
% least-significant bit at the right


if nargin == 1
    num_bits = ceil(log2(max(dec_nums)));
else
    num_bits = max(num_bits, ceil(log2(max(dec_nums))));
end

binmat = zeros(length(dec_nums), num_bits+1);
for i_bit = 0:num_bits
    bit_values = floor(dec_nums/(2^(num_bits - i_bit)));
    binmat(:, i_bit+1) = bit_values;
    dec_nums = dec_nums - bit_values*2^(num_bits-i_bit);
    
    
end