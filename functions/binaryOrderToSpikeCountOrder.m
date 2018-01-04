function [sc_ord, sc, bin_mat] = binaryOrderToSpikeCountOrder(n)
% function to get the reordering of binary words so that words are sorted
% by spike count (i.e. 000 001 010 100 011 110 101 111) rather than binary
% number order (000 001 010 011 100 101 110 111)

bin_rep = dec2bin(0:(2^n -1 ));

bin_mat = zeros(size(bin_rep));
for ii = 1:n
    bin_mat(:, ii) = str2num(bin_rep(:, ii));
end

sc = sum(bin_mat, 2);
[~, sc_ord] = sort(sc);