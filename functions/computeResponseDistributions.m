function [q_i, obs_ct] = computeResponseDistributions(input_wordnums, output_spikes)

% first, condition inputs
[m,n] = size(input_wordnums);
% second dimension should be time (long dimension)
if m > n
    input_wordnums = input_wordnums';
end
% if not already in number format (i.e. in "word" format) convert to word
% numbers
if n ~= 1 && m ~= 1
    numL0 = min(m,n);
    input_wordnums = (2.^(0:numL0-1))*input_wordnums;
    num_words = 2^numL0;
else
    num_words = 2^ceil(log2(max(input_wordnums)));
end

% second dimension should be time (long dimension)
[m2, n2] = size(output_spikes);
if m2 > n2
    output_spikes = output_spikes';
end
    
num_out = size(output_spikes, 1);

q_i = zeros(num_out, num_words);

obs_ct = zeros(1, num_words);
% compute response distribution
for i_out = 1:num_out
    for i_word = 1:num_words
        q_i(i_out, i_word) = mean(output_spikes(i_out, input_wordnums == (i_word-1)));
        obs_ct(i_word) = sum(input_wordnums == (i_word-1));
    end
    
end

q_i(isnan(q_i)) = 0;