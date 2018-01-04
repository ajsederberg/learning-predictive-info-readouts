function [perceptron_out, input_mat, proj_rule] = perceptronToProjectionRule(w, b)
% maps perceptron weights to a projection rule. Inputs are 0,1 and outputs
% are 0, 1. NOTE: CHANGED TO INCLUDE ZERO WORD ON 4/29/15. CHECK
% COMPATIBILITY

[m, n] = size(w);
if m == 1 || n == 1
    % make w a row vector
    w = reshape(w, 1, numel(w));
    num_inputs = length(w);
else
    num_inputs = n;
end
% be careful with this use of str2num. OK in 2014b. 
input_vec = dec2bin(0:(2^num_inputs - 1), num_inputs)';
input_mat = zeros(size(input_vec));
for i_col = 1:size(input_vec, 2)
    input_mat(:, i_col) = flipud(str2num(input_vec(:, i_col)));
end

perceptron_out = w*input_mat - b > 0;

proj_rule(perceptron_out) = '1';
proj_rule(~perceptron_out) = '0';