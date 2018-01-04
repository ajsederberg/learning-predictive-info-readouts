function [set_detWt, wt_dom] = findDeterminativeInput(ic_weights, learned_predI)

if size(ic_weights, 2) ~= size(learned_predI, 2)
    ic_weights = ic_weights';
    if size(ic_weights, 2) ~= size(learned_predI, 2)
        ic_weights = ic_weights';
        learned_predI = learned_predI';
    end
end

num_sets = size(learned_predI, 1);

xc_mat = corrcoef([learned_predI ; ic_weights]');

input_corr = xc_mat(1:num_sets, num_sets+1:end);

[~, set_detWt] = max(input_corr, [], 2);

sorted_corr = sort(input_corr, 2);
wt_dom = (sorted_corr(:, end) - sorted_corr(:, end-1))./sorted_corr(:, end);