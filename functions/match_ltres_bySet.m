function matched_set_inds = match_ltres_bySet(all_lt_res)
% WRITE THIS FUNCTION TO 


all_set_strings = cellfun(@(x) cellfun(@(y) y.i_set_string, x, 'uniformoutput', false),...
    all_lt_res, 'uniformoutput', false);


uniq_set_strings = all_set_strings{1};
for i_set2 = 2:length(all_set_strings)
    uniq_set_strings = union(uniq_set_strings, all_set_strings{i_set2});
end

matched_set_inds = zeros(length(uniq_set_strings), length(all_lt_res));


for i_set = 1:length(all_set_strings)
    
    [~, uniq_set_inds, all_set_inds] = intersect(uniq_set_strings, all_set_strings{i_set});
    
    matched_set_inds(uniq_set_inds, i_set) = all_set_inds; 


end

keep_inds = ~any(matched_set_inds == 0, 2);

matched_set_inds = matched_set_inds(keep_inds, :);