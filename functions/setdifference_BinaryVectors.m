function [is_new, old_in_new_inds, old_in_new_inds_dup0] = setdifference_BinaryVectors(old_parts, new_parts)
% Compares the binary row vectors in new_parts to the binary row vectors in
% old_parts. new_parts(is_new, :) are rows that do not exist in old_parts.
% new_part(nonzeros(old_in_new_inds), :) = UNIQUE[ old_parts ] 
% where "UNIQUE[]" is row-wise uniqueness in old_parts. 


[is_new, ~] = checkPartitionUniqueness_RPMethod(...
    old_parts, new_parts);
% this find the indices such that
% old_parts = new_parts(old_in_new_inds, :)
[~, old_in_new_inds] = checkPartitionUniqueness_RPMethod(...
    new_parts, old_parts);  

% check old_in_new_inds for duplicates and set those to
% zeros
[uniq_oini, ind_oini] = unique(old_in_new_inds);
old_in_new_inds_dup0 = 0*old_in_new_inds;
old_in_new_inds_dup0(ind_oini) = uniq_oini;