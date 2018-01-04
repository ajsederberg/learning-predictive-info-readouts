function [uniq_part_set, uniq_inds, full_inds] = returnUniquePartitionSet(full_partition_set)
% returns the set of unique partitions. Each row of full_partition_set
% corresponds to a partition. 

uniq_part_set = 0*full_partition_set;

num_sets = size(full_partition_set, 1);
% uniq_inds is the index of unique elements of full_partition_set
uniq_inds = zeros(num_sets, 1);
% full_inds is the index of each element of full_partition_set in the
% uniq_part_set
full_inds = zeros(num_sets, 1);

% enter first partition into set of unique partitions
u_ctr = 1;
uniq_part_set(1, :) = full_partition_set(1, :);
uniq_inds(u_ctr) = 1;
full_inds(1) = 1;
% check each subsequent partition for uniqueness against the existing set
for ii = 2:num_sets
    
    [is_new, old_ind] = ...
        checkPartitionUniqueness(uniq_part_set(1:u_ctr, :), full_partition_set(ii, :));
    if is_new
        u_ctr = u_ctr + 1;
        uniq_part_set(u_ctr, :) = full_partition_set(ii, :);
        uniq_inds(u_ctr) = ii;
        full_inds(ii) = u_ctr;
    else
        full_inds(ii) = old_ind;
    end
end

uniq_part_set = uniq_part_set(1:u_ctr, :);
uniq_inds = uniq_inds(1:u_ctr);