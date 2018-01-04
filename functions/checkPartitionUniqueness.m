function [is_new, old_ind] = checkPartitionUniqueness(partition_set, new_partition)
% Checks whether the new_partition is already in the partition set
% returns true/false (in is_new) as well as the index in the set (if the
% partition was found in the set, i.e. is_new = false)

if size(partition_set, 2) ~= length(new_partition)
    partition_set = partition_set';
    warning('partitions should be row vector. check.')
end

is_new = true;
ctr = 0;
while is_new && ctr < size(partition_set, 1)
    ctr = ctr + 1;
    if all(partition_set(ctr, :) == new_partition)
        is_new = false;
        old_ind = ctr;
    end
end

if is_new
    old_ind = [];
end