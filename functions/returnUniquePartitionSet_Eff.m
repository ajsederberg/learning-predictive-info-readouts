function [uniq_part_set, uniq_inds, full_inds] = returnUniquePartitionSet_Eff(full_partition_set)
% returns the set of unique partitions. Each row of full_partition_set
% corresponds to a partition. 



num_sets = size(full_partition_set, 1);



%%
full_p11 = 2*(full_partition_set - 0.5);
block_length = 100;
num_blocks = ceil(num_sets/block_length);

% keep_inds will have "true" wherever an element is unique within its
% block. 
keep_inds = false(num_sets, 1);
for i_block = 1:num_blocks
    block_inds = (1:block_length) + block_length*(i_block - 1);
    block_inds = block_inds(block_inds <= size(full_partition_set, 1));
    this_block_length = length(block_inds);
    
    part_p11 = full_p11(block_inds, :);
    part_p11_sq = part_p11*(part_p11');
    
    block_uniq_inds = false(this_block_length, 1);
    
    [i_eq, j_eq] = find(part_p11_sq == size(full_partition_set, 2));
    
    while ~isempty(i_eq)
        
        block_uniq_inds(i_eq(1)) = true;        
        equal_to_ieq = i_eq == i_eq(1);
        
        j_remove = j_eq(equal_to_ieq);
%         block_uniq_inds(j_remove) = false;
        
        i_was_removedj = ismember(i_eq, j_remove);
%         block_uniq_inds(i_was_removedj) = false;

        
        % remove all pairs eithe requal to the first component (i_eq(1)) or
        % equal to any of the pairs that were equal to the first component
        i_eq = i_eq(~i_was_removedj & ~equal_to_ieq);
        j_eq = j_eq(~i_was_removedj & ~equal_to_ieq);

        
    end
    
    keep_inds(block_inds(block_uniq_inds)) = true;
end

partial_partition_set = full_partition_set(keep_inds, :);

[uniq_part_set, partial_uniq_inds, partial_full_inds] = returnUniquePartitionSet(partial_partition_set);

full_kept_inds = find(keep_inds);
uniq_inds = full_kept_inds(partial_uniq_inds);


% return nan's -- this is not function yet
full_inds = nan(size(full_partition_set, 1), 1);
