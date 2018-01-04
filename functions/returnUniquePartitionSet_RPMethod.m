function [uniq_part_set, uniq_inds, full_inds] = returnUniquePartitionSet_RPMethod(full_partition_set)
% This returns the set of unique partitions using random projections of the
% full partitions set
% Won't work if the partitions are "too high" dimension -- don't know where
% that line is yet. 

rng('default'); % control rand: check for pathological cases. 
random_prior = sort(rand(size(full_partition_set, 2), 1), 'descend');

random_projweights = full_partition_set*random_prior;

% multiply everything by 0.1/minimum difference: make this an integer
% problem. 
scaling_fac = min(1e14, round(1/min(abs(diff(random_projweights)))));

random_projweights = floor(scaling_fac*random_projweights);

[~, uniq_inds, full_inds] = unique(random_projweights);

uniq_part_set = full_partition_set(uniq_inds, :);