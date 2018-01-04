function [is_new, old_ind, num_matches] = checkPartitionUniqueness_RPMethod(uniq_partition_set, test_partitions)
% This returns the set of unique partitions using random projections of the
% full partitions set
% Won't work if the partitions are "too high" dimension -- seems ok up to
% N = 7. 


rng('default'); % control rand: prevent (highly unlikely) pathological cases. 
random_prior = sort(rand(size(uniq_partition_set, 2), 1), 'descend');

random_projweights = uniq_partition_set*random_prior;

test_inds = test_partitions*random_prior;
  
scaling_fac = max(20000, min(1e6, round(0.1/min(abs(diff(random_projweights))))));

random_projweights = floor(scaling_fac*random_projweights);
test_inds = floor(scaling_fac*test_inds);

% find where test_inds is a member of random_projweights
% old_ind(i) has the location of test_inds(i) in random_projweights (0 if
% not in random_projweights)
[is_old, old_ind] = ismember( test_inds, random_projweights);

is_new = ~is_old; 

num_matches = zeros(size(test_partitions, 1), 1);
for ii = 1:length(num_matches)
    try num_matches(ii) = sum(test_partitions(ii, :) == uniq_partition_set(old_ind(ii), :));
    catch
        num_matches(ii) = 0;
    end
end

display([num2str(sum(num_matches == size(test_partitions, 2))) ' of ' num2str(size(test_partitions, 1)) ' matched'])

% if sum(num_matches == size(test_partitions, 2)) < size(test_partitions, 1)
%     keyboard
% end