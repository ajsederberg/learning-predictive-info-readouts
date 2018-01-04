function S = generateRandomSubsets(full_set, max_size, num_draws)
%GENERATERANDOMSUBSETS generates randomly drawn subsets of the set full_set
%with sizes from 1 to max_size. There are up to  num_draws of each set.

if length(max_size) == 2
    min_size = max_size(1);
    max_size = max_size(2);
else
    min_size = 1;
end


N = length(full_set);
if N < max_size
    S = [];
    return;
end
max_size = min(max_size, N - 1);

if N <= min_size
    S = [];
    return;
end
S = zeros(N*num_draws, max_size);

ctr = 0;
for s = min_size:max_size
    
    % first, if nchoosek(N, s) <= num_draws, just take all possible sets
    if nchoosek(N, s) <= num_draws
        ss = generateAllCombinations(full_set, s, num_draws);
    else
        % initialize matrix for num_draws sets of size s
        ss = zeros(num_draws, s);

        for k = 1:num_draws
            % randomly draw s numbers from the full_set of size N
            inds = randperm(N, s);
            % translate indices to IDs, sort from small to large
            ss(k, :) = sort(full_set(inds));
        end

        % turn each set of s elements into a base-N number
        num_ss = zeros(size(ss, 1), 1);
        for i_s = 1:s
            num_ss = num_ss + ss(:, i_s)*N^(i_s-1);
        end
        % identify unique sets
        [~, uniq_s] = unique(num_ss);
        % keep only the unique sets
        ss = ss(uniq_s, :);
    end 
    
    if numel(ss) == max_size
        ss = ss(:)';    % make sure it's a row vector if there's only one set drawn
    end
    S(ctr+1:(ctr + size(ss, 1)), 1:s) = ss;
    % keep track of how many sets have been added 
    ctr = ctr + size(ss, 1);
end

S = S(1:ctr, :);