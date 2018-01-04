function S = generateAllCombinations(labels, k, max_draws)
%GENERATERANDOMSUBSETS generates randomly drawn subsets of the set full_set
%with sizes from 1 to max_size. There are num_draws of each set.

N = length(labels);

% don't let it run forever
if nchoosek(N, k) > max_draws
    display('N choose k is too large. Choose a smaller k or N.')
    S = [];
    return;
end


if k == 1
    S = (1:N)';
    labels = reshape(labels, length(labels), 1);
elseif k == N
    S = 1:N;
else
    
    nextS = generateAllCombinations(2:N, k - 1, max_draws);
    partS = [ones(size(nextS, 1), 1) nextS];
    restS = generateAllCombinations(2:N, k, max_draws);
    S = [partS; restS];
end

S = labels(S);