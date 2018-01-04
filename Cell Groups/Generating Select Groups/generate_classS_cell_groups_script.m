% forgot to save script on home computer, writing code on imac
% this is for generating cell groups

% multimovie first

mmc = load('cell_list_20111103', 'cells_to_use');
cells_to_use = mmc.cells_to_use;

X = load('2011-11-03_multmovie_binned.mat', 'binned');
mmdata = X.binned;

[numReps, numFrames, numCells, numMovies] = size(mmdata);

aveFRbyMovieMM = squeeze(60*mean(mean(mmdata, 2), 1));
% discard from 'cells_to_use' those with mean firing rate of less than
% 0.25Hz
cells_to_use = cells_to_use(aveFRbyMovieMM(cells_to_use, 4) > 0.25);

% and from the first fishmovie
fmc = load('2011-02-25_barfishcheck_experiment_and_calcs_PNAS_paper', 'fishdata');
fmdata = fmc.fishdata;
aveFRbyMovieFM = squeeze(60*mean(mean(fmdata, 3), 2));


%%

load('multimovie_pairIpred')

dinputIpred = triu(inputIpred);
inputIpred = tril(inputIpred);

dinputIpred = dinputIpred + dinputIpred';
inputIpred = inputIpred + inputIpred';

[~, ord] = sort(sum(inputIpred), 'descend');

figure()
imagesc(inputIpred(ord,ord) - dinputIpred(ord,ord), [0 0.03])

final_cells_to_use = cells_to_use(ord(1:19));
%%
figure(),

for i_sp = 1:length(final_cells_to_use)
    subplot(4, 5, i_sp)
    imagesc(squeeze(mmdata(:, :, final_cells_to_use(i_sp), 4)), [0 1])
end

%% select sets with the most predictive information
% min_pI_thr = 0.0747;
% best_MMdata_sets = cell_sets(inputIpred4 > min_pI_thr, :);
% best_MMDS_inputPI = inputIpred4(inputIpred4 > min_pI_thr);
% best_MMDS_dinputPI = dinputIpred4(inputIpred4 > min_pI_thr);
% saveInParfor('multimovie_setsOfFour', best_MMdata_sets, best_MMDS_dinputPI, ...
%     best_MMDS_inputPI, min_pI_thr)
clear;
load('multimovie_quadIpred')

num_groups = 1000;


cell_IDs = unique(cell_sets);
for cells_in_group = 3:10
    
    keep_trying = true;
    
    if cells_in_group ~= 4
        sets_of_N = zeros(num_groups, cells_in_group);
    
        while keep_trying
            sets_of_N = generateRandomSubsets(cell_IDs, ...
                cells_in_group*[1 1], num_groups+100);
            if size(sets_of_N, 1) > num_groups
                sets_of_N = sets_of_N(1:num_groups, :);
            end
            subplot(1,1,1)
            histogram(sets_of_N, 'binmethod', 'integers')
            title(num2str(size(sets_of_N, 1)))
            keep_trying = strcmp(input('keep trying?y/n'), 'y');

        
        end
        inputIpred_gp = [];
        
    else
        while keep_trying
            rand_set = randperm(size(cell_sets, 1), num_groups);
            sets_of_N = cell_sets(rand_set, :);
            subplot(121)
            histogram(inputIpred4(rand_set))
            subplot(122)
            histogram(sets_of_N, 'binmethod', 'integers')
            inputIpred_gp = inputIpred4(rand_set);

            keep_trying = strcmp(input('keep trying?y/n'), 'y');
    
        end
        
    end
    
    save(['multimovie_setsOf' num2str(cells_in_group)], 'sets_of_N', 'inputIpred_gp')
    
end

%% choose pairs from 'fishmovie' similarly

load('fishmovie_pairIpred')

dinputIpred_FD = triu(inputIpred_FD);
inputIpred_FD = tril(inputIpred_FD);

dinputIpred_FD = dinputIpred_FD + dinputIpred_FD';
inputIpred_FD = inputIpred_FD + inputIpred_FD';

[~, ord_FD] = sort(sum(inputIpred_FD), 'descend');

figure()
imagesc(inputIpred_FD(ord_FD, ord_FD) - 2*dinputIpred_FD(ord_FD, ord_FD), [0 0.06])
%% sample all groups of four

final_cells_to_use_FD = ord_FD(1:25);
special_cells = final_cells_to_use_FD(1:3);

num_sets = nchoosek(length(final_cells_to_use_FD), 4);

inputIpred_FD4 = zeros(num_sets, 1);
dinputIpred_FD4 = zeros(num_sets, 1);
cell_sets_FD = zeros(num_sets, 4);
%%
ctr = 0;
tic;
for ii = 2:length(final_cells_to_use_FD)
    for jj = 1:(ii-1)
        for kk = 1:(jj-1)
            for mm = 1:(kk-1)
                
                ctr = ctr + 1;
                cell_sets_FD(ctr, :) = final_cells_to_use_FD([ii jj kk mm]);
                
                input_words = squeeze(fmdata(cell_sets_FD(ctr, :), :, :)) ~= 0;
                input_words = permute(input_words, [3 2 1]);
                input_words = reshape(input_words, size(input_words, 1)*size(input_words, 2), size(input_words, 3));
                [~, inpI] = computePredictiveInfo_FSS(input_words, input_words(:, 1), 1, nan);

                inputIpred_FD4(ctr) = inpI(1);
                dinputIpred_FD4(ctr) = inpI(2);
                
                if mod(ctr, 1000) == 0
                    display([num2str(ctr) ' done of ' num2str(num_sets) ' in ' num2str(toc, '%1.0f') ' seconds'])
                end
            end
        end
    end
end
%%
load('fishmovie_quadIpred', 'inputIpred_FD4', 'dinputIpred_FD4', 'cell_sets_FD')

%%

% put 4-cell groups into pair-pair matrix
n_fctu = length(final_cells_to_use_FD);
predI_prpr = zeros(n_fctu^2);

all_pairs = zeros(n_fctu^2, 2);
all_pairs(:, 1) = reshape(repmat(1:n_fctu, n_fctu, 1), n_fctu^2, 1);
all_pairs(:, 2) = reshape(permute(repmat(1:n_fctu, n_fctu, 1), [2 1]), n_fctu^2, 1);

sorted_cell_sets_FD = sort(cell_sets_FD, 2);
all_pair_inds = sort(final_cells_to_use_FD(all_pairs), 2);

for i_pair = 2:size(all_pair_inds, 1)
    for j_pair = 1:i_pair
        this_set = sort([all_pair_inds(i_pair, :) all_pair_inds(j_pair, :)]);
        if length(unique(this_set)) == 4
            [~, which_set] = intersect(sorted_cell_sets_FD, this_set, 'rows');
            
            predI_prpr(i_pair, j_pair) = inputIpred_FD4(which_set);
            
        end
    end
end

predI_prpr = predI_prpr + predI_prpr';
%%
load('fishmovie_quadIpred', 'inputIpred_FD4', 'dinputIpred_FD4', ...
    'cell_sets_FD', 'predI_prpr', 'all_pairs', 'all_pair_inds')

%%

makeMyFigure(20, 10);
subplot(221)
histogram(inputIpred_FD4)

subplot(222)
[~, pr_ord] = sort(sum(predI_prpr));
imagesc(predI_prpr(pr_ord, pr_ord), [0 0.15])

subplot(2, 2, 3)
histogram(all_pairs(pr_ord(26:487), :), 'binmethod', 'integers')

subplot(2, 2,4)
histogram(all_pairs(pr_ord(488:625), :), 'binmethod', 'integers')

% observation: all the best groups have some subset of the same three best
% cells in them. Best not to sample based on pred-I: instead select from
% the cells that initially made the cut. 
%% fishmovie groups

clear
load('fishmovie_quadIpred')
num_groups = 1000;

special_cells = [13 24 27];

cell_IDs = unique(cell_sets_FD);
for cells_in_group = 3:10
    
    keep_trying = true;
    
    if cells_in_group ~= 4
        sets_of_N = zeros(num_groups, cells_in_group);
    
        while keep_trying
            sets_of_N = generateRandomSubsets(cell_IDs, ...
                cells_in_group*[1 1], num_groups+100);
            if size(sets_of_N, 1) > num_groups
                sets_of_N = sets_of_N(1:num_groups, :);
            end
            subplot(1,1,1)
            histogram(sets_of_N, 'binmethod', 'integers')
            title(num2str(size(sets_of_N, 1)))
            keep_trying = strcmp(input('keep trying?y/n'), 'y');

        
        end
        inputIpred_gp = [];
        
    else
        while keep_trying
            rand_set = randperm(size(cell_sets_FD, 1), num_groups);
            sets_of_N = cell_sets_FD(rand_set, :);
            subplot(121)
            histogram(inputIpred_FD4(rand_set))
            subplot(122)
            histogram(sets_of_N, 'binmethod', 'integers')
            inputIpred_gp = inputIpred_FD4(rand_set);

            keep_trying = strcmp(input('keep trying?y/n'), 'y');
    
        end
        
    end
    
    save(['fishmovie_setsOf' num2str(cells_in_group)], 'sets_of_N', 'inputIpred_gp', 'special_cells')
    
end

