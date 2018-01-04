% Script to generate all necessary numbers for the manuscript. 
% make Figure 4: the aggregated learning results for sets of 4, 7, 10
set_sizes = [4 7 10];
% makeMyFigure(30, 20);
lt_symbols = {'o', '+', '^', 'd'};
learning_rule_labels = {'pair'};
set_colors = [0 0 0; 0.5 0.5 0.5; lines(1)]; 
i_lt = 1;

tau_future_vals = -9:2:3;

example_inds = [63 1 1];

frac_edges = linspace(0, 1.5, 51);
frac_bins = edges2bins(frac_edges);


learned_stimPredI = cell(length(set_sizes), 1);
learned_stimFR = cell(length(set_sizes), 1);
optimal_stimPredI = cell(length(set_sizes), 1);
optimal_stimFR = cell(length(set_sizes), 1);
word_stimPredI = cell(length(set_sizes), 1);
word_stimFR = cell(length(set_sizes), 1);
ic_stability = cell(length(set_sizes), 1);
set_processedRes = cell(length(set_sizes), 1);
set_PredI  = cell(length(set_sizes), 1);
for i_ss = 1:length(set_sizes)

set_size = set_sizes(i_ss);

[all_lt_resS, all_lt_res_infoS, processed_resultsS] = getSimulationInfoResults(set_size, ['fS' num2str(set_size)], {'vanilla'});
[all_lt_resA, all_lt_res_infoA, processed_resultsA] = getSimulationInfoResults(set_size, ['fA' num2str(set_size)], {'vanilla'});
%%
all_lt_res = cell(length(all_lt_resA),1);
all_lt_res_info = cell(length(all_lt_res),1);
pr_flds = {'prct_learned', 'learned_perc_ID', 'learned_pI', 'prct_maxFR', 'rule_sims', 'lt_fr_bins', ...
    'all_lt_res_hull_bin', 'all_lt_wb_hull', 'all_perc_sims', ...
    'all_lt_wb_hull_perc_ID', 'all_lt_perc_pI', 'all_lt_perc_fr', 'all_lt_last10dW'};

pr_flds_short = {'pI_high', 'wb_hull_all'};
clear processed_results

processed_results = processed_resultsA;

for i_lt = 1:length(all_lt_res)
    all_lt_res{i_lt} = [all_lt_resS{i_lt};all_lt_resA{i_lt}];
    all_lt_res_info{i_lt} = [all_lt_res_infoS{i_lt};all_lt_res_infoA{i_lt}];
    
    for i_pf = 1:length(pr_flds)
        processed_results.(pr_flds{i_pf}){i_lt, 1} = [processed_resultsS.(pr_flds{i_pf}){i_lt}; ...
            processed_resultsA.(pr_flds{i_pf}){i_lt}];
    end
    

end

for i_pfs = 1:length(pr_flds_short)
    processed_results.(pr_flds_short{i_pfs}) = [processed_resultsS.(pr_flds_short{i_pfs}); ...
        processed_resultsA.(pr_flds_short{i_pfs})];
end
%% load processed results
pr_flds = fieldnames(processed_results);
for i_prf = 1:length(pr_flds)
    eval([pr_flds{i_prf} ' = processed_results.' pr_flds{i_prf} ';']);
end

set_processedRes{i_ss} = processed_results;
%% load the perceptron dictionary

load(['aggregated results/fA' num2str(set_size) '_info_calcs.mat'], 'perc_lookup_struct')
%% Compute the stimulus information of the learned perceptrons
% take the selected learning type
lt_wb_hull_perc_ID = all_lt_wb_hull_perc_ID{i_lt};  % dictionary IDs of perceptron hull
lt_res_hull_bin = all_lt_res_hull_bin{i_lt};        % hull IDs of learning results
lt_res_info = all_lt_res_info{i_lt};                % info calculations (for dictionary perceptrons)
lt_learned_perc_ID = learned_perc_ID{i_lt};         % dictionary IDs of learned perceptrons

% Function to return selected row elements of a cell array of structures
info_val_fun = @(x_struct, y_inds, fld1, fld2) ...
    cellfun(@(z1, z2) z1.(fld1).(fld2)(z2, :), x_struct, y_inds, ...
    'uniformoutput', false);
% find the index of the optimal perceptrons in the perceptron dictionary
optimal_perc_ID = cellfun(@(x, y) x(y), ...
    lt_wb_hull_perc_ID, lt_res_hull_bin, 'uniformoutput', false);
% find the stimulus information for all learned perceptrons
[~, learned_stimPredI{i_ss}] = cellfun(@(x, y) findPerceptronsInStimInfo(x, y, perc_lookup_struct, 'percstim_predI_vt', tau_future_vals), ...
    lt_res_info, lt_learned_perc_ID, 'uniformoutput', false);
[~, learned_stimFR{i_ss}] = cellfun(@(x, y) findPerceptronsInStimInfo(x, y, perc_lookup_struct, 'nrb_fr_perc', []), ...
    lt_res_info, lt_learned_perc_ID, 'uniformoutput', false);

%% find the stimulus information for all optimal perceptrons
% lt_optperc_stim_info = info_val_fun(lt_res_info, optimal_perc_ID, 'X_stim', 'percstim_predI_vt');
[~, optimal_stimPredI{i_ss}] = cellfun(@(x, y) findPerceptronsInStimInfo(x, y, perc_lookup_struct, 'percstim_predI_vt', tau_future_vals), ...
    lt_res_info, optimal_perc_ID, 'uniformoutput', false);
[~, optimal_stimFR{i_ss}] = cellfun(@(x, y) findPerceptronsInStimInfo(x, y, perc_lookup_struct, 'nrb_fr_perc', []), ...
    lt_res_info, optimal_perc_ID, 'uniformoutput', false);
% tau_vec = cellfun(@(x) x.X_stim.tau_futures, lt_res_info, 'uniformoutput', false);
% find the word-word stimulus information for each group
% full_wordstim_info = cellfun( @(x) x.X_stim.wordstim_predI_vt, ...
%     lt_res_info(keep_groups), 'uniformoutput', false);
[~, word_stimPredI{i_ss}] = cellfun(@(x) findPerceptronsInStimInfo(x, [], [], 'wordstim_predI_vt', tau_future_vals), ...
    lt_res_info, 'uniformoutput', false);
word_stimFR{i_ss} = cellfun(@(x) sum(x.X_stim.nrb_fr_cells), lt_res_info, ...
    'uniformoutput', false);

set_PredI{i_ss} = cellfun(@(x) x.X_ww.word_pI(x.X_ww.tau_futures == 1), lt_res_info);

% get the std dev of weights over the last 10 simulations sets
ic_stability{i_ss} = cellfun(@(x) x', set_processedRes{i_ss}.all_lt_last10dW{1}, ...
    'uniformoutput', false);
end
%%
infoX4A = load('aggregated results/fA4_info_calcs.mat');
infoX4S = load('aggregated results/fS4_info_calcs.mat');

infoX4 = [infoX4S.info_ds_arr; infoX4A.info_ds_arr];


infoX10A = load('aggregated results/fA10_info_calcs.mat');
infoX10S = load('aggregated results/fS10_info_calcs.mat');

infoX10 = [infoX10S.info_ds_arr; infoX10A.info_ds_arr];
%% Figure 1 and 2 numbers

% More generally, a large majority (mean: 78%, SE: 16 % across N = 240 sets) 
% of linear readouts will have stimulus predictive information less than 1 bit/second. 
setOf4_stimPredIperS = cellfun(@(x) x.X_stim.percstim_predI_vt(:, ...
    x.X_stim.tau_futures == 1)*60, infoX4, 'uniformoutput', false);
full_PercSetCalculated = cellfun(@(x) length(x) == 149, setOf4_stimPredIperS);
setOf4_stimPredIperS = setOf4_stimPredIperS(full_PercSetCalculated);

mean_less_than_1bps = cellfun(@(x) mean(x < 1), setOf4_stimPredIperS);
display(['SETS OF 4: mean % of readouts under 1 bit/s: ' num2str(mean(mean_less_than_1bps)*100, '%1.0f') ...
    ', SE: ' num2str(std(mean_less_than_1bps)*100, '%1.0f')])


% Across many sets (N = 240, all cells represented in at least one set) of 
% four cells, the stimulus predictive information of a binary word is highly 
% correlated with the internal predictive information (Fig. 1G, r = 0.59, N = 240 sets), 
word4_stimPredI = cellfun(@(x) x.X_stim.wordstim_predI_vt(x.X_stim.tau_futures == 1), ...
    infoX4(full_PercSetCalculated));
word_intnrbPredI =  cellfun(@(x) x.X_stim.wordword_nrbpredI_vt(x.X_stim.tau_futures == 1), ...
    infoX4(full_PercSetCalculated));
r_stim_nrbin = OffDiag(corrcoef(word4_stimPredI, word_intnrbPredI));
display(['SETS OF 4: correlation btw stim and internal pred-I (nrb movie) : ' num2str(r_stim_nrbin, '%1.2f') ...
    ', N = ' num2str(length(word4_stimPredI))])

% Across sets for data recorded during the moving-bar movie, correlation 
% between internal and stimulus predictive information of linear readouts 
% is high (linear correlation: 0.71 +/- 0.32 (mean +/- SE across sets); 
% N=215 of 240 sets with significant correlation (p < 0.01)). 
setOf4_intnrbPredIperS = cellfun(@(x) x.X_stim.percword_nrbpredI_vt(:, ...
    x.X_stim.tau_futures == 1)*60, infoX4(full_PercSetCalculated), 'uniformoutput', false);
[readout_xc_nrb, readout_p_xc_nrb] = cellfun(@(x, y) corrcoef(x,y), ...
    setOf4_intnrbPredIperS, setOf4_stimPredIperS, 'uniformoutput', false);
readout_xc_nrb = cellfun(@(x) x(1,2), readout_xc_nrb);
readout_p_xc_nrb = cellfun(@(x) x(1,2), readout_p_xc_nrb);
display(['SETS OF 4: cross-set mean correlation (SE) is ' num2str(mean(readout_xc_nrb), '%1.2f') ...
    ' (' num2str(std(readout_xc_nrb), '%1.2f') '), with p < .01 for ' ...
    num2str(sum(readout_p_xc_nrb < 0.01)) ' of ' num2str(length(readout_p_xc_nrb)) ' sets'])



% The correlation between internal predictive information on the fish movie 
% and on the moving-bar movie was high (Fig 2B; r = 0.80 for sets of 4 cells; 
% = 0.87 for sets of 10 cells; p < 1e-8), as was the correlation between the 
% stimulus predictive information and the internal fish movie predictive 
% information (Fig. 2C; r= 0.55 for sets of 4 cells; 
% = 0.78 for sets of 10 cells; p < 1e-7). 
word_intPredI = cellfun(@(x) x.X_ww.word_pI(x.X_ww.tau_futures == 1), ...
    infoX4(full_PercSetCalculated));
r_stim_internal = OffDiag(corrcoef(word4_stimPredI, word_intPredI));
r_nrbint_internal = OffDiag(corrcoef(word_intnrbPredI, word_intPredI));
display(['SETS OF 4: correlation btw internal pred-I (nrb vs fish movie) : ' num2str(r_nrbint_internal, '%1.2f') ...
    ', N = ' num2str(length(word_intPredI))])
display(['SETS OF 4: correlation btw stim and internal pred-I (fish movie) : ' num2str(r_stim_internal, '%1.2f') ...
    ', N = ' num2str(length(word4_stimPredI))])

% sets of 10
setOf10_stimPredIperS = cellfun(@(x) x.X_stim.percstim_predI_vt(:, ...
    x.X_stim.tau_futures == 1)*60, infoX10, 'uniformoutput', false);
setOf10_intnrbPredIperS = cellfun(@(x) x.X_stim.percword_nrbpredI_vt(:, ...
    x.X_stim.tau_futures == 1)*60, infoX10, 'uniformoutput', false);
setOf10_predIperS = cellfun(@(x) x.X_pw.pI_perc(:, x.X_pw.tau_futures == 1)*60, ...
    infoX10, 'uniformoutput', false);
word10_intPredI = cellfun(@(x) x.X_ww.word_pI(x.X_ww.tau_futures == 1), ...
    infoX10);
word10_stimPredI = cellfun(@(x) x.X_stim.wordstim_predI_vt(x.X_stim.tau_futures == 1), ...
    infoX10);
word10_intnrbPredI =  cellfun(@(x) x.X_stim.wordword_nrbpredI_vt(x.X_stim.tau_futures == 1), ...
    infoX10);
r10_stim_internal = OffDiag(corrcoef(word10_stimPredI, word10_intnrbPredI));
r10_nrbint_internal = OffDiag(corrcoef(word10_intnrbPredI, word10_intPredI));
display(['SETS OF 10: correlation btw internal pred-I (nrb vs fish movie) : ' num2str(r10_nrbint_internal, '%1.2f') ...
    ', N = ' num2str(length(word10_stimPredI))])
display(['SETS OF 10: correlation btw stim and internal pred-I (fish movie) : ' num2str(r10_stim_internal, '%1.2f') ...
    ', N = ' num2str(length(word10_stimPredI))])

% Similar to the readouts of stimulus predictive information during the 
% moving-bar movie, most readouts of internal predictive information had 
% less than 1 bit/second predictive information (72% +/- 23% SE of all 
% readouts across sets), but a small fraction were highly informative 
% (Fig. 2D, single set example). For this particular set, the predictive 
% information of readouts of internal information are significantly correlated 
% with readouts of external information (Fig. 1E; r = 0.33, p < 1e-4). 
setOf4_predIperS = cellfun(@(x) x.X_pw.pI_perc(:, x.X_pw.tau_futures == 1)*60, ...
    infoX4(full_PercSetCalculated), 'uniformoutput', false);
mean_less_than_1bps = cellfun(@(x) mean(x < 1), setOf4_predIperS);
display(['SETS OF 4: mean % of internal pred-I readouts under 1 bit/s: ' num2str(mean(mean_less_than_1bps)*100, '%1.0f') ...
    ', SE: ' num2str(std(mean_less_than_1bps)*100, '%1.0f')])
ex_info = load('code for respository/information calculation results/manyICs_g4/set668_fA4.mat');
perc_pI_668 = 60*ex_info.X_pw.pI_perc(:, 2);
perc_stimpI_668 = 60*ex_info.X_stim.percstim_predI_vt(:, ex_info.X_stim.tau_futures == 1);
[xc_668, p_xc_668] = corrcoef(perc_pI_668, perc_stimpI_668);
display(['fA4(668): internal and stim corr coef: ' num2str(xc_668(1,2), '%1.2f') ', p = ' num2str(p_xc_668(1,2), '%1.2d')])


% Across all sampled sets of 4 cells, the linear correlation between 
% % internal predictive information and stimulus predictive information 
% was even higher than it was for the sample set (Fig. 2E, sets of 4: 
% mean: 0.73, SE: 0.24 across 240 sets). Typical correlations between 
% internal and external predictive information readouts are lower but 
% still highly significant for sets of 10 cells (correlation: 0.42, SE: 0.33).
[readout_xc_predI, readout_p_xc_predI] = cellfun(@(x, y) corrcoef(x,y), ...
    setOf4_predIperS, setOf4_stimPredIperS, 'uniformoutput', false);
readout_xc_predI = cellfun(@(x) x(1,2), readout_xc_predI);
readout_p_xc_predI = cellfun(@(x) x(1, 2), readout_p_xc_predI);
display(['SETS OF 4: cross-set mean correlation (stim vs. internal-fish readout) (SE) is ' num2str(mean(readout_xc_predI), '%1.2f') ...
    ' (' num2str(std(readout_xc_predI), '%1.2f') '), with p < .01 for ' ...
    num2str(sum(readout_p_xc_predI < 0.01)) ' of ' num2str(length(readout_xc_predI)) ' sets'])

[readout10_xc_predI, readout10_p_xc_predI] = cellfun(@(x, y) corrcoef(x,y), ...
    setOf10_predIperS, setOf10_stimPredIperS, 'uniformoutput', false);
readout10_xc_predI = cellfun(@(x) x(1,2), readout10_xc_predI);
readout10_p_xc_predI = cellfun(@(x) x(1, 2), readout10_p_xc_predI);
display(['SETS OF 10: cross-set mean correlation (stim vs. internal-fish readout) (SE) is ' num2str(mean(readout10_xc_predI), '%1.2f') ...
    ' (' num2str(std(readout10_xc_predI), '%1.2f') '), with p < .01 for ' ...
    num2str(sum(readout10_p_xc_predI < 0.01)) ' of ' num2str(length(readout10_xc_predI)) ' sets'])

%% Fig 3 numbers

% For small sets (4 cells), readouts learned under the pair rule conveyed 
% near-optimal predictive information (Fig. 3D). The average percent of the 
% optimal predictive information learned was 87% (13% SE across N = 360 sets). 
% Learned readouts had firing rates distributed across the range of readout 
% firing rates (Fig. 3E), with an average firing rate (mean: 3.4 Hz, SE: 1.2 Hz, 
% N = 360 sets) that is 68% of the maximum firing rate. We quantified the similarity 
% of learned readout rules to the optimal readout rule (Fig. 3C) using the 
% correlation between learned and optimal rules across input words, weighted 
% by the frequency of input words (Fig. 3F). Learned rules were highly similar, 
% but not identical, to optimal rules (Fig. 3F; learned mean: 0.71). 
% Although they did not precisely match the optimal readout structure, readouts 
% learned under the pair rule were efficient at representing predictive information. 

i_s4 = 1;
i_lt = 1;
averagePrctLearned = cellfun(@(x) mean(x), set_processedRes{i_s4}.prct_learned{i_lt}(full_PercSetCalculated));
display(['SETS OF 4: average percent learned: ' num2str(mean(averagePrctLearned), '%1.0f') ...
    ' (SE: ' num2str(std(averagePrctLearned), '%1.0f') ', N = ' num2str(length(averagePrctLearned)) ')' ])
averageFiringRate = cellfun(@(x, y) 60*max(y)*mean(x)/100, ...
    set_processedRes{i_s4}.prct_maxFR{i_lt}(full_PercSetCalculated), ...
    set_processedRes{i_s4}.all_lt_perc_fr{i_lt}(full_PercSetCalculated));
averageFiringRateFRAC = cellfun(@(x, y) mean(x), set_processedRes{i_s4}.prct_maxFR{i_lt}(full_PercSetCalculated));
display(['SETS OF 4: average FR: ' num2str(mean(averageFiringRate), '%1.2f') ...
    ' (SE: ' num2str(std(averageFiringRate), '%1.2f') ', N = ' num2str(length(averageFiringRate)) ')' ])
display(['SETS OF 4: average % FR: ' num2str(mean(averageFiringRateFRAC), '%1.0f') ...
    ' (SE: ' num2str(std(averageFiringRateFRAC), '%1.0f') ', N = ' num2str(length(averageFiringRateFRAC)) ')' ])

averageRuleSims = cellfun(@(x) mean(x), set_processedRes{i_s4}.rule_sims{i_lt}(full_PercSetCalculated));
display(['SETS OF 4: average rule similarity: ' num2str(mean(averageRuleSims), '%1.2f') ...
    ' (SE: ' num2str(std(averageRuleSims), '%1.2f') ', N = ' num2str(length(averageRuleSims)) ')' ])

%% Fig 4 numbers
% Compared to the optimal sampled readouts, we observed a small decrease in 
% readout efficiency (Fig. 4B), from 87% on average for groups of 4, to 80% 
% (SE over groups, 10%, N = 100 groups of 7) and 
% 79% (SE over groups, 8%, N = 100 groups of 10). Learned readouts for 
% these cells were much less similar to the optimal readout than the readouts 
% for groups of four were (0.61, 0.60 for groups of 7 and 10, respectively; 
% Fig. 4A). While there is a relationship between structural similarity to 
% the optimal rule and efficiency of the readout, many readouts have a high 
% degree of predictive information efficiency with low structural similarity 
% to the optimal rule (Fig. 4C-D; similar for sets of 4 (not shown)). Thus, 
% for sets of 4-10 cells, efficient readouts of predictive information could 
% be learned without finding the exact structure of the optimal readout. 


i_s7 = 2;
i_lt = 1;
averagePrctLearned = cellfun(@(x) mean(x), set_processedRes{i_s7}.prct_learned{i_lt});
display(['SETS OF 7: average percent learned: ' num2str(mean(averagePrctLearned), '%1.0f') ...
    ' (SE: ' num2str(std(averagePrctLearned), '%1.0f') ', N = ' num2str(length(averagePrctLearned)) ')' ])
averageFiringRate = cellfun(@(x, y) 60*max(y)*mean(x)/100, ...
    set_processedRes{i_s7}.prct_maxFR{i_lt}, ...
    set_processedRes{i_s7}.all_lt_perc_fr{i_lt});
averageFiringRateFRAC = cellfun(@(x, y) mean(x), set_processedRes{i_s7}.prct_maxFR{i_lt});
display(['SETS OF 7: average FR: ' num2str(mean(averageFiringRate), '%1.2f') ...
    ' (SE: ' num2str(std(averageFiringRate), '%1.2f') ', N = ' num2str(length(averageFiringRate)) ')' ])
display(['SETS OF 7: average % FR: ' num2str(mean(averageFiringRateFRAC), '%1.0f') ...
    ' (SE: ' num2str(std(averageFiringRateFRAC), '%1.0f') ', N = ' num2str(length(averageFiringRateFRAC)) ')' ])

averageRuleSims = cellfun(@(x) mean(x), set_processedRes{i_s7}.rule_sims{i_lt});
display(['SETS OF 7: average rule similarity: ' num2str(mean(averageRuleSims), '%1.2f') ...
    ' (SE: ' num2str(std(averageRuleSims), '%1.2f') ', N = ' num2str(length(averageRuleSims)) ')' ])

i_s10 = 3;
i_lt = 1;
averagePrctLearned = cellfun(@(x) mean(x), set_processedRes{i_s10}.prct_learned{i_lt});
display(['SETS OF 10: average percent learned: ' num2str(mean(averagePrctLearned), '%1.0f') ...
    ' (SE: ' num2str(std(averagePrctLearned), '%1.0f') ', N = ' num2str(length(averagePrctLearned)) ')' ])
averageFiringRate = cellfun(@(x, y) 60*max(y)*mean(x)/100, ...
    set_processedRes{i_s10}.prct_maxFR{i_lt}, ...
    set_processedRes{i_s10}.all_lt_perc_fr{i_lt});
averageFiringRateFRAC = cellfun(@(x, y) mean(x), set_processedRes{i_s10}.prct_maxFR{i_lt});
display(['SETS OF 10: average FR: ' num2str(mean(averageFiringRate), '%1.2f') ...
    ' (SE: ' num2str(std(averageFiringRate), '%1.2f') ', N = ' num2str(length(averageFiringRate)) ')' ])
display(['SETS OF 10: average % FR: ' num2str(mean(averageFiringRateFRAC), '%1.0f') ...
    ' (SE: ' num2str(std(averageFiringRateFRAC), '%1.0f') ', N = ' num2str(length(averageFiringRateFRAC)) ')' ])

averageRuleSims = cellfun(@(x) mean(x), set_processedRes{i_s10}.rule_sims{i_lt});
display(['SETS OF 10: average rule similarity: ' num2str(mean(averageRuleSims), '%1.2f') ...
    ' (SE: ' num2str(std(averageRuleSims), '%1.2f') ', N = ' num2str(length(averageRuleSims)) ')' ])


%% Figure 5: stimulus info
% compute predI/spike
learned_pIps = cellfun(@(xS, yS) cellfun(@(x,y) x(:, tau_future_vals == 1)./y, ...
    xS, yS, 'uniformoutput', false), ...
    learned_stimPredI, learned_stimFR, 'uniformoutput', false);
%%
optimal_pIps = cellfun(@(xS, yS) cellfun(@(x,y) x(:, tau_future_vals == 1)./y, ...
    xS, yS, 'uniformoutput', false), ...
    optimal_stimPredI, optimal_stimFR, 'uniformoutput', false);
%%
total_pIps = cellfun(@(xS, yS) cellfun(@(x,y) x(:, tau_future_vals == 1)./y, ...
    xS, yS, 'uniformoutput', false), ...
    word_stimPredI, word_stimFR, 'uniformoutput', false);
%%
% get predI
learned_stimpI = cellfun(@(xS) cellfun(@(x) x(:, tau_future_vals == 1), ...
    xS, 'uniformoutput', false), ...
    learned_stimPredI, 'uniformoutput', false);
optimal_stimpI = cellfun(@(xS) cellfun(@(x) x(:, tau_future_vals == 1), ...
    xS, 'uniformoutput', false), optimal_stimPredI, 'uniformoutput', false);
total_stimpI = cellfun(@(xS) cellfun(@(x) x(:, tau_future_vals == 1), ...
    xS, 'uniformoutput', false), word_stimPredI, 'uniformoutput', false);


fracTotalPredI = cellfun(@(xS, yS) cellfun(@(x,y) x(:, tau_future_vals == 1)/y(tau_future_vals == 1), ...
    xS, yS, 'uniformoutput', false), ...
    learned_stimPredI, word_stimPredI, 'uniformoutput', false);

fracOptPredI = cellfun(@(xS, yS) cellfun(@(x, y) ...
    x(:, tau_future_vals == 1)./y(:, tau_future_vals == 1), ...
    xS, yS, 'uniformoutput', false), ...
    learned_stimPredI, optimal_stimPredI, 'uniformoutput', false);

fracTotalPredIps = cellfun(@(xS, yS) cellfun(@(x,y) x/y, ...
    xS, yS, 'uniformoutput', false), ...
    learned_pIps, total_pIps, 'uniformoutput', false);

fracOptPredIps = cellfun(@(xS, yS) cellfun(@(x, y) x./y, ...
    xS, yS, 'uniformoutput', false), ...
    learned_pIps, optimal_pIps, 'uniformoutput', false);
% find the index of the learned readout corresponding to the median
% internal pred-I
learned_internalPredI = cellfun(@(x)  cellfun(@(y) median(y(1:end-1)), x.learned_pI{1}), ...
    set_processedRes, 'uniformoutput', false);

[~, learned_internalPredI_ICord] = cellfun(@(x)  cellfun(@(y) sort(y), ...
        x.learned_pI{1}, 'uniformoutput', false), ...
        set_processedRes, 'uniformoutput', false);

benchmark_prctl = 50;

learned_internalPredI_ICind = cellfun(@(x)  cellfun(@(y) y(round(length(y)*benchmark_prctl/100)), ...
    x, 'uniformoutput', false), ...
    learned_internalPredI_ICord, 'uniformoutput', false);

had200_ICs = cellfun(@(x) cellfun(@(y) length(y) == 200, x), learned_pIps, ...
    'uniformoutput', false);
had2k_ICs = cellfun(@(x) ~x, had200_ICs, 'uniformoutput', false);
had200_ICs{1} = had200_ICs{1} & full_PercSetCalculated;
median_fTPI = cellfun(@(x, inds, ind200) cellfun(@(y, y_inds) y(y_inds), x(ind200), inds(ind200)), ...
    fracTotalPredI, learned_internalPredI_ICind, had200_ICs, 'uniformoutput', false);
median_fOPI = cellfun(@(x, inds, ind200) cellfun(@(y, y_inds) y(y_inds), x(ind200), inds(ind200)), ...
    fracOptPredI, learned_internalPredI_ICind, had200_ICs, 'uniformoutput', false);


median_fTPIps = cellfun(@(x, inds, ind200) cellfun(@(y, y_inds) y(y_inds), x(ind200), inds(ind200)), ...
    fracTotalPredIps, learned_internalPredI_ICind, had200_ICs, 'uniformoutput', false);
median_fOPIps = cellfun(@(x, inds, ind200) cellfun(@(y, y_inds) y(y_inds), x(ind200), inds(ind200)), ...
    fracOptPredIps, learned_internalPredI_ICind, had200_ICs, 'uniformoutput', false);

y_T = cellfun(@(x) nanmean(min(x, 5)), median_fTPI);
dy_T = cellfun(@(x) nanstd(x)/sqrt(length(x)-1), median_fTPI);

y_O = cellfun(@(x) nanmean(min(x, 5)), median_fOPI);
dy_O = cellfun(@(x) nanstd(x)/sqrt(length(x)-1), median_fOPI);



y_Tps = cellfun(@(x) nanmean(min(x, 5)), median_fTPIps);
dy_Tps = cellfun(@(x) nanstd(x)/sqrt(length(x)-1), median_fTPIps);

y_Ops = cellfun(@(x) nanmean(min(x, 5)), median_fOPIps);
dy_Ops = cellfun(@(x) nanstd(x)/sqrt(length(x)-1), median_fOPIps);


display([num2str(benchmark_prctl) 'percentile learned perceptron (benchmarked on internal info):' ...
    'efficiency (pI/spike) relative to optimal is ' num2str(y_Ops', '%1.2f, ') '(mean)' ...
    num2str(dy_Ops', '%1.2f, ') '(SE)'])
display([num2str(benchmark_prctl) 'percentile learned perceptron (benchmarked on internal info):' ...
    'efficiency (pI/spike) relative to cell set is ' num2str(y_Tps', '%1.2f, ') '(mean)' ...
    num2str(dy_Tps', '%1.2f, ') '(SE)'])

display([num2str(benchmark_prctl) 'percentile learned perceptron (benchmarked on internal info):' ...
    'info rate (pI/s) relative to optimal is ' num2str(y_O', '%1.2f, ') '(mean)' ...
    num2str(dy_O', '%1.2f, ') '(SE)'])
display([num2str(benchmark_prctl) 'percentile learned perceptron (benchmarked on internal info):' ...
    'info rate (pI/s) relative to cell set is ' num2str(y_T', '%1.2f, ') '(mean)' ...
    num2str(dy_T', '%1.2f, ') '(SE)'])


%% for each group, how many learned perceptrons were more efficient than the optimal 
% perceptron? 

prct_thr = 0.95;

numICs_overPctThr_fracTotalStimPI = cellfun(@(x, y) cellfun(@(x2) sum(x2 >= prct_thr), ...
    x(y)), fracTotalPredI, had200_ICs, 'uniformoutput', false);
numICs_overPctThr_fracTotalStimPIps = cellfun(@(x, y) cellfun(@(x2) sum(x2 >= prct_thr), ...
    x(y)), fracTotalPredIps, had200_ICs, 'uniformoutput', false);

numICs_overPctThr_fracOptStimPI = cellfun(@(x, y) cellfun(@(x2) sum(x2 >= prct_thr), ...
    x(y)), fracOptPredI, had200_ICs, 'uniformoutput', false);
numICs_overPctThr_fracOptStimPIps = cellfun(@(x, y) cellfun(@(x2) sum(x2 >= prct_thr), ...
    x(y)), fracOptPredIps, had200_ICs, 'uniformoutput', false);


display([num2str(cellfun(@(x) sum(x > 0), numICs_overPctThr_fracTotalStimPI)', '%1.0f, ') ' of ' ...
    'sets had ICs leading to readouts with > ' num2str(100*prct_thr, '%1.0f') ...
    '% of the total stimulus pred-I'])
display([num2str(cellfun(@(x) sum(x > 0), numICs_overPctThr_fracTotalStimPIps)', '%1.0f, ') ' of ' ...
    'sets had ICs leading to readouts with > ' num2str(100*prct_thr, '%1.0f') ...
    '% of the total stimulus pred-I/spike'])

display([num2str(cellfun(@(x) sum(x > 0), numICs_overPctThr_fracOptStimPI)', '%1.0f, ') ' of ' ...
    'sets had ICs leading to readouts with > ' num2str(100*prct_thr, '%1.0f') ...
    '% of the optimal stimulus pred-I'])
display([num2str(cellfun(@(x) sum(x > 0), numICs_overPctThr_fracOptStimPIps)', '%1.0f, ') ' of ' ...
    'sets had ICs leading to readouts with > ' num2str(100*prct_thr, '%1.0f') ...
    '% of the optimal stimulus pred-I/spike'])

display([num2str(cellfun(@(x) sum(x), had200_ICs)', '%1.0f, ') ' sets in each group'])

display('******** 2K IC results')

fracICs_overPctThr_fracTotalStimPI = cellfun(@(x, y) cellfun(@(x2) mean(x2 >= prct_thr), ...
    x(y)), fracTotalPredI, had2k_ICs, 'uniformoutput', false);
fracICs_overPctThr_fracTotalStimPIps = cellfun(@(x, y) cellfun(@(x2) mean(x2 >= prct_thr), ...
    x(y)), fracTotalPredIps, had2k_ICs, 'uniformoutput', false);

fracICs_overPctThr_fracOptStimPI = cellfun(@(x, y) cellfun(@(x2) mean(x2 >= prct_thr), ...
    x(y)), fracOptPredI, had2k_ICs, 'uniformoutput', false);
fracICs_overPctThr_fracOptStimPIps = cellfun(@(x, y) cellfun(@(x2) mean(x2 >= prct_thr), ...
    x(y)), fracOptPredIps, had2k_ICs, 'uniformoutput', false);

extract_fun = @(x2) cell2mat(cellfun(@(x) [mean(x) std(x)], x2(1), 'uniformoutput', false));

display([num2str(extract_fun(fracICs_overPctThr_fracTotalStimPI), '%1.2f, ') ' of ' ...
    ' ICs leading to readouts with > ' num2str(100*prct_thr, '%1.0f') ...
    '% of the total stimulus pred-I'])
display([num2str(extract_fun(fracICs_overPctThr_fracTotalStimPIps), '%1.2f, ') ' of ' ...
    ' ICs leading to readouts with > ' num2str(100*prct_thr, '%1.0f') ...
    '% of the total stimulus pred-I/spike'])

display([num2str(extract_fun(fracICs_overPctThr_fracOptStimPI), '%1.2f, ') ' of ' ...
    ' ICs leading to readouts with > ' num2str(100*prct_thr, '%1.0f') ...
    '% of the optimal stimulus pred-I'])
display([num2str(extract_fun(fracICs_overPctThr_fracOptStimPIps), '%1.2f, ') ' of ' ...
    ' ICs leading to readouts with > ' num2str(100*prct_thr, '%1.0f') ...
    '% of the optimal stimulus pred-I/spike'])

display([num2str(cellfun(@(x) sum(x), had2k_ICs)', '%1.0f, ') ' sets in each group'])

