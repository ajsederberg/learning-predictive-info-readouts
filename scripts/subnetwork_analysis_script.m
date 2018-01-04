% Subnetwork analysis script



%% load aggregated info calculations
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


%% find the set of 4 with the median readout-readout correlation
[~, ind_order] = sort(readout_xc_predI);
set_ind = find(readout_xc_predI == median(readout_xc_predI(1:end-1)));
all_inds = find(full_PercSetCalculated);

absolute_set_ordered_indices = all_inds(ind_order);
ww_info_ordered = word_intPredI(ind_order);

possible_inds = find(ww_info_ordered > 0.125 & ww_info_ordered < 0.135);

%% find the set of 10 with the median readout-readout correlation
[~, ind_order10] = sort(readout10_xc_predI);

absolute_set_ordered_indices10 = (ind_order10);
ww10_info_ordered = word10_intPredI(ind_order10);

possible_inds10 = find(ww10_info_ordered > 0.245 & ww10_info_ordered < 0.255)

%% make a simple scatter plot of readout pred-I for the stimulus info and the internal info
%  try 123, 70, 150

for i_mid10 = [9 197] % 99 125 197]
    absolute_set_index10 = absolute_set_ordered_indices10(i_mid10);
    ex_info = infoX10{absolute_set_index10};
    y = 60*ex_info.X_pw.pI_perc(:, 2);
    x = 60*ex_info.X_stim.percstim_predI_vt(:, ex_info.X_stim.tau_futures == 1);

    ww_info = num2str(ex_info.X_ww.word_pI(1), '%1.2f');
    ww_stiminfo = num2str(ex_info.X_stim.wordstim_predI_vt(ex_info.X_stim.tau_futures == 1), '%1.2f');
    mb = polyfit(x, y, 1);

    makeMyFigure(8,8);
    plot(x, y, 'o')
    hold on
    x_vec = linspace(0, max(x), 11);
    plot(x_vec, polyval(mb, x_vec), 'k')
    [xc , p_xc] = corrcoef(x, y);
    axis tight
    axis square
    xlabel('stimulus pred-I (bit/s)')
    ylabel('readout-word pred-I (bit/s)')
    title({['r = ' num2str(xc(1,2), '%1.2f')]; ['set ' num2str(ex_info.set_ID)]; ...
        ['word-word info: ' ww_info ' bits/bin']; ...
        ['word-stim info: ' ww_stiminfo ' bits/bin']})
    set(gca, 'box', 'off', 'color', 'none')
%     print(gcf, '-dpdf', ['scripts/Figure Plots/readoutReadoutScatter_g10_set' num2str(ex_info.set_ID)])
end

%% color point by perceptron weight
num_colors = 10;
cmap = (jet(num_colors));
for i_mid10 = 99 %[9 99 125 197]
    absolute_set_index10 = absolute_set_ordered_indices10(i_mid10);
    ex_info = infoX10{absolute_set_index10};
    y = 60*ex_info.X_pw.pI_perc(:, 2);
%     y = 60*ex_info.X_stim.percword_nrbpredI_vt(:, ex_info.X_stim.tau_futures == 1);

    x = 60*ex_info.X_stim.percstim_predI_vt(:, ex_info.X_stim.tau_futures == 1);

    ww_info = num2str(ex_info.X_ww.word_pI(1), '%1.2f');
    mb = polyfit(x, y, 1);

    
    for i_wt = 1:10 %[1 4 8 9] %
        color_inds = ceil(num_colors*(eps + ex_info.X_stim.perc_ws(:, i_wt))/1.11);

        makeMyFigure(8,8);
        hold on

        for i_cmap = 1:num_colors
            inds = color_inds == i_cmap;
            plot(x(inds), y(inds), 'o', 'color', cmap(i_cmap, :))
%             x_vec = linspace(0, max(x), 11);
%             plot(x_vec, polyval(mb, x_vec), 'k')

        end
        [xc , p_xc] = corrcoef(x, y);

        axis tight
        axis square
        xlabel('stimulus pred-I (bit/s)')
        ylabel('readout-word pred-I (bit/s)')
        title({['r = ' num2str(xc(1,2), '%1.2f')]; ['set ' num2str(ex_info.set_ID) ', w_{' num2str(i_wt) '}']; ...
            ['word-word info: ' ww_info ' bits/bin']})
        set(gca, 'box', 'off', 'color', 'none')
%     print(gcf, '-dpdf', ['scripts/Figure Plots/readoutReadoutScatter_g10_set' num2str(ex_info.set_ID)])
    end
end
%% make a simple scatter plot of readout pred-I for the stimulus info and the internal info
%  try 123, 70, 150

for i_mid = [119 ] %119 154]
    absolute_set_index = absolute_set_ordered_indices(i_mid);
    ex_info = infoX4{absolute_set_index};
    y = 60*ex_info.X_pw.pI_perc(:, 2);
    x = 60*ex_info.X_stim.percstim_predI_vt(:, ex_info.X_stim.tau_futures == 1);

    ww_info = num2str(ex_info.X_ww.word_pI(1), '%1.2f');
    mb = polyfit(x, y, 1);

    makeMyFigure(8,8);
    plot(x, y, 'o')
    hold on
    x_vec = linspace(0, max(x), 11);
    plot(x_vec, polyval(mb, x_vec), 'k')
    [xc , p_xc] = corrcoef(x, y);
    axis tight
    axis square
    xlabel('stimulus pred-I (bit/s)')
    ylabel('readout-word pred-I (bit/s)')
    title({['r = ' num2str(xc(1,2), '%1.2f')]; ['set ' num2str(ex_info.set_ID)]; ...
        ['word-word info: ' ww_info ' bits/bin']})
    set(gca, 'box', 'off', 'color', 'none')
%     print(gcf, '-dpdf', ['scripts/Figure Plots/readoutReadoutScatter_set' num2str(ex_info.set_ID)])
end


%% color point by perceptron weight, plot with STA's alongside
num_colors = 5;
cmap = (jet(num_colors));

set_size = 4;
num_cols = 3;

time_range_bins = -15:1:5;



makeMyFigure(num_cols*15, set_size*12);
for i_mid = possible_inds' %[9 99 125 197]
    clf;
    if set_size == 4
        absolute_set_index = absolute_set_ordered_indices(i_mid);
        ex_info = infoX4{absolute_set_index};
    else
        absolute_set_index = absolute_set_ordered_indices10(i_mid);        
        ex_info = infoX10{absolute_set_index};
    end
    
    ww_info = num2str(ex_info.X_ww.word_pI(1), '%1.2f');
    wstim_info = num2str(ex_info.X_stim.wordstim_predI_vt(ex_info.X_stim.tau_futures == 1), '%1.2f');
       
    y = ex_info.X_pw.pI_perc(:, 2)/str2double(ww_info);
    x = ex_info.X_stim.percstim_predI_vt(:, ex_info.X_stim.tau_futures == 1)/str2double(wstim_info);

    mb = polyfit(x, y, 1);

    cell_set = ex_info.X_stim.cell_IDs;
    comp_type = 'CDM_pv';
    learned_ws = eye(set_size);
    [nrb_stim, nrb_data] = barmovie_activity_readout_function(learned_ws, cell_set, comp_type);
    [sta_velocity, v_sta_arr, velocity_bins, time_range_bins] = computeVelocitySTA_NRBdata(nrb_stim, nrb_data, time_range_bins);
    [sta_position, p_sta_arr, position_bins, time_range_bins] = computePositionSTA_NRBdata(nrb_stim, nrb_data, time_range_bins);

    wt_info_centroid = zeros(set_size, 2);
    for i_wt = 1:set_size %[1 4 8 9] %
        color_inds = ceil(num_colors*(eps + ex_info.X_stim.perc_ws(:, i_wt))/1.11);

        subplot(set_size, num_cols, num_cols*(i_wt - 1) + 1)
        hold on

        for i_cmap = 1:num_colors
            inds = color_inds == i_cmap;
            plot(x(inds), y(inds), 'o', 'color', cmap(i_cmap, :))
%             x_vec = linspace(0, max(x), 11);
%             plot(x_vec, polyval(mb, x_vec), 'k')

            if i_cmap == num_colors
                wt_info_centroid(i_wt, :) = mean([x(inds), y(inds)]) - mean([x(~inds), y(~inds)]);
                
                xInOut = [mean(x(inds)) mean(x(~inds))];
                yInOut = [mean(y(inds)) mean(y(~inds))];
                xInOut = [xInOut(1) max(xInOut) xInOut(2)];
                yInOut = [yInOut(1) min(yInOut) yInOut(2)];
                
                plot(xInOut, yInOut,...
                    'ko-', 'color', cmap(i_cmap, :), 'linewidth', 1.5)

            end

        end
        [xc , p_xc] = corrcoef(x, y);

        axis tight
        axis square
        xlabel('stimulus pred-I (bit/s)')
        ylabel('readout-word pred-I (bit/s)')
        title({['r = ' num2str(xc(1,2), '%1.2f')]; ['set ' num2str(ex_info.set_ID) ', cell ' num2str(cell_set(i_wt)) ' weight']; ...
            ['fish word-word info: ' ww_info ' bits/bin'];['word-stim info: ' wstim_info] })
        set(gca, 'box', 'off', 'color', 'none')
%     print(gcf, '-dpdf', ['scripts/Figure Plots/readoutReadoutScatter_g10_set' num2str(ex_info.set_ID)])

        % plot spike-triggered position and velocity distributions
        subplot(set_size, num_cols, num_cols*(i_wt - 1) + 2)
        imagesc(time_range_bins*1000/60, position_bins, sta_position{i_wt})
        hold on
        plot(time_range_bins*1000/60, p_sta_arr{i_wt, 1}, 'r', 'linewidth', 1.5)
        plot(time_range_bins*1000/60, p_sta_arr{i_wt, 1} + [p_sta_arr{i_wt, 2}; -p_sta_arr{i_wt, 2}] , 'r', 'linewidth', 1)
        axis square
        colorbar
        set(gca, 'ydir', 'normal')
        title(['cell ' num2str(cell_set(i_wt)) ' ST position dist'])
        
        subplot(set_size, num_cols, num_cols*(i_wt - 1) + 3)
        imagesc(time_range_bins*1000/60, velocity_bins, sta_velocity{i_wt})
        hold on
        plot(time_range_bins*1000/60, v_sta_arr{i_wt, 1}, 'r', 'linewidth', 1.5)
        plot(time_range_bins*1000/60, v_sta_arr{i_wt, 1} + [v_sta_arr{i_wt, 2}; -v_sta_arr{i_wt, 2}] , 'r', 'linewidth', 1)
        axis square
        colorbar      
        set(gca, 'ydir', 'normal')
        title(['cell ' num2str(cell_set(i_wt)) ' ST velocity dist'])
    end
    print(gcf, '-dpng', ['scripts/Figure Plots/revisions/readoutscatter_withSTA_g' ...
        num2str(set_size) '_set' num2str(ex_info.set_ID)])
    
    
    %% makea  second figure, showing summary information
    makeMyFigure(30, 15);
    subplot(1, 3, 1)
    plot(1:2, wt_info_centroid, 'linewidth', 1.5)
    ylabel('cell weight impact on pred-I')
    set(gca, 'xtick', 1:2, 'xticklabel', {'external', 'internal'}, 'color', 'none')
    axis square
    xlim([0.5 2.5])
    legend(num2str((1:set_size)', 'cell %1.0f'), 'location', 'best')
    
    cell_map = lines(set_size);
    subplot(1, 3, 2)
    hold on
    for i_cell = 1:set_size
        plot(time_range_bins*1000/60, p_sta_arr{i_cell, 1}, ...
            'color', cell_map(i_cell, :),  'linewidth', 1.5)
%         plot(time_range_bins*1000/60, ...
%             p_sta_arr{i_cell, 1} + [p_sta_arr{i_cell, 2}; -p_sta_arr{i_cell, 2}], ...
%             'color', cell_map(i_cell, :), 'linewidth', 1)
    end
    ylabel('spike-triggered position')
    set(gca, 'color', 'none')
    axis square
    xlim(time_range_bins([1 end])*1000/60)
    
    subplot(1, 3, 3)
    hold on
    for i_cell = 1:set_size
        plot(time_range_bins*1000/60, v_sta_arr{i_cell, 1}, ...
            'color', cell_map(i_cell, :),  'linewidth', 1.5)
%         plot(time_range_bins*1000/60, ...
%             v_sta_arr{i_cell, 1} + [v_sta_arr{i_cell, 2}; -v_sta_arr{i_cell, 2}], ...
%             'color', cell_map(i_cell, :), 'linewidth', 1)
    end
    ylabel('spike-triggered velocity')
    set(gca, 'color', 'none')
    axis square
    xlim(time_range_bins([1 end])*1000/60)
    %%
    print(gcf, '-dpng', ['scripts/Figure Plots/revisions/readoutsummary_withSTA_g' ...
        num2str(set_size) '_set' num2str(ex_info.set_ID)])
    close gcf
end


%%

figure()
hold on
plot(word_intPredI, readout_xc_predI, 'ko')
plot(word10_intPredI, readout10_xc_predI, 'bo')

%% 
param.data_source = 'fishmovieA';
param.movie_type = 'fish';
param.numL0 = 4;
param.training_reps = 1;
[~, train_s] = loadDataForLearningSimulation(668, param);


x_start = 4000;
x_length = 2*1141;

makeMyFigure(40, 8);
imagesc(1:size(train_s, 2), [], train_s)
colormap(flipud(gray))
xlim(x_start + [1 x_length])
set(gca, 'ytick', 1:4)
set(gca, 'xtick', x_start + [x_length/2 x_length], 'xticklabel', [1 2])
i_trial = 1;

test_mat = squeeze(test_s(:, :, i_trial));

% print(gcf, '-dpdf', 'scripts/Figure Plots/trainingdata_snapshot')



%% load raw spike times for this set of cells

 ex_info = infoX10{197};
 
 cell_set = ex_info.X_stim.cell_IDs;
 
 x_ST = load('2011-02-25_barfishcheck_experiment_and_calcs_PNAS_paper.mat', ...
     'SpikeTimes', 'cells_to_use', 'block_begin_time', 'block_end_time', ...
     'repbardata', 'repbarstim');
 
 useSpikeTimes = x_ST.SpikeTimes(x_ST.cells_to_use);
 
 cell_set_SpikeTimes = useSpikeTimes(cell_set);
 
 
%% load NRB spike times

cell_IDs = cell_set; %[13 16 27 49] ; % AKA set 668_fA4
set_size = length(cell_IDs);
comp_type = 'CDM_pv';
uniq_visited_percs = zeros(1, set_size);
[nrb_stim, nrb_data, learned_y_out] = barmovie_activity_readout_function(uniq_visited_percs, cell_IDs, comp_type);

i_seg = 10;

x_stim = nrb_stim{i_seg};
t_vec = (1:length(x_stim))/60;


makeMyFigure(20, 8);

subplot(3, 1, 1)
plot(t_vec, x_stim, 'k', 'linewidth', 1)
hold on
plot(max(t_vec) - [2 1], min(x_stim)*[1 1], 'r', 'linewidth', 1)
axis tight
set(gca, 'color', 'none', 'box', 'off', 'xtick', [], 'ytick', [])
subplot(3, 1, [2 3])
imagesc(t_vec, 1:set_size, nrb_data{i_seg}')
set(gca, 'ytick', [1 set_size], 'xtick', [])
colormap(flipud(gray))
% print(gcf, '-dpdf', 'scripts/Figure Plots/nrb_spike_ex')


%% repeating bar example

cell_IDs = [13 16 27 49] ; % AKA set 668_fA4
set_size = length(cell_IDs);


i_seg = 10;
x_stim = x_ST.repbarstim;

t_vec = (1:length(x_stim))/60;

spike_mat = x_ST.repbardata(cell_IDs, :, :);


makeMyFigure(15, 4);

% subplot(5, 1, 1)
plot(t_vec, x_stim, 'k', 'linewidth', 1)
hold on
plot(max(t_vec) - [2 1], min(x_stim)*[1 1], 'r', 'linewidth', 1)
axis tight
set(gca, 'color', 'none', 'box', 'off', 'xtick', [], 'ytick', [])

% print(gcf, '-dpdf', 'scripts/Figure Plots/repbar_stim_ex')


makeMyFigure(15, 15);

for i_sp = 1:length(cell_IDs)
    subplot(length(cell_IDs), 1, i_sp);
    imagesc(t_vec, [], squeeze(spike_mat(i_sp, :, :)))
    set(gca, 'ytick', [1 size(spike_mat, 2)], 'xtick', [])
    colormap(flipud(gray))
    ylabel(['RGC ' num2str(i_sp)])
end
% print(gcf, '-dpdf', 'scripts/Figure Plots/repbar_spike_ex')
 %%
 
Fs_sample = 10000;
time_period = 30*Fs_sample;
start_time = x_ST.block_begin_time{2}(1) + 20*Fs_sample;
% start_time = min(cellfun(@min, cell_set_SpikeTimes));
% stop_time = max(cellfun(@max, cell_set_SpikeTimes));


inRangeSpikes = cellfun(@(x) double(x( x > start_time & x < (start_time + time_period))), ...
    cell_set_SpikeTimes, 'uniformoutput', false);
inRangeSpikes = cellfun(@(x) reshape(x, 1, numel(x)), inRangeSpikes, 'uniformoutput', false);
makeMyFigure(20, 20);

zoom_x = [23.95 24.25];
spikeTimeBins = zoom_x(1):(1/60):zoom_x(2);

x_bin = spikeTimeBins(5) + [0 1/60];
cell_square = [1 1 11 11 1] - 0.5;
x_sq_inds = [1 2 2 1 1];
subplot(2, 1, 1)
hold on
for i_cell = 1:length(inRangeSpikes)
    t_values = inRangeSpikes{i_cell} - start_time;
    y_values = i_cell + [0*t_values; 1 + 0*t_values] - 0.5;
    plot([t_values;t_values]/Fs_sample, y_values, 'k')
end
set(gca, 'color','none')
set(gca, 'ytick', 1:length(inRangeSpikes))

plot([24 25], [0 0],'k', 'linewidth', 1.5)
plot(zoom_x(x_sq_inds), cell_square, 'k')
subplot(2, 3, 4)
hold on
for i_cell = 1:length(inRangeSpikes)
    t_values = inRangeSpikes{i_cell} - start_time;
    y_values = i_cell + [0*t_values; 1 + 0*t_values] - 0.5;
    plot([t_values;t_values]/Fs_sample, y_values, 'k')
end
set(gca, 'color','none')
set(gca, 'ytick', 1:length(inRangeSpikes))
xlim(zoom_x)

plot(zoom_x(end) - [0 6/60], [0 0],'k', 'linewidth', 1.5)

plot(x_bin(x_sq_inds), cell_square, 'k')
ylim([0 11])



binnedSpikeTimes = cell2mat(cellfun(@(x) histcounts((x - start_time)/Fs_sample, spikeTimeBins), ...
    inRangeSpikes, 'uniformoutput', false)');

subplot(2, 3, [5 6])
text(0, 0.5, num2str(flipud(binnedSpikeTimes~= 0)))
set(gca, 'visible', 'off')

% print(gcf, '-dpdf', 'scripts/Figure Plots/spikeTimes_to_Words_ex')

%% what can we learn about feature predictiveness by looking at the diversity of 
% cell subgroup predictive info values?

which_info = infoX10;

dWpI = cellfun(@(x) x.X_ww.dword_pI(1), which_info);

wN_stimIvT = (cellfun(@(x) x.X_stim.wordstim_predI_vt, ...
    which_info, 'uniformoutput', false));
stimPredI_T = cellfun(@(x) x.X_stim.tau_futures, which_info, 'uniformoutput', false);

[wN_maxStimI, w4MaxSPI_ind] = cellfun(@(x) max(x), wN_stimIvT);

wN_stimPredI = cellfun(@(x, y) x(y == 1), wN_stimIvT, stimPredI_T);


wN_intnrbPredI =  cellfun(@(x) x.X_stim.wordword_nrbpredI_vt(...
    x.X_stim.tau_futures == 1), which_info);
wN_nrbFR = cellfun(@(x) sum(x.X_stim.nrb_fr_cells), which_info);

wN_intPredI = cellfun(@(x) x.X_ww.word_pI(x.X_ww.tau_futures == 1), which_info);
wN_fishFR = cellfun(@(x) sum(x.X_ww.word_fr), which_info);

% % if plot in bits/spike, divide by firing rates
% wN_maxStimI = wN_maxStimI./wN_nrbFR;
% wN_stimPredI = wN_stimPredI./wN_nrbFR;
%%
num_bins = 4;
% ext_predI_binedges = [0 linspace(max(dWpI), max(wN_stimPredI), num_bins)];
% ext_info_binedges = [0 linspace(max(dWpI), max(wN_maxStimI), num_bins)];
ext_predI_binedges = prctile(wN_stimPredI, linspace(0, 100, num_bins + 1)); %   linspace(0, max(wN_stimPredI), num_bins+1);
ext_info_binedges = prctile(wN_maxStimI, linspace(0, 100, num_bins + 1)); %linspace(0, max(wN_maxStimI), num_bins+1);

fr_nrb_binedges = prctile(wN_nrbFR, linspace(0, 100, num_bins + 1));
int_nrb_binedges = [0 linspace(max(dWpI), max(wN_intnrbPredI), num_bins)];
int_fish_binedges = [0 linspace(max(dWpI), max(wN_intPredI), num_bins)];


stimI_cmap = jet(num_bins);
% stimI_cmap(1, :) = 0;

makeMyFigure(30, 30);
colormap(stimI_cmap)
subplot(221)
hold on
for i_pIb = 1:num_bins
    inds = wN_stimPredI > ext_predI_binedges(i_pIb) & ...
        wN_stimPredI <= ext_predI_binedges(i_pIb + 1);
%     plot(wN_intnrbPredI(inds), wN_intPredI(inds), '.', 'color', stimI_cmap(i_pIb, :), ...
%         'markersize', 15);
    
    plot(wN_intnrbPredI(inds)./wN_nrbFR(inds), wN_intPredI(inds)./wN_fishFR(inds), '.', 'color', stimI_cmap(i_pIb, :), ...
        'markersize', 15);
end
eqline
set(gca, 'fontsize', 12, 'color', 'none')
axis square
ch = colorbar;
ch.Ticks = linspace(0, 1, num_bins + 1);
ch.TickLabels = num2str(ext_predI_binedges', '%1.1f');
ch.Label.String = 'ext. pred info (bits/spk)';
xlabel('internal pred-I (bar)')
ylabel('internal pred-I (movie)')

subplot(222)
hold on
for i_pIb = 1:num_bins
    inds = wN_nrbFR > fr_nrb_binedges(i_pIb) & ...
        wN_nrbFR <= fr_nrb_binedges(i_pIb + 1);
    plot(wN_intnrbPredI(inds)./wN_nrbFR(inds), wN_stimPredI(inds), '.', 'color', stimI_cmap(i_pIb, :), ...
        'markersize', 15);
end
eqline
set(gca, 'fontsize', 12, 'color', 'none')
axis square
ch = colorbar;
ch.Ticks = linspace(0, 1, num_bins + 1);
ch.TickLabels = num2str(60*fr_nrb_binedges', '%1.1f');
ch.Label.String = 'nrb fr';
xlabel('internal pred-I (bar)')
ylabel('external pred-I (bar)')


subplot(223)
hold on
for i_pIb = 1:num_bins
    inds = wN_intnrbPredI > int_nrb_binedges(i_pIb) & ...
        wN_intnrbPredI <= int_nrb_binedges(i_pIb + 1);
    plot(wN_stimPredI(inds), wN_maxStimI(inds), '.', 'color', stimI_cmap(i_pIb, :), ...
        'markersize', 15);

end
eqline
set(gca, 'fontsize', 12, 'color', 'none')
axis square
ch = colorbar;
ch.Ticks = linspace(0, 1, num_bins + 1);
ch.TickLabels = num2str(60*int_nrb_binedges', '%1.1f');
ch.Label.String = 'int. (fish) pred info';
xlabel('external pred-I (bar)')
ylabel('external stim info (bar)')

subplot(224)
hold on
for i_pIb = 1:num_bins
    inds = wN_intPredI > int_fish_binedges(i_pIb) & ...
        wN_intPredI <= int_fish_binedges(i_pIb + 1);
    plot(wN_intnrbPredI(inds)./wN_nrbFR(inds), wN_stimPredI(inds), '.', 'color', stimI_cmap(i_pIb, :), ...
        'markersize', 15);
end
eqline
set(gca, 'fontsize', 12, 'color', 'none')
axis square
ch = colorbar;
ch.Ticks = linspace(0, 1, num_bins + 1);
ch.TickLabels = num2str(60*int_fish_binedges', '%1.1f');
ch.Label.String = 'int. (fish) info';
xlabel('internal pred-I (bar)')
ylabel('external pred-I (bar)')
