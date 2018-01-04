function [outX, file_name] = computeBarmovieStimulusPredI(perc_ws, tau_futures, set_num, comp_type, data_source)
% Returns the stimulus predictive information calculations for all
% perceptrons in the matrix learned_ws in a struct with field: 
% 
% outX.percstim_predI_vt = readout-stimulus predictive info;
% outX.percword_nrbpredI_vt = readout-word internal predictive info on NRB set;
% outX.wordstim_predI_vt = word-stimulus predictive info;
% outX.cellstim_predI_vt = individual cell-stimulus predictive info; 
% outX.cell_IDs = cell_IDs;
% outX.set_num = set_num;
% outX.tau_futures = tau_futures; 
% outX.perc_ws = perc_ws;
% outX.nrb_fr_perc = nrb_fr_perc; 
% outX.nrb_fr_cells = nrb_fr_cells;
% outX.w_labels = w_labels;

% get the IDs of cells is set "set_num"
numL0 = size(perc_ws, 2);
X_groups = load(['Cell Groups/' data_source '_setsOf' num2str(numL0) '.mat']);
cell_IDs = X_groups.sets_of_N(set_num, :);

X = load('2011-02-25_barfishcheck_experiment_and_calcs_PNAS_paper.mat', ...
    'nonrepbardata', 'nonrepbarstim');
nonrepbardata = X.nonrepbardata;
nonrepbarstim = X.nonrepbarstim;

% find appropriate boundaries for equipartition bins from nonrepbarstim
% there will be many 0 values that are not stimuli
nzbnd_nrbstim = nonzeros(nonrepbarstim);
nzbnd_nrbstim(nzbnd_nrbstim > 260) = nan;
[cdf_stim, cdf_edges] = histcounts(nzbnd_nrbstim, 'binmethod', 'integer', 'normalization', 'cdf');

int_stim_bins = edges2bins(cdf_edges);

use_bins = cdf_stim < 0.99 & cdf_stim > 0.01;
cdf_stim = cdf_stim(use_bins);
int_stim_bins = int_stim_bins(use_bins);
num_bins = 37;
eq_bin_edges = interp1([0 cdf_stim 1], [0 int_stim_bins 260.5], linspace(0, 1, num_bins+1));
eq_hist = histcounts(nonrepbarstim, eq_bin_edges);

eq_bins = (edges2bins(eq_bin_edges));
% check! bars should be about the same height
figure()
bar(eq_bins(2:end)', eq_hist(2:end)')

if strcmp(comp_type, 'CDM_pv') || strcmp(comp_type, 'SEP_pv')
    %%
    nan_nrbstim = nonrepbarstim;
    nan_nrbstim(nonrepbarstim>260) = nan;
    nan_nrbstim(nonrepbarstim == 0) = nan;
    
    v_nrbstim = diff(nan_nrbstim, 1, 2);
    
    [cdf_stim, cdf_edges] = histcounts(v_nrbstim, 'binmethod', 'integer', 'normalization', 'cdf');

    int_stim_bins = edges2bins(cdf_edges);

    use_bins = cdf_stim < 0.99 & cdf_stim > 0.01;
    cdf_stim = cdf_stim(use_bins);
    int_stim_bins = int_stim_bins(use_bins);
    num_bins = 37;
    eq_vbin_edges = interp1([0 cdf_stim 1], [int_stim_bins(1)-1 int_stim_bins int_stim_bins(end)+1], linspace(0, 1, num_bins+1));
    eq_vbin_edges = unique(floor(eq_vbin_edges) + 0.5);
%     eq_hist = histcounts(v_nrbstim, eq_vbin_edges);

end
%%


nrb_data = cell(size(nonrepbardata, 2), 1);
nrb_stim = cell(size(nonrepbardata, 2), 1); 
if strcmp(comp_type, 'CDM_pv') || strcmp(comp_type, 'SEP_pv')
    nrb_binned_stim = cell(size(nonrepbardata, 2), 2); 
else
    nrb_binned_stim = cell(size(nonrepbardata, 2), 1); 
end
learned_y_out = cell(size(nonrepbardata, 2), 1);

if all(isnan(perc_ws(:)))   % compute word-number representation for the token-word information
    full_word_record = reshape(permute(nonrepbardata(cell_IDs, :, :), [3 2 1]), ...
        size(nonrepbardata, 2)*size(nonrepbardata,3), length(cell_IDs)  );
    
    word_record_numbers = full_word_record*(2.^(0:1:(numL0-1))');
    [word_labels, wl_inds] = unique(word_record_numbers);

    word_tokens = full_word_record(wl_inds, :);
    keep_words = sum(word_tokens, 2) > 0 & sum(word_tokens, 2) < 4;
    word_labels = word_labels(keep_words);
    word_tokens = word_tokens(keep_words, :);
    
    outX.word_tokens = word_tokens;
    
    perc_ws = nan(length(word_labels), numL0);
end

for i_trial = 1:size(nonrepbardata, 2)
    t_vec = nonrepbarstim(i_trial, :) ~= 0;
    nrb_data{i_trial} = squeeze(nonrepbardata(cell_IDs, i_trial, t_vec))';
    nrb_stim{i_trial} = nonrepbarstim(i_trial, t_vec);
    
    % get the bins corresponding the the stim
    [~, ~, stim_bin_inds] = histcounts(nrb_stim{i_trial}, eq_bin_edges);
    nrb_binned_stim{i_trial, 1} = stim_bin_inds';
    if strcmp(comp_type, 'CDM_pv') || strcmp(comp_type, 'SEP_pv')
        [~, ~, stim_vbin_inds] = histcounts(diff(nrb_stim{i_trial}), eq_vbin_edges);        
        nrb_binned_stim{i_trial, 2} = stim_vbin_inds';
        % drop last value so the vectors are the same length
        nrb_binned_stim{i_trial, 1} = nrb_binned_stim{i_trial, 1}(1:end-1);
        nrb_data{i_trial} = nrb_data{i_trial}(1:end-1, :);
    end

    
    cell_data = squeeze(nonrepbardata(cell_IDs, i_trial, t_vec));
    
    if all(isnan(perc_ws(:)))
        % code for computing token-word information
        trial_data = cell_data'*(2.^(0:1:(numL0-1))');
        
        % for each token in the full dataset, construct a readout
        y_out = zeros(size(trial_data, 1), size(word_tokens,1));
        for i_wt = 1:size(word_tokens,1)
            y_out(trial_data == word_labels(i_wt), i_wt) = 1;
        end
        learned_y_out{i_trial}  = y_out;
    else
        learned_y_out{i_trial} = double(perc_ws*cell_data > 1)';
    end
    if strcmp(comp_type, 'CDM_pv') || strcmp(comp_type, 'SEP_pv')
        learned_y_out{i_trial} = learned_y_out{i_trial}(1:end-1, :);
    end
end

%% now compute some curves
percstim_predI_vt = zeros(size(perc_ws, 1), length(tau_futures));
percword_nrbpredI_vt = zeros(size(perc_ws, 1), length(tau_futures));
wordword_nrbpredI_vt = zeros(1, length(tau_futures));
wordstim_predI_vt = zeros(1, length(tau_futures));
cellstim_predI_vt = zeros(size(perc_ws, 2), length(tau_futures));
dpercstim_predI_vt = zeros(size(perc_ws, 1), length(tau_futures));
dpercword_nrbpredI_vt = zeros(size(perc_ws, 1), length(tau_futures));
dwordword_nrbpredI_vt = zeros(1, length(tau_futures));
dwordstim_predI_vt = zeros(1, length(tau_futures));
dcellstim_predI_vt = zeros(size(perc_ws, 2), length(tau_futures));

% set this higher than the lenght of tau_futures if you want all values
% caluclated; set lower if you want a limited number calculated. 
max_num_tau_futures = 25;
if length(tau_futures) > max_num_tau_futures
    % compute for tau_futures, but compute for tau_futures == 1 first. If that
    % is < 0.02
    outof_order_tau_futures = [1 setdiff(tau_futures, 1)];

    for ii = 1:length(tau_futures)
        tau_future = outof_order_tau_futures(ii);

        if ii == 1 
            % info between readout and stimulus

            Ipred_ystim = computePredictiveInfo_FSS_forNRB(learned_y_out, nrb_binned_stim, tau_future, comp_type);

            % info between readout and word
            Ipred_ynrb = computePredictiveInfo_FSS_forNRB(learned_y_out, nrb_data, tau_future, 'CDM');
            
            % internal info - word-word
            if tau_future > 0
                Ipred_wordnrb = computePredictiveInfo_FSS_forNRB(nrb_data, nrb_data, tau_future, 'CDM');
            else
                Ipred_wordnrb = [0;0];
            end
            % info between word/cells and stimulus
            Ipred_wordstim = computePredictiveInfo_FSS_forNRB(nrb_data, nrb_binned_stim, tau_future, comp_type);

            
            percstim_predI_vt(:, ii) = Ipred_ystim(1:end-1, 1);
            percword_nrbpredI_vt(:, ii) = Ipred_ynrb(1:end-1, 1);
            wordstim_predI_vt(ii) = Ipred_wordstim(end, 1);
            wordword_nrbpredI_vt(ii) = Ipred_wordnrb(end, 1);
            cellstim_predI_vt(:, ii) = Ipred_wordstim(1:end-1, 1);

            dpercstim_predI_vt(:, ii) = Ipred_ystim(1:end-1, 2);
            dpercword_nrbpredI_vt(:, ii) = Ipred_ynrb(1:end-1, 2);
            dwordstim_predI_vt(ii) = Ipred_wordstim(end, 2);
            dwordword_nrbpredI_vt(ii) = Ipred_wordnrb(end, 2);
            dcellstim_predI_vt(:, ii) = Ipred_wordstim(1:end-1, 2);
        else
                % for second and later, do not compute all tau_future for the
                % perceptrons with very low stimulus information. 
            hi_pI_inds = percstim_predI_vt(:, 1) > 0.015;

            part_ly_out = cellfun(@(x) x(:, hi_pI_inds), learned_y_out, 'uniformoutput', false);

            Ipred_ystim = computePredictiveInfo_FSS_forNRB(part_ly_out, nrb_binned_stim, tau_future, comp_type);

            % info between readout and word
            Ipred_ynrb = computePredictiveInfo_FSS_forNRB(part_ly_out, nrb_data, tau_future, 'CDM');
            % info between word/cells and stimulus
            Ipred_wordstim = computePredictiveInfo_FSS_forNRB(nrb_data, nrb_binned_stim, tau_future, comp_type);
            % internal info - word-word
            if tau_future > 0
                Ipred_wordnrb = computePredictiveInfo_FSS_forNRB(nrb_data, nrb_data, tau_future, 'CDM');
            else
                Ipred_wordnrb = [0;0];
            end

            percstim_predI_vt(hi_pI_inds, ii) = Ipred_ystim(1:end-1, 1);
            percword_nrbpredI_vt(hi_pI_inds, ii) = Ipred_ynrb(1:end-1, 1);
            wordword_nrbpredI_vt(ii) = Ipred_wordnrb(end, 1);
            wordstim_predI_vt(ii) = Ipred_wordstim(end, 1);
            cellstim_predI_vt(:, ii) = Ipred_wordstim(1:end-1, 1);

            dpercstim_predI_vt(hi_pI_inds, ii) = Ipred_ystim(1:end-1, 2);
            dpercword_nrbpredI_vt(hi_pI_inds, ii) = Ipred_ynrb(1:end-1, 2);
            dwordword_nrbpredI_vt(ii) = Ipred_wordnrb(end, 2);
            dwordstim_predI_vt(ii) = Ipred_wordstim(end, 2);
            dcellstim_predI_vt(:, ii) = Ipred_wordstim(1:end-1, 2);
        end 
        
        display(['stim progress: ' num2str(ii) ' of ' num2str(length(tau_futures))])
    end

    [~, ord] = sort(outof_order_tau_futures);
    percstim_predI_vt = percstim_predI_vt(:, ord);
    percword_nrbpredI_vt = percword_nrbpredI_vt(:, ord);
    wordstim_predI_vt = wordstim_predI_vt(ord);
    wordword_nrbpredI_vt = wordword_nrbpredI_vt(ord);
    cellstim_predI_vt = cellstim_predI_vt(:, ord);
    
    dpercstim_predI_vt = dpercstim_predI_vt(:, ord);
    dpercword_nrbpredI_vt = dpercword_nrbpredI_vt(:, ord);
    dwordstim_predI_vt = dwordstim_predI_vt(ord);
    dwordword_nrbpredI_vt = dwordword_nrbpredI_vt(ord);
    dcellstim_predI_vt = dcellstim_predI_vt(:, ord);
    
    
else
    
    
    
    % if there are only max_num_tau_futures tau_futures or fewer, compute for all 

    for ii = 1:length(tau_futures)
        tau_future = tau_futures(ii);


        Ipred_ystim = computePredictiveInfo_FSS_forNRB(learned_y_out, nrb_binned_stim, tau_future, comp_type);

        % info between readout and word
        Ipred_ynrb = computePredictiveInfo_FSS_forNRB(learned_y_out, nrb_data, tau_future, comp_type(1:3));
        % info between word/cells and stimulus
        Ipred_wordstim = computePredictiveInfo_FSS_forNRB(nrb_data, nrb_binned_stim, tau_future, comp_type);

        % internal info - word-word
        if tau_future > 0
            Ipred_wordnrb = computePredictiveInfo_FSS_forNRB(nrb_data, nrb_data, tau_future, comp_type(1:3));
        else
            Ipred_wordnrb = [0 0;0 0];
        end
        
        
        
        percstim_predI_vt(:, ii) = Ipred_ystim(1:end-1, 1);
        percword_nrbpredI_vt(:, ii) = Ipred_ynrb(1:end-1, 1);
        wordword_nrbpredI_vt(ii) = Ipred_wordnrb(end, 1);
        wordstim_predI_vt(ii) = Ipred_wordstim(end, 1);
        cellstim_predI_vt(:, ii) = Ipred_wordstim(1:end-1, 1);

        dpercstim_predI_vt(:, ii) = Ipred_ystim(1:end-1, 2);
        dpercword_nrbpredI_vt(:, ii) = Ipred_ynrb(1:end-1, 2);
        dwordword_nrbpredI_vt(ii) = Ipred_wordnrb(end, 2);
        dwordstim_predI_vt(ii) = Ipred_wordstim(end, 2);
        dcellstim_predI_vt(:, ii) = Ipred_wordstim(1:end-1, 2);
  
    end


    
    


end
%%
bin_size = 1000/60;
nrb_fr_perc = mean(cell2mat(learned_y_out))';
nrb_fr_cells = mean(cell2mat(nrb_data))';
%%
if nargout == 2
    makeMyFigure(20,20);

    hold on
    plot(bin_size*tau_futures, diag(1./nrb_fr_perc)*percstim_predI_vt, 'linewidth', 2);
    % plot(tau_futures, diag(1./nrb_fr_cells)*cellstim_predI_vt, 'b')
    plot(bin_size*tau_futures, wordstim_predI_vt/sum(nrb_fr_cells), 'k', 'linewidth', 2);
    hold off
    set(gca, 'fontsize', 14)
    xlabel('dt (ms)')
    ylabel('I( position(t+dt), perceptron(t)) (bits per spike)')

    %%


    file_name = ['stimulus info results/' data_source '_' ...
        '_set' num2str(set_num) '_g' num2str(length(cell_IDs)) ...
        '_perceptrons_stimulus_info_' comp_type];
    print(gcf, '-dpdf', file_name);

    close all

    saveInParfor(file_name, ...
        percstim_predI_vt, percword_nrbpredI_vt, wordstim_predI_vt, cellstim_predI_vt, wordword_nrbpredI_vt, ...
        dpercstim_predI_vt, dpercword_nrbpredI_vt, dwordstim_predI_vt, dcellstim_predI_vt, dwordword_nrbpredI_vt, ...
        cell_IDs, set_num, tau_futures, perc_ws, nrb_fr_perc, nrb_fr_cells);
end
outX.file_tag = ['stiminfo_set' num2str(set_num) '_' data_source([1 end]) num2str(numL0)];
outX.percstim_predI_vt = percstim_predI_vt;
outX.dpercstim_predI_vt = dpercstim_predI_vt;

outX.percword_nrbpredI_vt = percword_nrbpredI_vt;
outX.dpercword_nrbpredI_vt = dpercword_nrbpredI_vt;

outX.wordstim_predI_vt = wordstim_predI_vt;
outX.dwordstim_predI_vt = dwordstim_predI_vt;

outX.wordword_nrbpredI_vt = wordword_nrbpredI_vt;
outX.dwordword_nrbpredI_vt = dwordword_nrbpredI_vt;

outX.cellstim_predI_vt = cellstim_predI_vt; 
outX.dcellstim_predI_vt = dcellstim_predI_vt; 

outX.cell_IDs = cell_IDs;
outX.set_num = set_num;
outX.tau_futures = tau_futures; 
outX.perc_ws = perc_ws;
outX.nrb_fr_perc = nrb_fr_perc; 
outX.nrb_fr_cells = nrb_fr_cells;
