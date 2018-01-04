function [nrb_stim, nrb_data, learned_y_out] = barmovie_activity_readout_function(learned_ws, cell_IDs, comp_type)
% this returns a cell array with stimulus information (bar position) for
% each trial, the cell set rasters for each trial, and the perceptron
% readout for all perceptrons in the 'learned_ws' matrix. Note that
% 'cell_IDs' and the columns of 'learned_ws' MUST correspond. 

X = load('2011-02-25_barfishcheck_experiment_and_calcs_PNAS_paper.mat', 'nonrepbardata', ...
    'nonrepbarstim');
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


if strcmp(comp_type, 'CDM_pv')
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
if strcmp(comp_type, 'CDM_pv')
    nrb_binned_stim = cell(size(nonrepbardata, 2), 2); 
else
    nrb_binned_stim = cell(size(nonrepbardata, 2), 1); 
end
learned_y_out = cell(size(nonrepbardata, 2), 1);

for i_trial = 1:size(nonrepbardata, 2)
    t_vec = nonrepbarstim(i_trial, :) ~= 0;
    nrb_data{i_trial} = squeeze(nonrepbardata(cell_IDs, i_trial, t_vec))';
    nrb_stim{i_trial} = nonrepbarstim(i_trial, t_vec);
    
    % get the bins corresponding the the stim
    [~, ~, stim_bin_inds] = histcounts(nrb_stim{i_trial}, eq_bin_edges);
    nrb_binned_stim{i_trial, 1} = stim_bin_inds';
    if strcmp(comp_type, 'CDM_pv')
        [~, ~, stim_vbin_inds] = histcounts(diff(nrb_stim{i_trial}), eq_vbin_edges);        
        nrb_binned_stim{i_trial, 2} = stim_vbin_inds';
        % drop last value so the vectors are the same length
        nrb_binned_stim{i_trial, 1} = nrb_binned_stim{i_trial, 1}(1:end-1);
        nrb_data{i_trial} = nrb_data{i_trial}(1:end-1, :);
    end

    
    cell_data = squeeze(nonrepbardata(cell_IDs, i_trial, t_vec));
    learned_y_out{i_trial} = double(learned_ws*cell_data > 1)';
    if strcmp(comp_type, 'CDM_pv')
        learned_y_out{i_trial} = learned_y_out{i_trial}(1:end-1, :);
    end
end

