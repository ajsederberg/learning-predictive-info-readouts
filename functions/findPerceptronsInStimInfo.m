function [stim_info_inds, field_vals] = findPerceptronsInStimInfo(X_infostruct, perc_w_inds, perc_dict, field_name, tau_future_vals)
% Helper function to look up the indices of the perceptrons in X_stim that
% correspond to the weights in perc_ws. 

if ~isempty(perc_w_inds)
    % First use the dictionary to look up the weights of the perceptrons listed
    % in perc_w_inds
    set_size = size(perc_dict.percX.unique_perceptrons, 2);
    if set_size ==4
        perc_rules = perc_dict.percX.unique_partitions(perc_w_inds, :);
    else
        perc_rules = perc_dict.percX.final_partitions(perc_w_inds, :);
    end

    % now match those to the perceptrons for which stim info was calculated
    stim_rules = perceptronToProjectionRule(X_infostruct.X_stim.perc_ws, 1);
    % this finds a list of indices (length = # of inds in perc_w_inds) for
    % which of the stim-info claculations corresponds to the dictionary
    % perceptron. 
    [~, stim_info_inds] = checkPartitionUniqueness_RPMethod(stim_rules, perc_rules);
else
    % this is a way to return the word-word info at the specified values of
    % tau_future_vals
    stim_info_inds = 1;
end
if nargin > 3
% return values in the field named by field_name. Options for field names
% are: percstim_predI_vt, percword_nrbpredI_vt, nrb_fr_perc, perc_ws
    if ~isempty(tau_future_vals)
        [~, output_inds, input_inds] = intersect(tau_future_vals, X_infostruct.X_stim.tau_futures);
        field_vals = nan(length(stim_info_inds), length(tau_future_vals));
    else
        field_vals = nan(length(stim_info_inds), 1);
        input_inds = 1;
        output_inds = 1;
    end
    if ~any(stim_info_inds == 0)
        field_vals(:, output_inds) = X_infostruct.X_stim.(field_name)(stim_info_inds, input_inds);
    else
        field_vals(stim_info_inds ~= 0, output_inds) =  ...
            X_infostruct.X_stim.(field_name)(stim_info_inds(stim_info_inds~= 0), input_inds);
    end
end