function [p_v_mat_arr, v_sta_arr, velocity_bins, time_range_bins] = computeVelocitySTA_NRBdata(nrb_stim, nrb_data, time_range_bins)
% computes the spike-triggered velocity distribution given raw stimulus
% data (nrb_stim) and single-bit readouts in cell array nrb_data. nrb_data
% could either be full of nTx(# perceptron readouts) OR nTx(# cells) binary
% spike mats. For either, the STA/spike-triggered distributions are also
% computed for columns of nrb_data{ii} individually. 

num_readouts = size(nrb_data{1}, 2);

p_v_mat_arr = cell(num_readouts, 1);
v_sta_arr = cell(num_readouts, 2);

% find reasonable min/max velocities from nrb_stim. Taking median discards
% weird outliers. 
min_velocity = floor(median(cellfun(@(x) min(diff(x)), nrb_stim))-3);
max_velocity = ceil(median(cellfun(@(x) max(diff(x)), nrb_stim))+3);

velocity_bins = min_velocity:1:max_velocity;

% have to break this down on a trial-by-trial basis, as data were not
% acquired in a single stretch
num_trials = length(nrb_stim);

for i_perc = 1:num_readouts

    trial_hist_arr = zeros(num_trials, length(velocity_bins), length(time_range_bins));
    for i_trial = 1:num_trials
        % have to compute velocity by taking diff
        trial_velocities = diff(nrb_stim{i_trial});
        % drop the first bin to compensate for the velocity calculation. By
        % dropping the first bin, rather than the last bin, the velocity
        % estimate is shifted back relative to the spiking activity. 
        trial_spikes = find(nrb_data{i_trial}(2:end, i_perc));

        % discard any spikes that occur too close to start or end of trial
        tr_sp_min = trial_spikes + time_range_bins(1);
        tr_sp_max = trial_spikes + time_range_bins(end);
        keep_spikes = tr_sp_min >= 1 & tr_sp_max <= length(trial_velocities);
        trial_spikes = trial_spikes(keep_spikes);
    %     display([num2str(mean(keep_spikes)*100, '%1.0f') '% of ' ...
    %         num2str(length(keep_spikes)) ' output spikes fell within range'])

        if ~isempty(trial_spikes)
            spike_vel_snippets = trial_velocities(...
                repmat(trial_spikes, 1, length(time_range_bins)) + ...
                repmat(time_range_bins, length(trial_spikes), 1));

        %     trial_v_hist = zeros(length(velocity_bins), length(time_range_bins));
            for i_trb = 1:length(time_range_bins)
                trial_hist_arr(i_trial, :, i_trb) = histcounts(spike_vel_snippets(:, i_trb), ...
                    'binlimits', velocity_bins([1 end]), 'binmethod', 'integers');

            end
        end
    end
    % get counts over all trials
    ct_v_mat = squeeze(sum(trial_hist_arr, 1));
    % normalize by counts to get probabilities 
    p_v_mat = ct_v_mat*diag(1./sum(ct_v_mat));

    p_v_mat_arr{i_perc} = p_v_mat;

    % take average velocity at each time bin to get the STA for velocity
    v_sta_arr{i_perc} = velocity_bins*p_v_mat;
    vSqVal = (velocity_bins.^2)*p_v_mat;
    v_sta_arr{i_perc, 2} = sqrt(vSqVal - v_sta_arr{i_perc, 1}.^2);
end