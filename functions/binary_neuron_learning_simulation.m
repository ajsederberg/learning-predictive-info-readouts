function [trialW, trialFR] = binary_neuron_learning_simulation(i_set, param)
% FUNCTIONALIZED from super_simple_binary_neuron_script_V2.m
% data_source: 'model' or 'retina'
% i_set: index of set (range depends on data source)
% learning type: 'vanilla', 'bcmlike', 'homeostatic'
% compute_info: 0 or 1; whether to calculate pred-I after each traning
% chunk. 
% param: all model parameters. Should have fields:
%           - numL1: integer, number of initial conditions tested
%           - w0max: upper bound for initial weights 
%           - wmaxes : (loop variable, vanilla or homeostatic) maximum
%           weights
%           - alphas_ltd: (loop variable) LTD strengths
%           - target_frs: (loop variable, bcmlike only) target firing rates
%           - NT: number of training repetitions to lump into each "trial".
%                   for retina data: needs to be a divisor of 102
%                   for model data: probably needs to be a multiple of 10

% load parameters
binsize = 1000/60; % bin size is always 1000/60
numL1 = param.numL1;
learning_type = param.learning_type;
data_source = param.data_source;

% %% 1: Load training data 

[train_s, test_s, i_set_string] = ...
    loadDataForLearningSimulation(i_set, param);

numL0 = size(test_s, 1);
num_training_trials = size(train_s, 3);

% save for speed analyses
total_training_spikes = sum(sum(train_s, 3), 2);
total_training_bins = numel(train_s)/numL0;

%%
tic
 
% set the three for-loop variables
if strcmp(learning_type, 'vanilla') ...
        || strcmp(learning_type, 'triplet') ...
        || strcmp(learning_type, 'triploh') ...
        || strcmp(learning_type, 'vanilloh')
    ivar1name = 'wmax';
    ivar1 = param.wmaxes;
    c1 = 100;
    
    ivar2name = 'ltd';
    ivar2 = param.alphas_ltd;
    c2 = 10;
    
    if strcmp(learning_type, 'triplet') || strcmp(learning_type, 'triploh')
        ivar3name = 'tauy';
        ivar3 = param.tau_y;
        c3 = 1;
    else
        ivar3name = 'noise';
        ivar3 = param.noise_sigmas;
        c3 = 100;
    end
    
end

 
% index for looping through different target firing rates or alphas_ltd,
% etc. 
for i_v3 = 1:length(ivar3)
    if strcmp(learning_type, 'triplet') || strcmp(learning_type, 'triploh')
        noise_sig_val = 0;
        % tau_y is the time constant for the second spike interaction in
        % the triplet rule. it is in units of time bins (i.e. tau_1 = 1 is
        % 16.7 ms)
        tau_y = ivar3(i_v3);
    else
        % this has been changed to 0*ivar3 to eliminate the possibility of
        % adding noise. 
        noise_sig_val = 0*ivar3(i_v3);
    end
for i_v1 = 1:length(ivar1)
    % assign named variable for better readability later in the script
    wmaxes = ivar1;

for i_v2 = 1:length(ivar2)
    alphas_ltd = ivar2;
    
    file_string = [ivar1name num2str(ivar1(i_v1)*c1) '_' ...
        ivar2name num2str(ivar2(i_v2)*c2) '_' ...
        ivar3name num2str(ivar3(i_v3)*c3) '_set' i_set_string]; 
    
    if isfield(param, 'return_file_string_only') && param.return_file_string_only
        trialW = file_string;
        trialFR = [];
        return;
    end
    %%

    b = 1;
  
    rng(numL0);   
    if isfield(param, 'w0_uniform') && param.w0_uniform 
        w = (w0max-w0min)*(rand(numL1, numL0)) + w0min;    
        w(w < w0min) = w0min;    % initial condition is positive weights
        w(w > w0max) = w0max;   % initial condition is below w0max
    elseif isfield(param, 'w0_values')
        % option to explicitly set w0_values
        w = param.w0_values;
        
    else
        w = 0.7*b + 0.25*b*randn(numL1, numL0);
        w(w < 0) = 0;
        w(w > wmaxes(i_v1)) = wmaxes(i_v1);
        % any rows that sum to less than 1.5, add a boost
        w_too_small = sum(w, 2) < 1.5;
        w(w_too_small, :) = w(w_too_small, :) + 1.5/numL0;
    end
        
    
    % Save the weights at each snapshot. Also save firing rate (Hz)
    trialW = zeros(num_training_trials + 1, numL1, numL0);
    trialFR = zeros(num_training_trials + 1, numL1);

    for i_trial = 1:num_training_trials
        pause(0.01)
        % store w at the start of the trial
        trialW(i_trial, :, :) = w;
        % start of each trial, compute firing rate over test segment
        syn_noise = randn(numL1, size(test_s, 2))*noise_sig_val;
        save_x = w*test_s + syn_noise > b;
        trialFR(i_trial, :) = mean(save_x, 2);

        % load dynamics parameters, minimum w strength and rate of learning
        wmin = param.wmin;
        alphaLTP = param.alphaLTP;

        %% now implement learning
        if i_trial > 20 
%             w_diff = abs(diff(trialW((i_trial-1):i_trial, :, :), 1, 1));
%             w_diff2 = abs(diff(trialW((i_trial-2):i_trial, :, :), 2, 1));
%             regl_diff = squeeze(w_diff(end, :, :))./(1+squeeze(trialW(i_trial, :, :)));
%             in_steadystate = all(regl_diff < 0.001, 2);

%             % saturation check
%             atMIN = abs(wmin - squeeze(...
%                 mean(trialW((i_trial-2):i_trial, :, :), 1))) < 0.001; 
%             atWMAX = abs(wmaxes(i_v1) - squeeze(...
%                 mean(trialW((i_trial-2):i_trial, :, :), 1))) < 0.001; 
%             in_steadystate = all( atMIN | atWMAX, 2);

            % check if this trial is equal to previous 10 trial average
            mean_last10Trials = squeeze(mean(trialW(i_trial + (-10:1:-1), :, :), 1));
            last_trial = squeeze(trialW(i_trial, :, :));
            diff_from_mean = abs(mean_last10Trials - last_trial);
            std_last10Trials = squeeze(std(trialW(i_trial + (-10:1:0), :, :), [], 1));
            in_steadystate = all(diff_from_mean < 0.005 & std_last10Trials < 0.005, 2);

        else 
            in_steadystate = false(numL1, 1);
        end
%         if sum(in_steadystate) > 0
%             dummy_var = 0;
%         end
        
        % this bipasses the initial conditions that have reached steady
        % state. 
        int_w = w(~in_steadystate, :);       
        x = zeros(sum(~in_steadystate), 1);    % start off

        
    if sum(~in_steadystate) > 0

        if strcmp(learning_type, 'triplet') || strcmp(learning_type, 'triploh')
                % This modulates LTP by exp(-timeSinceLastSPike/~100 ms). If
                % set to Inf initially, no LTP will ever occur. Set to nan and
                % only set to a real number after the first spike occurs.            
                    
            if sum(in_steadystate) < length(in_steadystate)
                % If some initial conditions have entered steady_state, drop
                % these from timeSinceLastSpike
                timeSinceLastSpike = nan(numL1, 1);
                timeSinceLastSpike = timeSinceLastSpike(~in_steadystate);
            else

                timeSinceLastSpike = nan(numL1, 1);
            end
        end

        for t=1:size(train_s, 2)
            inp = train_s(:, t, i_trial);
                
            % for learning type "vanilla" or "constantLTD"
            if strcmp(learning_type, 'vanilla') || strcmp(learning_type, 'vanilloh') ||...
                    strcmp(learning_type, 'antiSTDP') || strcmp(learning_type, 'constantLTD')
                [x, int_w] = binaryNeuronPair(x, int_w, inp, noise_sig_val, b, learning_type, ...
                    alphaLTP, alphas_ltd(i_v2), wmaxes(i_v1), wmin);

            elseif strcmp(learning_type, 'triplet') || strcmp(learning_type, 'triploh')
                % if current value of x is 1, timeSinceLastSpike is 1
                timeSinceLastSpike(x == 1) = 1;
                % otherwise add 1 to timeSinceLastSpike
                timeSinceLastSpike(x == 0) = timeSinceLastSpike(x == 0) + 1;
                % key difference between this and vanilla is the last
                % variable, 'timeSinceLastSpike'
                [x, int_w] = binaryNeuronPair(x, int_w, inp, noise_sig_val, b, learning_type, ...
                    alphaLTP, alphas_ltd(i_v2), wmaxes(i_v1), wmin, timeSinceLastSpike, tau_y);    

            end
            
        end
        % update w, ave_fr
        w(~in_steadystate, :) = int_w;   

    end


    end
    % save weights - end of simulation
    trialW(num_training_trials + 1, :, :) = w;
    
    
    % compute firing rate
    syn_noise = randn(numL1, size(test_s, 2))*noise_sig_val;
    test_proj = w*test_s + syn_noise > b;
    trialFR(num_training_trials + 1, :) = mean(test_proj, 2);


    
    %%
    h = makeMyFigure(30, 30); 
    cmap = lines(numL0);
    subplot(3, 1, 1)
    cla;
    hold on
    for i_l0 = 1:numL0
        plot(squeeze(trialW(:, :, i_l0)), 'color', cmap(i_l0, :));
    end
    hold off
    title('w during learning')
    ylabel('weight 1')
    xlabel('trial')
    
    subplot(3, 2, 5)
    max_fr = mean(sum(test_s, 1));
    histogram(trialFR(end, :)/max_fr)
    xlabel('firing rate - end sim')
    
    subplot(3, 2, 6)
    [~, w_ord] = sort(trialFR(end, :));
    imagesc(w(w_ord, :)', [0 wmaxes(i_v1)])
    colorbar
    title('w strength')

    % plot readout similarity to final rule 
    end_W = repmat((trialW(end, :, :)), [size(trialW, 1) 1 1]);
    w_difference_sim = mean(abs(trialW - end_W), 3);
    subplot(3, 1, 2)
    plot(w_difference_sim)
    title('L1 weight difference from final weights')
    
    %% save and print results snapshot
    if param.saveandprint == 1
        print(h, '-dpdf', ['learning simulation results/' data_source '_' learning_type '_' file_string]);

        if strcmp(learning_type, 'triplet') || strcmp(learning_type, 'triploh')
            
            wmax = wmaxes(i_v1);
            alphaLTD = alphas_ltd(i_v2);

            saveInParfor(['learning simulation results/' data_source '_' learning_type '_' file_string '.mat'], ...
                trialW, trialFR,  noise_sig_val, wmax, alphaLTD, tau_y, i_set, i_set_string, ...
                num_training_trials, total_training_bins, total_training_spikes, param);
        else
            
            wmax = wmaxes(i_v1);
            alphaLTD = alphas_ltd(i_v2);

            saveInParfor(['learning simulation results/' data_source '_' learning_type '_' file_string '.mat'], ...
                trialW, trialFR, noise_sig_val, wmax, alphaLTD, i_set, i_set_string, ...
                num_training_trials, total_training_bins, total_training_spikes, param);
        end

    end
toc
end
end

end