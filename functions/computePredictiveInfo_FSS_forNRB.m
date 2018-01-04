function [xyIpred, inputIpred, inputEntropy, sample_problems] = computePredictiveInfo_FSS_forNRB(x_words, y_values, tau_future, comp_type, nreps)
% wrapper function to simplify scripts. uses SEP's code for computing
% I(x;y). Key difference between this function and the one for "fishdata"
% (i.e. computePredictiveInfo and computePredictiveInf_FSS) is that this
% assumes it will receive cell array inputs in neural words and in
% stimulus values. Entries in a cell correspond to different trials, such
% that these time series may not be directly concatenated. 


% this has been modified from computePredictiveInfo_FSS to compute the
% predictive information about the stimulus, not future activity

% INPUTS:
% data_x = x data, running from 1 to bin n_x, must be integers
% data_y = y data, running from 1 to bin n_y, must be integers
% n_x = number of x bins, one dimension of Pjoint
% n_y = number of y bins, other dimension of Pjoint
% fracs = the data fractions at which to sample Pjoint, [ex: [1 0.9 0.8 0.5]]
%           one usually doesn't want to fit fractions below ~ 0.7
%           include fraction 0.5 for error estimate on I
% nreps = the number of fraction repetitions to use for finite size scaling
% varargin = seed [optional]


% [I,I_err,IR,IR_err,{I_fss,IR_fss,P,PR}] = calc_info_P_joint(data_x,data_y,n_x,n_y,fracs,nreps,varargin)
%

% OUTPUTS:
% I = information in the Pjoint distribution, finite size scaled
% I_err_std = from the 50% data fraction: std(I_samples_50%)./sqrt(2)
% I_err_frac = |I_samples_100% - I_infinity|/I_samples_100%
% IR, IR_err_std, IR_frac: same as I, but for the shuffled data (random shuffled of bin
% assignments, tests sampling problems.  If IR +/- IR_err*2 ~= 0, worry.


if strcmp(comp_type, 'SEP')
%% step 1: get the tau_future offset between neural words and stimulus tokens, 
% then combine across the cell array into a time series. 

data_x_cellarr = cell(length(x_words), 1);
data_y_cellarr = cell(length(x_words), 1);

for i_trial = 1:length(x_words)
    % data_x is the readout of each learned perceptron at time t
    if tau_future > 0
        d_x = x_words{i_trial}(1:(end-tau_future), :) + 1;
    else
        d_x = x_words{i_trial}((1-tau_future):(end), :) + 1;
    end
    data_x_cellarr{i_trial} = d_x;
    
    % data_y is the stimulus at time t + tau_future. 
    if tau_future > 0
        d_y = y_values{i_trial,1}((1+tau_future):end, :);
    else
        d_y = y_values{i_trial,1}(1:(end+tau_future), :);
    end
    
    data_y_cellarr{i_trial} = d_y;
    
    
    if numel(data_y_cellarr{i_trial}) == 0
        data_y_cellarr{i_trial} = [];
    end

end

data_x = cell2mat(data_x_cellarr);
n_x = max(data_x(:));



% data_y is the (#timepoints) by MM representation of
% the stimulus. This needs to be ocnverted to integers as well. This only
% has to be done once, so we do it here. 
data_y = cell2mat(data_y_cellarr);


conv_factor = ((max(1+data_y(:))).^(0:(size(data_y,2)-1)))';

data_y_nums = data_y*conv_factor;
[~, ~, data_y] = unique(data_y_nums);
n_y = max(data_y);


n_y = max(data_y);

fracs = [1 0.9 0.8 0.5];
if nargin > 4 && isnan(nreps)
    only_compute_input_predI = true;
else
    only_compute_input_predI = false;
end

if nargin <= 4 || isnan(nreps)
    
    nreps = 50;
end

if ~only_compute_input_predI
    xyIpred = zeros(size(data_x, 2) + 1, 3);
    sample_problems = zeros(size(data_x, 2), 1);
    for i_out = 1:size(data_x, 2)
        data_x1 = data_x(:, i_out);
        [I,I_err_std,I_err_frac,IR,IR_err_std, ~]= calc_info_P_joint(data_x1,data_y,n_x,n_y,fracs,nreps);

        xyIpred(i_out, :) = [I I_err_std I_err_frac];
        % If IR +/- IR_err*2 ~= 0, worry.
        sample_problems(i_out) = IR > IR_err_std*2 || IR < -IR_err_std*2;
    end
    
    
    
    % compute for full word, for up to 20 cells (or 20 perceptrons)
    if size(data_x, 2) <= 20
        conv_factor = (2.^(0:1:(size(data_x, 2)-1)))';

        input_nums = (data_x - 1)*conv_factor;
        [~, ~, data_x] = unique(input_nums);
        n_x = max(data_x(:));

        [I,I_err_std,I_err_frac,IR,IR_err_std, ~]= calc_info_P_joint(data_x,data_y,n_x,n_y,fracs,nreps);

        xyIpred(end, :) = [I I_err_std I_err_frac];
        % If IR +/- IR_err*2 ~= 0, worry.
        sample_problems(end) = IR > IR_err_std*2 || IR < -IR_err_std*2;
    end

    
    
else
    xyIpred = [];
    sample_problems = [];
end
% if nargout > 1
%     [I,I_err_std,I_err_frac, IR,IR_err_std, ~]= calc_info_P_joint(data_x1,data_x2,n_x,n_x,fracs,nreps);
%     inputIpred = [I I_err_std I_err_frac];
%     if only_compute_input_predI
%         sample_problems = IR > IR_err_std*2 || IR < -IR_err_std*2;
%     end
% end
% 
% if nargout > 2
%     [I,I_err_std,I_err_frac, DUMMYA, DUMMYB, DUMMYC]= calc_info_P_joint(data_x1,data_x1,n_x,n_x,fracs,nreps);
%     inputEntropy = [I I_err_std I_err_frac];
% end

inputEntropy = [];
inputIpred = [];

elseif strcmp(comp_type, 'SEP_pv')
% %     UPDATE THIS CODE TO COMPUTE INFO ABOUT POSITION AND VELOCITY
% % NEEDS SDEBUGGING STARTING HERE
%% step 1: get the tau_future offset between neural words and stimulus tokens, 
% then combine across the cell array into a time series. 
data_x_cellarr = cell(length(x_words), 1);
data_y_cellarr = cell(length(x_words), 1);
%%
for i_trial = 1:length(x_words)
    % data_x is the readout of each learned perceptron at time t
    if tau_future > 0
        d_x = x_words{i_trial}(1:(end-tau_future-1), :) + 1;
    else
        d_x = x_words{i_trial}((1-tau_future):(end-1), :) + 1;
    end
    data_x_cellarr{i_trial} = d_x;
    
    % data_y is the stimulus at time t + tau_future. 
    if tau_future > 0
        d_y = y_values{i_trial,1}((1+tau_future):(end-1), :);
        v_y = y_values{i_trial,2}((1+tau_future):(end-1), :);
    else
        d_y = y_values{i_trial,1}(1:(end+tau_future-1), :);
        v_y = y_values{i_trial,2}(1:(end+tau_future-1), :);
    end
    
    data_y_cellarr{i_trial} = [d_y v_y];
    
    
    if numel(data_y_cellarr{i_trial}) == 0
        data_y_cellarr{i_trial} = [];
    end
    
    
end
%% 
% data_x is spiking patterns, either (#time points) by (# readouts) or
% (#timepoints) by (#cells per word). These have to be converted to
% integers prior to being passed to calc_info_P_joint
data_x = cell2mat(data_x_cellarr);
n_x = max(data_x(:));

% data_y is the (#timepoints) by 2 (position, velocity) representation of
% the stimulus. This needs to be ocnverted to integers as well. This only
% has to be done once, so we do it here. 
data_y = cell2mat(data_y_cellarr);



conv_factor = [1; max(data_y(:))];

data_y_nums = data_y*conv_factor;
[~, ~, data_y] = unique(data_y_nums);
n_y = max(data_y);


%%

fracs = [1 0.9 0.8 0.5];
if nargin > 4 && isnan(nreps)
    only_compute_input_predI = true;
else
    only_compute_input_predI = false;
end

if nargin <= 4 || isnan(nreps)
    
    nreps = 50;
end

if ~only_compute_input_predI
    xyIpred = zeros(size(data_x, 2) + 1, 3);
    sample_problems = zeros(size(data_x, 2) + 1, 1);
    for i_out = 1:size(data_x, 2)
        data_x1 = data_x(:, i_out);
        [I,I_err_std,I_err_frac,IR,IR_err_std, ~]= calc_info_P_joint(data_x1,data_y,n_x,n_y,fracs,nreps);

        xyIpred(i_out, :) = [I I_err_std I_err_frac];
        % If IR +/- IR_err*2 ~= 0, worry.
        sample_problems(i_out) = IR > IR_err_std*2 || IR < -IR_err_std*2;
    end
    
    
    
    % compute for full word
    conv_factor = (2.^(0:1:(size(data_x, 2)-1)))';

    input_nums = (data_x - 1)*conv_factor;
    [~, ~, data_x] = unique(input_nums);
    n_x = max(data_x(:));

    [I,I_err_std,I_err_frac,IR,IR_err_std, ~]= calc_info_P_joint(data_x,data_y,n_x,n_y,fracs,nreps);

    xyIpred(end, :) = [I I_err_std I_err_frac];
    % If IR +/- IR_err*2 ~= 0, worry.
    sample_problems(end) = IR > IR_err_std*2 || IR < -IR_err_std*2;
    
    
    
else
    xyIpred = [];
    sample_problems = [];
end
% if nargout > 1
%     [I,I_err_std,I_err_frac, IR,IR_err_std, ~]= calc_info_P_joint(data_x1,data_x2,n_x,n_x,fracs,nreps);
%     inputIpred = [I I_err_std I_err_frac];
%     if only_compute_input_predI
%         sample_problems = IR > IR_err_std*2 || IR < -IR_err_std*2;
%     end
% end
% 
% if nargout > 2
%     [I,I_err_std,I_err_frac, DUMMYA, DUMMYB, DUMMYC]= calc_info_P_joint(data_x1,data_x1,n_x,n_x,fracs,nreps);
%     inputEntropy = [I I_err_std I_err_frac];
% end

inputEntropy = [];
inputIpred = [];

elseif strcmp(comp_type, 'CDM')
    % slightly different format required than for SEP code: x and y need to be binary 
    % we assume here that x_words IS ALREADY A BINARY WORD, but y_values
    % may be binary or not
    num_bits_y = max(cellfun(@(x) max(ceil(log2(x(:)))), y_values));
    data_x_cellarr = cell(length(x_words), 1);
    data_y_cellarr = cell(length(x_words), 1);

    for i_trial = 1:length(x_words)
        % data_x is the readout of each learned perceptron at time t: this
        % is a 0 or 1 and is already in the correct format for CDM
        % calculations
        if tau_future > 0
            d_x = x_words{i_trial}(1:(end-tau_future), :);
        else
            d_x = x_words{i_trial}((1-tau_future):end, :);
        end
        data_x_cellarr{i_trial} = d_x;

        
        
        % data_y is at time t + tau_future. May need to
        % convert to binary representation. 
        if tau_future > 0
            d_y = y_values{i_trial}((1+tau_future):end, :);
        else
            d_y = y_values{i_trial}(1:(end+tau_future), :);
        end
        % this converts to binary words if y was given in base-10
        if max(d_y) > 1 
            data_y_cellarr{i_trial} = dec2binMATRIX(d_y, num_bits_y);
        else
            data_y_cellarr{i_trial} = d_y;
        end
        
        if numel(data_y_cellarr{i_trial}) == 0
            data_y_cellarr{i_trial} = [];
        end
    end

    data_x = cell2mat(data_x_cellarr); 
    data_y = cell2mat(data_y_cellarr);
    
    xyIpred = zeros(size(data_x, 2)+1, 2);
    inputEntropy = zeros(size(data_x, 2)+1, 2);
    for i_out = 1:size(data_x, 2)
        try
            [xyInfo, xEntropy, ~] = computeCDMInfo(data_x(:, i_out), data_y);
            xyIpred(i_out, :) = xyInfo;
            inputEntropy(i_out, :) = xEntropy;
        catch
            
            xyIpred(i_out, :) = nan;
            inputEntropy(i_out, :) = nan;            
            
        end
    end
    
    try
        [xyInfo, xEntropy] = computeCDMInfo(data_x, data_y);
        xyIpred(end, :) = xyInfo;
        inputEntropy(end, :) = xEntropy;  
    catch
        xyIpred(end, :) = nan;
        inputEntropy(end, :) = nan;  
    end
  
    inputIpred = [];
    sample_problems = [];
    
elseif strcmp(comp_type, 'CDM_pv')
    % slightly different format required than for SEP code: x and y need to be binary 
    % we assume here that x_words IS ALREADY A BINARY WORD, but y_values
    % may be binary or not 
    % CDM_pv computes I(output_t; {position_t+1, velocity_t+1})
    num_bits_y = max(cellfun(@(x) max(ceil(log2(x))), y_values));
    data_x_cellarr = cell(length(x_words), 1);
    data_y_cellarr = cell(length(x_words), 1);
%%
    for i_trial = 1:length(x_words)
        % data_x is the readout of each learned perceptron at time t: this
        % is a 0 or 1 and is already in the correct format for CDM
        % calculations
        if tau_future > 0
            d_x = x_words{i_trial}(1:(end-tau_future-1), :);
        else
            d_x = x_words{i_trial}((1-tau_future):(end-1), :);
        end
        data_x_cellarr{i_trial} = d_x;

        
        
        % data_y is at time t + tau_future. May need to
        % convert to binary representation. 
        if tau_future > 0
            d_y = y_values{i_trial,1}((1+tau_future):(end-1), :);
            v_y = y_values{i_trial,2}((1+tau_future):(end-1), :);
        else
            d_y = y_values{i_trial,1}(1:(end+tau_future-1), :);
            v_y = y_values{i_trial,2}(1:(end+tau_future-1), :);
        end
        % this converts to binary words if y was given in base-10
        if max(d_y) > 1 
            d_mat = dec2binMATRIX(d_y, num_bits_y(1));
            v_mat = dec2binMATRIX(v_y, num_bits_y(1));
            data_y_cellarr{i_trial} =  [d_mat v_mat];
        else
            data_y_cellarr{i_trial} = [d_y v_y];
        end
        
        if numel(data_y_cellarr{i_trial}) == 0
            data_y_cellarr{i_trial} = [];
        end
    end
%%
    data_x = cell2mat(data_x_cellarr); 
    data_y = cell2mat(data_y_cellarr);
    
    xyIpred = zeros(size(data_x, 2)+1, 2);
    inputEntropy = zeros(size(data_x, 2)+1, 2);
    
    try 
        [xyInfo, xEntropy] = computeCDMInfo(data_x, data_y);
        xyIpred(end, :) = xyInfo;
        inputEntropy(end, :) = xEntropy;    
    catch
        xyIpred(end, :) = nan;
        inputEntropy(end, :) = nan;
    end
    
    
    
    for i_out = 1:size(data_x, 2)
    
        try 
            [xyInfo, xEntropy, ~] = computeCDMInfo(data_x(:, i_out), data_y);
            xyIpred(i_out, :) = xyInfo;
            inputEntropy(i_out, :) = xEntropy;
        catch
            xyIpred(i_out, :) = nan;
            inputEntropy(i_out, :) = nan;
            
        end
        
    end 
    

    inputIpred = [];
    sample_problems = [];
end