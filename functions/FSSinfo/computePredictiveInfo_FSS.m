function [projIpred, inputIpred, inputEntropy, sample_problems] = computePredictiveInfo_FSS(input_words, output_proj, tau_future, nreps)
% wrapper function to simplify scripts. uses SEP's code for computing
% I(x;y). Takes inputs in the same format as my 'computePredictiveInfo' and
% spits out in the same format as well. 

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




% this converts input_words to data_x by computing the base-2
% representating of the input word, then taking only the unique words. 
num_cells = size(input_words, 2);

conv_factor = (2.^(0:1:(num_cells-1)))';

input_nums = input_words*conv_factor;
[~, ~, data_x] = unique(input_nums);
n_x = max(data_x);

% output_proj is a vector of 0's and 1's for the output y = f(x_i). Add one
% so that data_y consists of 1's and 2's. 
data_y = output_proj + 1;
n_y = 2;

% finally, we want the predictive info, so take data_x at time t+tau_future
% and data_y at time t
data_x1 = data_x((1+tau_future):end);
data_x2 = data_x(1:(end-tau_future));
data_y2 = data_y(1:(end-tau_future));

fracs = [1 0.9 0.8 0.5];
if nargin >=4 && isnan(nreps)
    only_compute_input_predI = true;
else
    only_compute_input_predI = false;
end

if nargin < 4 || isnan(nreps)
    
    nreps = 50;
end

if ~only_compute_input_predI
    [I,I_err_std,I_err_frac,IR,IR_err_std, ~]= calc_info_P_joint(data_x1,data_y2,n_x,n_y,fracs,nreps);

    projIpred = [I I_err_std I_err_frac];
    % If IR +/- IR_err*2 ~= 0, worry.
    sample_problems = IR > IR_err_std*2 || IR < -IR_err_std*2;
else
    projIpred = [];
    sample_problems = [];
end
if nargout > 1
    [I,I_err_std,I_err_frac, IR,IR_err_std, ~]= calc_info_P_joint(data_x1,data_x2,n_x,n_x,fracs,nreps);
    inputIpred = [I I_err_std I_err_frac];
    if only_compute_input_predI
        sample_problems = IR > IR_err_std*2 || IR < -IR_err_std*2;
    end
end

if nargout > 2
    [I,I_err_std,I_err_frac, DUMMYA, DUMMYB, DUMMYC]= calc_info_P_joint(data_x1,data_x1,n_x,n_x,fracs,nreps);
    inputEntropy = [I I_err_std I_err_frac];
end