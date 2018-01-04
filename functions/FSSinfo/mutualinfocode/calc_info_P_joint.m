function [I,I_err_std,I_err_frac,IR,IR_err_std,IR_err_frac,varargout] = calc_info_P_joint(data_x,data_y,n_x,n_y,fracs,nreps,varargin)
% [I,I_err,IR,IR_err,{I_fss,IR_fss,P,PR}] = calc_info_P_joint(data_x,data_y,n_x,n_y,fracs,nreps,varargin)
% This function calculates the information in the joint distribution of 
% the x and y data, using finite size scaling.
%
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
%
% OUTPUTS:
% I = information in the Pjoint distribution, finite size scaled
% I_err_std = from the 50% data fraction: std(I_samples_50%)./sqrt(2)
% I_err_frac = |I_samples_100% - I_infinity|/I_samples_100%
% IR, IR_err_std, IR_frac: same as I, but for the shuffled data (random shuffled of bin
% assignments, tests sampling problems.  If IR +/- IR_err*2 ~= 0, worry.

% {...} = varargout = a number of optional output arguments:
% I_fss = the data from all the samples taken of P at different data fractions
% IR_fss = same as above, but for the randomly shuffled data
% P = the joint distribution, for all the data
% PR = same, but for the shuffled data

if nargin == 6
    seed = int32(rem(now,1)*1000000);
else
    seed = varargin{1}; %#ok<*ASGSL>
end

nfracreps = nreps;
ndata = numel(data_x);

fracs = sort(fracs,'descend');
if ~ismember(fracs,0.5)
    fracs = [fracs 0.5];
    fracs = sort(fracs,'descend');
end
nstop = find(fracs == 0.5);

nfracs = length(fracs);

%*% check inputs to c-code
if numel(data_x) ~= numel(data_y)
    error('sepFunction:argCheck','data_x and data_y must be the same length');
end

if min(data_x) < 1 || min(data_y) < 1
    error('sepFunction:argCheck','data values must start at 1');
end

if max(data_x) > n_x || max(data_y) > n_y
    error('sepFunction:argCheck','data values must not exceed n_x or n_y');
end

if sum(rem(data_x,1)) ~= 0 || sum(rem(data_y,1)) ~= 0
    error('sepFunction:argCheck','data must be integers');
end

if nargout < 6
    error('sepFunction:argCheck','need at least 6 output arguments');
end

n_xvals = n_x;
n_yvals = n_y;

[Ip,IpR,Pp,PpR]=calc_info_pjoint(data_x,data_y,fracs);

temp = polyfit(1./fracs(1:nstop),mean(Ip(:,1:nstop),1),1);
I_fss = Ip;
I = temp(2);
I_err_std = std(Ip(:,end))/sqrt(2);
I_err_frac = abs(Ip(1,1)-I)/Ip(1,1);

temp = polyfit(1./fracs(1:nstop),mean(IpR(:,1:nstop),1),1);
IR_fss = IpR;
IR = temp(2);
IR_err_std = std(IpR(:,end))/sqrt(2);
IR_err_frac = abs(IpR(1,1)-IR)/IpR(1,1);

varargout(1) = {I_fss};
varargout(2) = {IR_fss};
varargout(3) = {Pp};
varargout(4) = {PpR};