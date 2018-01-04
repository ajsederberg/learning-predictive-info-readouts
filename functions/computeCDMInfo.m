function [xyInfo, xEntropy, yEntropy] = computeCDMInfo(x, y, nMC)
% wrapper function to simplify scripts. uses CDM entropy estimates to
% estimate mutual info of x and y along with entropies. 

if nargin < 3
    nMC = 99;   % was 999
end


% H(x)
m = size(x, 2); % get the dimension
[nn, ocnts] = words2nnOcnts(x);
opts.nMC = nMC; % control the # of samples from posterior
opts.verbose = false;
[S_x, S_x_var] = computeH_CDM(nn, ocnts, m, opts);
xEntropy(1) = S_x;
xEntropy(2) = sqrt(S_x_var);

% H(y)
m = size(y, 2); % get the dimension
[nn, ocnts] = words2nnOcnts(y);
opts.nMC = nMC; % control the # of samples from posterior
opts.verbose = false;
[S_y, S_y_var] = computeH_CDM(nn, ocnts, m, opts);
yEntropy(1) = S_y;
yEntropy(2) = sqrt(S_y_var);

%H(output, w_{t+1})
m = size([y x], 2); % get the dimension
[nn, ocnts] = words2nnOcnts([y x]);
opts.nMC = nMC; % control the # of samples from posterior
opts.verbose = false;
[S_xy, S_xy_var] = computeH_CDM(nn, ocnts, m, opts);

% information
xyInfo(1) = S_x + S_y - S_xy;
xyInfo(2) = sqrt(S_x_var + S_y_var + S_xy_var);


