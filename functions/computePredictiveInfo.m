function [projIpred, inputIpred] = computePredictiveInfo(input_words, output_proj, tau_future, nMC)
% wrapper function to simplify scripts. uses CDM entropy estimates to
% estimate word-word predictive information from inputs and output-word
% predictive information. (projIpred <= inputIpred)

if nargin < 4
    nMC = 99;   % was 999
end
% offset projection and words by tau_future (# of frames into
% future)
postwordCDM = [output_proj(1:end-tau_future) input_words((1+tau_future):end, :)];
inputinputCDM    = [input_words(1:end-tau_future, :) input_words((1+tau_future):end, :)];


% H(w_t)
m = size(input_words, 2); % get the dimension
[nn, ocnts] = words2nnOcnts(input_words);
opts.nMC = nMC; % control the # of samples from posterior
opts.verbose = false;
[S_inputs_CDM, S_inputs_CDM_var] = computeH_CDM(nn, ocnts, m, opts);

% H(output)
m = size(output_proj, 2); % get the dimension
[nn, ocnts] = words2nnOcnts(output_proj); 
opts.nMC = nMC; % control the # of samples from posterior
opts.verbose = false;
[S_post_CDM, S_post_CDM_var] = computeH_CDM(nn, ocnts, m, opts);


%H(output, w_{t+1})
m = size(postwordCDM, 2); % get the dimension
[nn, ocnts] = words2nnOcnts(postwordCDM);
opts.nMC = nMC; % control the # of samples from posterior
opts.verbose = false;
[S_postword_CDM, S_postword_CDM_var] = computeH_CDM(nn, ocnts, m, opts);

% information
projIpred(1) = S_inputs_CDM + S_post_CDM - S_postword_CDM;
projIpred(2) = sqrt(S_inputs_CDM_var + S_post_CDM_var + S_postword_CDM_var);

if nargout == 2
    
    % H(w_t, w_{t+1})    
    m = size(inputinputCDM, 2); % get the dimension
    [nn, ocnts] = words2nnOcnts(inputinputCDM);
    opts.nMC = nMC; % control the # of samples from posterior
    opts.verbose = false;
    [S_inputinput, S_inputinput_var] = computeH_CDM(nn, ocnts, m, opts);

    inputIpred(1) = 2*S_inputs_CDM - S_inputinput;
    inputIpred(2) = sqrt(2*S_inputs_CDM_var + S_inputinput_var);
end
