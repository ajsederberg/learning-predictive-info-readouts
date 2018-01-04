function [Xout, file_name] = computeInternalPredI_perceptrons(perc_ws, tau_futures, set_num, param)
% script to calculate predictive info landscape for the set of supplied
% perceptrons. 

cellset_dir = 'information calculation results/';


numL0 = param.numL0;
    
Xout = [];
[~, full_word_record] = loadDataForLearningSimulation(set_num, param);

full_word_record = full_word_record';
word_fr = mean(full_word_record);

if ~isempty(perc_ws)
    if strcmp(perc_ws, 'words')
        % compute the spikeword - full word information 
        word_record_numbers = full_word_record*(2.^(0:1:(numL0-1))');
        [word_labels, wl_inds] = unique(word_record_numbers);
        
        word_tokens = full_word_record(wl_inds, :);
        keep_words = sum(word_tokens, 2) > 0 & sum(word_tokens, 2) < 4;
        word_labels = word_labels(keep_words);
        word_tokens = word_tokens(keep_words, :);
        
        perc_out = zeros(size(full_word_record, 1), length(word_labels));
        for i_wt = 1:length(word_labels)
            perc_out(word_record_numbers == word_labels(i_wt), i_wt) = 1;
        end
        fr_perc = mean(perc_out)';
        
        file_tag =  ['tokenWord_IIp_set' num2str(set_num) '_' param.data_source([1 end]) num2str(numL0)];
        num_perceptrons_samples = size(perc_out, 2);
        pI_percT  = zeros(num_perceptrons_samples, length(tau_futures));
        dpI_percT = zeros(num_perceptrons_samples, length(tau_futures));
        
        Xout.word_tokens = word_tokens;
        
    else
        % if perceptrons are supplied, then compute internal pred-I for those
        file_tag =  ['percWord_IIp_set' num2str(set_num) '_' param.data_source([1 end]) num2str(numL0)];

        num_perceptrons_samples = size(perc_ws, 1);
        pI_percT  = zeros(num_perceptrons_samples, length(tau_futures));
        dpI_percT = zeros(num_perceptrons_samples, length(tau_futures));

        b = 1;
        perc_out = double(full_word_record*perc_ws' > b);
        fr_perc = mean(perc_out)';   
    end
else
    % if no perceptrons are supplied, then compute internal word-word info
    num_perceptrons_samples = 1;
    word_pI = zeros(1, length(tau_futures));
    dword_pI = zeros(1, length(tau_futures));
    file_tag =  ['wordWord_IIp_set' num2str(set_num) '_' param.data_source([1 end]) num2str(numL0)];

end



% compute predictive information of all perceptrons

for i_tau = 1:length(tau_futures)
    tau_future = tau_futures(i_tau);
    if ~isempty(perc_ws)
        parfor ii = 1:num_perceptrons_samples
            pI_out = computePredictiveInfo(full_word_record, perc_out(:, ii), tau_future);
            pI_percT(ii, i_tau) = pI_out(1);
            dpI_percT(ii, i_tau) = pI_out(2);
        end
    
    else
        % compute word-word pred-I only
        [~, pI_word] = computePredictiveInfo(full_word_record, 0*full_word_record(:, 1), tau_future);

        word_pI(i_tau) = pI_word(1);
        dword_pI(i_tau) = pI_word(2);
            
    end
end
display(['set ' num2str(set_num) ' perceptrons done, ' datestr(now, 'HH:MM')])
%%

if nargout == 2
    file_name = [cellset_dir file_tag];
    if ~isempty(file_name)
        saveInParfor(file_name, set_num, fr_perc, pI_percT, dpI_percT, ...
                tau_futures, param, perc_ws, b);
    else
        saveInParfor(file_name, set_num, word_pI, dword_pI, word_fr, ...
                tau_futures, param);        

    end
end

% Everything except the inputs goes in the output struct to minimize file
% size
Xout.set_num = set_num; 
Xout.param = param; 
Xout.file_tag = file_tag;
Xout.tau_futures = tau_futures;
if ~isempty(perc_ws)
    Xout.fr_perc = fr_perc;
    Xout.pI_perc = pI_percT; 
    Xout.dpI_perc = dpI_percT;
    
                        
    % close the pool to prevent hang-ups (may not be a problem in
    % 2016a, was a problem in 2014b and before.
    delete(gcp('nocreate'))
else
    Xout.word_pI = word_pI;
    Xout.dword_pI = dword_pI;
    Xout.word_fr = word_fr;            
end