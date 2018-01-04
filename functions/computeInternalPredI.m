function [Xout, file_name] = computeInternalPredI(sampled_partitions, tau_futures, set_num, param)
% script to calculate predictive info landscape 

cellset_dir = 'cell sets multi and fishmovie/';

numL0 = param.numL0;
    
num_partitions_sampled = size(sampled_partitions, 1);
fr_proj = zeros(num_partitions_sampled, 1);
pI_proj  = zeros(num_partitions_sampled, length(tau_futures));
dpI_proj = zeros(num_partitions_sampled, length(tau_futures));

[~, full_word_record] = loadDataForLearningSimulation(set_num, param);

spikingbin_times = find(sum(full_word_record, 1) ~= 0);
word_record = full_word_record(:, spikingbin_times);
conv_factor = 2.^(0:(numL0-1));
word_numbers = conv_factor*full_word_record;
word_nums = nonzeros(word_numbers);

% this is a (# non-zero spiking words) by (number of spiking times)
% representation of full-word spikinga ctivity. The spike times are in
% 'spikingbing_tiems'
sparse_wordmat = sparse(word_nums, 1:length(word_nums), 1, ...
        2^numL0 - 1, length(word_nums));

tic;


for i_tau = 1:length(tau_futures)
    tau_future = tau_futures(i_tau);

    sliced_sampled_partitions = sampled_partitions(:, 2:end)';
    
    parfor i_sps = 1:num_partitions_sampled


        % sparse wordmat has a 1 whenever a particular pattern is observed.
        % This collapses that into a single 0/1 based on whether that
        % pattern is a 0/1 in the sampled_partition vector. 
        % This drops the first word (the 0000 word)
        wij_proj = sliced_sampled_partitions(:, i_sps)'; %sampled_partitions(i_sps, 2:end);
        proj_active = wij_proj*sparse_wordmat;
        proj_bintimes = spikingbin_times(proj_active == 1);


        % use CDM estimator to get word-proj predictive information
        wordsCDM = zeros(max(spikingbin_times), numL0);
        wordsCDM(spikingbin_times, :) = word_record';

        projCDM = zeros(max(spikingbin_times), 1);
        projCDM(proj_bintimes) = 1;

        fr_proj(i_sps) = mean(projCDM);

        pI_out = computePredictiveInfo(wordsCDM, projCDM, tau_future);
        pI_proj(i_sps, i_tau) = pI_out(1);
        dpI_proj(i_sps, i_tau) = pI_out(2);



    end


end

        
% close the pool to prevent hang-ups (may not be a problem in
% 2016a, was a problem in 2014b and before.
delete(gcp('nocreate'))

display(['set ' num2str(set_num) ' partitions done, ' datestr(now, 'HH:MM')])
   
file_tag =  ['partitionWord_IIp_set' num2str(set_num) '_' ...
    param.data_source([1 end]) num2str(numL0)];

if nargout == 2
    file_name = [cellset_dir file_tag];
    saveInParfor(file_name, fr_proj, pI_proj, dpI_proj, set_num, tau_futures);
end

% Everything except the inputs go in the output struct - this is for the
% sake of smaller datasets.
Xout.set_num = set_num; 
Xout.param = param; 
Xout.file_tag = file_tag;
Xout.tau_futures = tau_futures;
Xout.fr_proj = fr_proj;
Xout.pI_proj = pI_proj; 
Xout.dpI_proj = dpI_proj;
