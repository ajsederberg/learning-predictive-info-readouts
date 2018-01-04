function [X_tokenw, X_tokenstim] = script4d_cellset_predictiveinfo_script(numL0_list, dir_name)
% Script to measure token information. Based on the efficient cellset predI
% script. 


data_sources = {'fishmovieA', 'fishmovieS'}; %, 'multimovieS', 'multimovieA'};



% Calculations required:
% computeInternalPredI function:
%  1) Word-word internal pred I, with lags 1:10
%  2) Perceptron-word internal pred I, for all perceptrons visited by
%  learning (1 lag only)
%  3) Partition-word internal pred I, for sampled partitions
% computeBarmovieStimulusPredI: 
%  4) Word-(p,v) stimulus information, with lags -10:1:10 [for fA and fS
%  sets only)
%  5) Perceptron-(p,v) stimulus information, for final learned perceptrons,
%  with lags -10:1:10. [for fA and fS sets only]


if nargin < 2
    dir_name = ['information calculation results/datetime_' datestr(now, 'yyyymmdd') '/'];
end
if ~exist(dir_name, 'dir')
    mkdir(dir_name);
else
    warning('results directory already exists. continuing anyway.')
end


for numL0 = numL0_list
    

    for i_ds = 1:length(data_sources) 

        data_source = data_sources{i_ds};

        X_res = load(['aggregated results/' data_source([1 end]) num2str(numL0) '_learning_results.mat']);

        results_sets = cellfun(@(x) x.set_ID, X_res.res_ds_arr);
    
        clear X_res
        i_set_list = unique(results_sets)';
        
        for i_set = i_set_list
            param = setLearningParameters('vanilla', i_set, ...
                data_source, 'numL0', numL0, 'saveandprint', false);

            file_tag = ['tokeninfo_set' num2str(i_set) '_' data_source([1 end]) num2str(numL0)];
            
            if ~exist([dir_name file_tag '.mat'], 'file')

            internal_tau_futures = [0 1];

            X_tokenw = computeInternalPredI_perceptrons('words', internal_tau_futures, i_set, param);
            
            % now compute stimulus information 
            if strcmp(data_sources{i_ds}(1:4), 'fish')
                % 4) word-(p,v) stimulus information
                % 5) perceptron-(p,v) stimulus information for final learned
                % perceptrons, with lags
                stim_tau_futures = 0:1;
                comp_type = 'CDM_pv';

                X_tokenstim = computeBarmovieStimulusPredI(nan(1, numL0), stim_tau_futures, i_set, comp_type, data_source);


            else
                X_tokenstim = [];
            end


            save([dir_name file_tag], 'X_tokenw', 'X_tokenstim');

            close all;
            end

        end
        
        
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
end

                

% 
% learned_freq: [0.5550 0.1050 0.1250 0.1100 0.1050]
% learned_percword_predI_fishmovie: [5x1 double]
% learned_percword_predI_nrb: [5x1 double]
% learned_percstim_predI_nrb: [5x21 double]
% learned_perc_fr_fishmovie: [5x1 double]
% learned_perc_fr_nrb: [5x1 double]
% optimal_percword_predI_fishmovie: [5x1 double]
% optimal_percword_predI_nrb: [5x1 double]
% optimal_percstim_predI_nrb: [5x21 double]
% optimal_perc_fr_fishmovie: [5x1 double]
% optimal_perc_fr_nrb: [5x1 double]
% group_predI_fishmovie: [0.0578 0.0106]
% group_predI_nrb: [1x21 double]
% group_fr_fishmovie: 0.0488
% group_fr_nrb: 0.0964
% set_num: 10