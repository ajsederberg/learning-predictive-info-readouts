function [X_ww, X_stim] = script4b_cellset_predictiveinfo_script(numL0_list, i_set_list, dir_name)
% Simplified script to measure the word-word internal and the word-stimulus predictive
% information for all cell sets.
% Written as a function to run for different set sizes, but can easily go
% back to being a script if numL0_list is set. numL0_list is a subset of
% [3 4 5 6 7 8 9 10]. 


data_sources = {'fishmovieA', 'fishmovieS'}; %, 'multimovieS', 'multimovieA'};



% Calculations required:
% computeInternalPredI function:
%  1) Word-word internal pred I, with lags 1:10

% computeBarmovieStimulusPredI: 
%  4) Word-(p,v) stimulus information, with lags -10:1:10 [for fA and fS
%  sets only)

if nargin < 3
    dir_name = ['information calculation results/datetime_' datestr(now, 'yyyymmdd') '/'];
    if ~exist(dir_name, 'dir')
        mkdir(dir_name);
    else
        warning('results directory already exists. continuing anyway.')
    end
else
    display(['saving results to ' dir_name])
end
for numL0 = numL0_list    
    for i_set = i_set_list
         for i_ds = 1:length(data_sources)

            data_source = data_sources{i_ds};


            param = setLearningParameters('vanilla', i_set, ...
                data_source, 'numL0', numL0, 'saveandprint', false);

            file_tag = ['wordword_set' num2str(i_set) '_' data_source([1 end]) num2str(numL0)];

            % this checks for an existing file in the directory. run_ww and
            % run_stim keep track of whether the ww and stim info
            % calcualtions have already been run
            if exist([dir_name file_tag '.mat'], 'file')
                already_run = load([dir_name file_tag]);
                run_ww = ~isfield(already_run, 'X_ww');
                run_stim = ~isfield(already_run, 'X_stim');
                
            else
                run_ww = true;
                run_stim = true;
            end
            
            if run_ww    
                % 1) word-word- internal pred I, with lags 1:10
                ww_internal_tau_futures = 1:10;
                X_ww = computeInternalPredI_perceptrons([], ww_internal_tau_futures, i_set, param);
            else
                % if ww info was already run, load from saved file. 
                X_ww = already_run.X_ww;
            end
            
            if run_stim

                if strcmp(data_sources{i_ds}(1:4), 'fish')
                    % 4) word-(p,v) stimulus information
                    % 5) perceptron-(p,v) stimulus information for final learned
                    % perceptrons, with lags
                    stim_tau_futures = -10:1:10;
                    comp_type = 'CDM_pv';
                    one_perc = ones(1, numL0);
                    X_stim = computeBarmovieStimulusPredI(one_perc, stim_tau_futures, i_set, comp_type, data_source);


                else
                    X_stim = [];
                end
            else
                % if stim info was already run, load from saved file. 
                X_stim = already_run.X_stim;
            end
            if run_stim || run_ww
                if exist([dir_name file_tag], 'file')
                    save([dir_name file_tag], 'X_ww', 'X_stim', '-append');
                    display(['appended for set ' file_tag])
                else
                    save([dir_name file_tag], 'X_ww', 'X_stim');
                    display(['saved set ' file_tag])
                end
            end
            close all;
        end
    end
end