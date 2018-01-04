function [X_ww, X_pw, X_stim] = script4c_cellset_predictiveinfo_script(numL0_list, dir_name, old_dir_names)
% Efficient improvement on Script to measure the word-word internal and the word-stimulus predictive
% information for all cell sets.
% Written as a function to run for different set sizes, but can easily go
% back to being a script if numL0_list is set. numL0_list is a subset of
% [3 4 5 6 7 8 9 10]. 
% dir_name is the target directory to save things in
% old_dir_names is where 

data_sources = {'fishmovieA', 'fishmovieS'};


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

if nargin < 3
    if nargin < 2
        dir_name = ['information calculation results/datetime_' datestr(now, 'yyyymmdd') '/'];
    end
    if ~exist(dir_name, 'dir')
        mkdir(dir_name);
    else
        warning('results directory already exists. continuing anyway.')
    end
    old_dir_names = '';
else
    if ~exist(dir_name, 'dir')
        mkdir(dir_name);
    else
        warning('results directory already exists. continuing anyway.')
    end
            
end
for numL0 = numL0_list
    
    percX = load(['Cell Groups/perceptron_list_size' num2str(numL0)]);
    final_percs = percX.final_perceptrons;
    uniq_visited_percs = percX.unique_perceptrons;
    
 
    % save this in the directory in case the originals change. 
    saveInParfor([dir_name 'sampledPercs_setsOf' num2str(numL0)], ...
        percX);
    
    % load the partitions/perceptrons for which info calculations will be
    % performed
    if numL0 < 6
        new_parts = percX.unique_partitions;
        new_percs = percX.unique_perceptrons;
    else
        new_parts = percX.final_partitions;
        new_percs = percX.final_perceptrons;
    end
    

             
    for i_ds = 1:length(data_sources) 

        data_source = data_sources{i_ds};

        X_res = load(['aggregated results/' data_source([1 end]) num2str(numL0) '_learning_results.mat']);

        results_sets = cellfun(@(x) x.set_ID, X_res.res_ds_arr);
        results_wfinal = cellfun(@(x) squeeze(x.trialW(end, :, :)), X_res.res_ds_arr, 'uniformoutput', false);
    
        clear X_res
        i_set_list = unique(results_sets)';
        
        for i_set = i_set_list
            param = setLearningParameters('vanilla', i_set, ...
                data_source, 'numL0', numL0, 'saveandprint', false);

            file_tag = ['set' num2str(i_set) '_' data_source([1 end]) num2str(numL0)];
            
            % load learning results and find the unique learned perceptrons
            % 	stimulus information will only be calculated for these and
            %   for the optimal hull perceptrons. 
            % If the file already exists in the target directory, then don't re-run it
            if ~exist([dir_name file_tag '.mat'], 'file')
                
                % check the old_dir_names list for any files corresponding
                % to this cell set. Directories are checked in order, and
                % only the first matching file is taken. 
                i_sd = 1;
                while i_sd <= length(old_dir_names)
                    old_dir_name = old_dir_names{i_sd};
                    file_name = [old_dir_name file_tag '.mat'];
                    try X_info0 = load(file_name);
                        try X_perc0 = load([old_dir_name 'sampledPercParts_setsOf' num2str(numL0) '.mat']);
                        catch 
                            X_perc0 = load([old_dir_name 'sampledPercs_setsOf' num2str(numL0) '.mat']);
                        end
                        i_sd = length(old_dir_names) + 1;
                    catch
                        i_sd = i_sd + 1;
                        X_info0 = [];
                    end
                end
                
                
                % If X_info0 does not have X_ww, then compute X_ww
                % (word-word information)
                if ~isfield(X_info0, 'X_ww')
                    % 1) word-word- internal pred I, with lags 1:10
                    ww_internal_tau_futures = 1:10;
                    X_ww = computeInternalPredI_perceptrons([], ww_internal_tau_futures, i_set, param);
                else
                    X_ww = X_info0.X_ww;
                    
                end
                % 2) perc-word internal pred I, for all perceptrons visited by
                % learning                
                if ~isfield(X_info0, 'X_pw')
                    internal_tau_futures = [0 1];
                    if numL0 < 6
                        X_pw = computeInternalPredI_perceptrons(uniq_visited_percs, internal_tau_futures, i_set, param);
                    else
                        X_pw = computeInternalPredI_perceptrons(final_percs, internal_tau_futures, i_set, param);
                    end
                else
                    
                    % what you want to do here: look at the partitions in
                    % X_perc0 and find any new ones in the new percX sets
                    % of visited/final perceptrons. Then run the pred-I
                    % calculation for the new ones, and merge things so
                    % they are in the exact same order as the current percX
                    % partitions/perceptrons. Same general prescription for
                    % the partitions pred-I and the stimulus pred-I 
                    
                    if numL0 < 6
                        old_parts = X_perc0.percX.unique_partitions;
                    else
                        old_parts = X_perc0.percX.final_partitions;
                    end
%                     [is_new, ~] = checkPartitionUniqueness_RPMethod(...
%                         old_parts, new_parts);
%                     % this find the indices such that
%                     % old_parts = new_parts(old_in_new_inds, :)
%                     [~, old_in_new_inds] = checkPartitionUniqueness_RPMethod(...
%                         new_parts, old_parts);  
%                     % check old_in_new_inds for duplicates and set those to
%                     % zeros
%                     [uniq_oini, ind_oini] = unique(old_in_new_inds);
%                     old_in_new_inds = 0*old_in_new_inds;
%                     old_in_new_inds(ind_oini) = uniq_oini;
                    [is_new, old_in_new_inds, old_in_new_inds_dup0] = setdifference_BinaryVectors(...
                        old_parts, new_parts);

                    old_fr_perc = X_info0.X_pw.fr_perc;
                    old_pI_perc = X_info0.X_pw.pI_perc; 
                    old_dpI_perc = X_info0.X_pw.dpI_perc;
                    
                    if any(is_new)
                        internal_tau_futures = X_info0.X_pw.tau_futures;

                        % run on the new perceptrons
                        new_X_pw = computeInternalPredI_perceptrons(new_percs(is_new, :), ...
                            internal_tau_futures, i_set, param);
                        
                        full_fr_perc = zeros(size(new_percs, 1), 1);
                        full_pI_perc = zeros(size(new_percs, 1), 2);
                        full_dpI_perc = zeros(size(new_percs, 1), 2);
                        
                        
                        full_fr_perc(is_new) = new_X_pw.fr_perc;
                        full_pI_perc(is_new, :) = new_X_pw.pI_perc;
                        full_dpI_perc(is_new, :) = new_X_pw.dpI_perc;
                        
                        % zeros in the index indicate an accidental
                        % duplicate
                        old_duplicates = old_in_new_inds_dup0 == 0;
                        full_fr_perc(old_in_new_inds(~old_duplicates)) = ...
                            old_fr_perc(~old_duplicates);
                        full_pI_perc(old_in_new_inds(~old_duplicates), :) = ...
                            old_pI_perc((~old_duplicates), :);
                        full_dpI_perc(old_in_new_inds(~old_duplicates), :) = ...
                            old_dpI_perc((~old_duplicates), :);
                        
                        % X_pw is the new X_pw for all fields except
                        % fr/pI/dpI_perc. 
                        X_pw = new_X_pw;
                        X_pw.fr_perc = full_fr_perc;
                        X_pw.pI_perc = full_pI_perc;
                        X_pw.dpI_perc = full_dpI_perc;
                    else
                        internal_tau_futures = X_info0.X_pw.tau_futures;
                            % if there are no new perceptrons, reassign the
                            % old X_pw to the new X_pw. 
                        X_pw = X_info0.X_pw;
                        % ensure that the order of perceptrons stays the
                        % same as the in the percX structure. 
                        X_pw.fr_perc  = X_pw.fr_perc(old_in_new_inds);
                        X_pw.pI_perc  = X_pw.pI_perc(old_in_new_inds, :);
                        X_pw.dpI_perc = X_pw.dpI_perc(old_in_new_inds, :);
                        
                    end
    
                end
                
                
            
                stim_tau_futures = -9:2:5;

                if ~isfield(X_info0, 'X_stim')
%                     % find the learned perceptron weights
%                     set_results_wfinal = cell2mat(results_wfinal(results_sets == i_set));
% 
%                     % also find the perceptron hull weights
%                     fr_edges = linspace(0, max(X_pw.fr_perc), 201);
%                     [~, ~, ~, full_perceptron_hull_wb] = computePerceptronHull_withWeights(fr_edges, X_pw.fr_perc, X_pw.pI_perc, new_percs);
% 
%                     % combine and find optimal . 
%                     learned_and_hull_wb = [set_results_wfinal; full_perceptron_hull_wb];
%                     learned_and_hull_partitions = perceptronToProjectionRule(learned_and_hull_wb, 1);
%                     [~, set_learnedhull_inds] = returnUniquePartitionSet(learned_and_hull_partitions);
%                     stim_wb = learned_and_hull_wb(set_learnedhull_inds, :);

                    % compute for all perceptrons
                    stim_wb = new_percs;


                    % now compute stimulus information 
                    if strcmp(data_sources{i_ds}(1:4), 'fish')
                        % 4) word-(p,v) stimulus information
                        % 5) perceptron-(p,v) stimulus information for final learned
                        % perceptrons, with lags
                        comp_type = 'CDM_pv';

                        X_stim = computeBarmovieStimulusPredI(stim_wb, stim_tau_futures, i_set, comp_type, data_source);


                    else
                        X_stim = [];
                    end
                else
                    % Check perceptron uniqueness, run for specific
                    % perceptrons missing from the original run. 
                    oldX_stim = X_info0.X_stim;

                    oldStim_parts = perceptronToProjectionRule(oldX_stim.perc_ws, 1);
                    [is_new_stimP, old_in_new_inds_stimP, old_in_new_inds_stimP_dup0] = setdifference_BinaryVectors(oldStim_parts, new_parts);
                    
                    % Check tau_futures hasn't changed 
                    old_tau_futures = oldX_stim.tau_futures;
                    [~, tau_f_inds] = intersect(stim_tau_futures, old_tau_futures);

                    [missing_tau_futures, missing_tau_inds] = setdiff(stim_tau_futures, old_tau_futures);
                    
                    % fields that need to be updated: 'percstim_predI_vt', 
                    % percword_nrbpredI_vt, perc_ws, nrb_fr_perc
                    
                    if any(is_new_stimP) || ~isempty(missing_tau_futures)
                        % run for only the new percs
                        new_stim_wb = new_percs(is_new_stimP, :);
                        comp_type = 'CDM_pv';
                        
                        if ~isempty(new_stim_wb)
                            newX_stim = computeBarmovieStimulusPredI(new_stim_wb, stim_tau_futures, i_set, comp_type, data_source);
                        else 
                            % in case we are running this only to fill in
                            % the missing tau_futures, still run this for
                            % the wordstim calculations
                            newX_stim = computeBarmovieStimulusPredI(zeros(1, numL0), stim_tau_futures, i_set, comp_type, data_source);
                        end
                        %%
                        X_stim = oldX_stim;
                        % update new fields
                        % % word-word info should come from the new
                        % computation
                        X_stim.wordstim_predI_vt = newX_stim.wordstim_predI_vt;
                        X_stim.wordword_nrbpredI_vt = newX_stim.wordword_nrbpredI_vt;
                        X_stim.cellstim_predI_vt = newX_stim.cellstim_predI_vt;
                        X_stim.tau_futures = newX_stim.tau_futures;
                        %%
                        % % initialize full fields
                        X_stim.percstim_predI_vt = zeros(length(is_new_stimP), length(stim_tau_futures));
                        X_stim.percword_nrbpredI_vt = zeros(length(is_new_stimP), length(stim_tau_futures));
                        X_stim.perc_ws = zeros(length(is_new_stimP), numL0);
                        X_stim.nrb_fr_perc = zeros(length(is_new_stimP), 1);
                        %%
                        % % put the new values in
                        if any(is_new_stimP)
                            X_stim.percstim_predI_vt(is_new_stimP, :) = newX_stim.percstim_predI_vt;
                            X_stim.percword_nrbpredI_vt(is_new_stimP, :) = newX_stim.percword_nrbpredI_vt;
                            X_stim.perc_ws(is_new_stimP, :) = newX_stim.perc_ws;
                            X_stim.nrb_fr_perc(is_new_stimP) = newX_stim.nrb_fr_perc;
                        end
                        %%
                        % % put the old values in 
                        uniq_old = old_in_new_inds_stimP_dup0 ~= 0;
                        old_ord = old_in_new_inds_stimP(uniq_old);

%%
                        if ~isempty(missing_tau_futures)
                            % compute stimulus info for the missing
                            % tau_futures for the old perceptrons (was
                            % calculated above for the new perceptrons
                            % already)
                            old_percs = new_percs(~is_new_stimP, :);
                            oldX_stim_missing_tau = computeBarmovieStimulusPredI(old_percs, ...
                                missing_tau_futures, i_set, comp_type, data_source);
                            
                            
                            X_stim.percstim_predI_vt(~is_new_stimP, missing_tau_inds) = oldX_stim_missing_tau.percstim_predI_vt;
                            X_stim.percword_nrbpredI_vt(~is_new_stimP, missing_tau_inds) = oldX_stim_missing_tau.percword_nrbpredI_vt;
                        end
                        
                        %%
                        X_stim.percstim_predI_vt(old_ord, tau_f_inds) = oldX_stim.percstim_predI_vt(uniq_old, :);
                        X_stim.percword_nrbpredI_vt(old_ord, tau_f_inds) = oldX_stim.percword_nrbpredI_vt(uniq_old, :);
                        X_stim.perc_ws(old_ord, :) = oldX_stim.perc_ws(uniq_old, :);
                        X_stim.nrb_fr_perc(old_ord) = oldX_stim.nrb_fr_perc(uniq_old);  
                        
                        %%
                    elseif ~all(old_in_new_inds_stimP == sort(old_in_new_inds_stimP))
                        % verify that everything is in the same perceptron
                        % order
                        old_ord = old_in_new_inds_stimP(old_in_new_inds_stimP_dup0 ~= 0);

                        X_stim = oldX_stim;
                        X_stim.percstim_predI_vt = oldX_stim.percstim_predI_vt(old_ord, :);
                        X_stim.percword_nrbpredI_vt = oldX_stim.percword_nrbpredI_vt(old_ord, :);
                        X_stim.perc_ws = oldX_stim.perc_ws(old_ord, :);
                        X_stim.nrb_fr_perc = oldX_stim.nrb_fr_perc(old_ord);    
                    else
                        X_stim = oldX_stim;
                        
                    end
                    
                   
                    
                    
                end

                save([dir_name file_tag], 'X_ww', 'X_pw', 'X_stim');

                close all;
                
            else
                display([dir_name file_tag '.mat exists.'])
            end
        end
        
        
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
end