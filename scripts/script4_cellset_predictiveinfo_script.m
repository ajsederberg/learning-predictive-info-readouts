function [X_ww, X_pw, X_partitionword, X_stim] = script4_cellset_predictiveinfo_script(numL0_list, i_set_list, compute_for_partitions, dir_name)
% Script to measure the word-word internal and the word-stimulus predictive
% information for all cell sets.
% Written as a function to run for different set sizes, but can easily go
% back to being a script if numL0_list is set. numL0_list is a subset of
% [3 4 5 6 7 8 9 10]. 


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

if nargin < 4
    dir_name = ['information calculation results/datetime_' datestr(now, 'yyyymmdd') '/'];
    if ~exist(dir_name, 'dir')
        mkdir(dir_name);
    else
        warning('results directory already exists. continuing anyway.')
    end
    old_dir_name = '';
else
    if ~exist(dir_name, 'dir')
        mkdir(dir_name);
    else
        % check for existing partitions/perceptrons in the directory. 
        old_dir_name = dir_name;
        % make new directory
        dir_name = ['information calculation results/datetime_' datestr(now, 'yyyymmdd') '/'];
        mkdir(dir_name);

    end    
    
end
for numL0 = numL0_list
    
    percX = load(['Cell Groups/perceptron_list_size' num2str(numL0)]);
    final_percs = percX.final_perceptrons;
    uniq_visited_percs = percX.unique_perceptrons;
    
    partX = load(['Cell Groups/partition_list_size' num2str(numL0)]);
    sampled_partitions = partX.sampled_partitions;
    
    % save this in the directory in case the originals change. 
    saveInParfor([dir_name 'sampledPercParts_setsOf' num2str(numL0)], ...
        percX, partX);
    
    for i_set = i_set_list
             
        for i_ds = 1:length(data_sources)

            data_source = data_sources{i_ds};


            param = setLearningParameters('vanilla', i_set, ...
                data_source, 'numL0', numL0, 'saveandprint', false);

            file_tag = ['set' num2str(i_set) '_' data_source([1 end]) num2str(numL0)];

            if ~exist([dir_name file_tag], 'file')
                
                if exist([old_dir_name file_tag '.mat'], 'file')
                    X_info0 = load([old_dir_name file_tag]);
                    X_perc0 = load([old_dir_name 'sampledPercParts_setsOf' num2str(numL0) '.mat']);
                else
                    X_info0 = [];
                end
                
                
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
                        new_parts = percX.unique_partitions;
                        new_percs = percX.unique_perceptrons;
                    else
                        old_parts = X_perc0.percX.final_partitions;
                        new_parts = percX.final_partitions;
                        new_percs = percX.final_perceptrons;
                        
                    end
                    [is_new, ~] = checkPartitionUniqueness_RPMethod(...
                        old_parts, new_parts);
                    % this find the indices such that
                    % old_parts = new_parts(old_in_new_inds, :)
                    [~, old_in_new_inds] = checkPartitionUniqueness_RPMethod(...
                        new_parts, old_parts);  

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
                        
                        [~, new_in_new_inds] = checkPartitionUniqueness_RPMethod(...
                            new_parts, new_parts(is_new, :));
                        full_fr_perc(new_in_new_inds) = new_X_pw.fr_perc;
                        full_pI_perc(new_in_new_inds, :) = new_X_pw.pI_perc;
                        full_dpI_perc(new_in_new_inds, :) = new_X_pw.dpI_perc;
                        
                        % zeros in the index indicate an accidental
                        % duplicate
                        old_duplicates = old_in_new_inds == 0;
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
                
                
                if compute_for_partitions
                % 3) partition-word internal pred I, for sampled partitions
                    X_partitionword = computeInternalPredI(sampled_partitions, internal_tau_futures, i_set, param);
                else
                    try 
                        X_partitionword = X_info0.X_partitionword;
                    catch
                        X_partitionword = [];
                    end
                end


                if strcmp(data_sources{i_ds}(1:4), 'fish')
                    % 4) word-(p,v) stimulus information
                    % 5) perceptron-(p,v) stimulus information for final learned
                    % perceptrons, with lags
                    if numL0 <= 7
                        stim_tau_futures = -9:2:3;
                    else
                        stim_tau_futures = -9:2:3;
                    end
                    comp_type = 'CDM_pv';
                    if numL0 < 6
                        X_stim = computeBarmovieStimulusPredI(uniq_visited_percs, stim_tau_futures, i_set, comp_type, data_source);
                    else
                        % 'final_percs' is a subset of all visited percs,
                        % neglecting intermediate points. 
                        X_stim = computeBarmovieStimulusPredI(final_percs, stim_tau_futures, i_set, comp_type, data_source);
                    end

                else
                    X_stim = [];
                end
            
            
                save([dir_name file_tag], 'X_ww', 'X_pw', 'X_partitionword', 'X_stim');

                close all;
            else
                display([file_tag ' exists'])
            end
        end
    end

        

end