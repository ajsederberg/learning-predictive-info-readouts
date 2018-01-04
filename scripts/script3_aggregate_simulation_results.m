function script3_aggregate_simulation_results(set_size, generate_partitions, data_sources)
%% Script to process the simulation results and save in the aggregated folder. 

if ~generate_partitions
% Creating one file for each data source of each set size. 

if nargin < 3
    data_sources = {'fishmovieS', 'fishmovieA', 'multimovieS', 'multimovieA'};
end

for i_ds = 1:length(data_sources)
    
    result_list = dir('learning simulation results');
    
    group_tag = [data_sources{i_ds}([1 end]) num2str(set_size)];

    
    % look for the sets that match the datasource
    result_matches_datasource = arrayfun(...
        @(x) ~isempty(regexp(x.name, [data_sources{i_ds} '_\w*.mat'], 'match')), ...
        result_list);

    
    datasource_file_list = result_list(result_matches_datasource);
    % create an empty cell array to hold results
    % Use cellfun in analyses to search for results with particular
    % attributes. 
    res_ds_arr = cell(length(datasource_file_list), 1);    
    for i_dfl = 1:length(datasource_file_list)
        X = load(['learning simulation results/' datasource_file_list(i_dfl).name]);
        
        set_ID = regexp(datasource_file_list(i_dfl).name, 'set(\w*)_', 'tokens');
        X.set_ID = str2double(set_ID{1});
        res_ds_arr{i_dfl} = X;
        
    end
    save(['aggregated results/' group_tag '_learning_results'], 'res_ds_arr', '-v7.3')
    
    

    
    all_learned_perceptrons_arr = cellfun(@(x) x.trialW, res_ds_arr, ...
        'uniformoutput', false);
    
    


    %% add these to the list of known perceptrons. 
    perceptron_list_file = ['Cell Groups/perceptron_list_size' num2str(set_size) '.mat'];
    if exist(perceptron_list_file, 'file')
        Y = load(perceptron_list_file);

        w_perc = Y.unique_perceptrons;
        part_perc = Y.unique_partitions;
        
        w_percFC = Y.final_perceptrons;     % for final points of learning simuluation
        part_percFC = Y.final_partitions;   

    
    else
        w_perc = [];
        part_perc = [];
        
        w_percFC = [];
        part_percFC = [];
    end
%%  This is a list of all visited perceptrons
    
    all_learned_perceptrons = cellfun( ...
        @(x) reshape(x, size(x,1)*size(x,2), size(x, 3)), ...
        all_learned_perceptrons_arr, 'uniformoutput', false);

    
    
    all_learned_partitions = cellfun(@(x) perceptronToProjectionRule(x, 1), ...
        all_learned_perceptrons, 'uniformoutput', false);

    %% This is a separate list of the final perceptrons only
    all_final_perceptrons = cellfun( ...
        @(x) squeeze(x(end, :, :)), ...
        all_learned_perceptrons_arr, 'uniformoutput', false);
    all_final_partitions = cellfun(@(x) perceptronToProjectionRule(x, 1), ...
        all_final_perceptrons, 'uniformoutput', false);
%% find the list of all unique perceptron/partitions
    
    [unique_learned_partitions, uniq_inds] = (cellfun(@(x) returnUniquePartitionSet_RPMethod(x), ...
        all_learned_partitions, 'uniformoutput', false));

    % append the already-known perceptrons to the unique learned perceptrons
    redundant_list_unique_perceptrons = [cell2mat(cellfun(@(x,y) x(y, :), all_learned_perceptrons, uniq_inds, ...
        'uniformoutput', false)); w_perc];
       
    % append the already-known partitions to the unique learned partitions
    % then flatten the cell array and find the unique set of perceptrons
    [unique_partitions, uber_uniq_inds] = returnUniquePartitionSet_RPMethod(...
        [cell2mat(unique_learned_partitions); part_perc]);

    unique_perceptrons = redundant_list_unique_perceptrons(uber_uniq_inds, :);

%% find the list of all unique FINAL (end of learning) perceptron/partitions    
    [unique_final_partitions, uniq_final_inds] = (cellfun(@(x) returnUniquePartitionSet_RPMethod(x), ...
        all_final_partitions, 'uniformoutput', false));

    % append the already-known perceptrons to the unique learned perceptrons
    redundant_list_final_perceptrons = [cell2mat(cellfun(@(x,y) x(y, :), all_final_perceptrons, uniq_final_inds, ...
        'uniformoutput', false)); w_percFC];
       
    % append the already-known partitions to the unique learned partitions
    % then flatten the cell array and find the unique set of perceptrons
    [final_partitions, uber_uniq_final_inds] = returnUniquePartitionSet_RPMethod(...
        [cell2mat(unique_final_partitions); part_percFC]);

    final_perceptrons = redundant_list_final_perceptrons(uber_uniq_final_inds, :);
    
    %%
    save(perceptron_list_file, 'unique_partitions', 'unique_perceptrons', ...
        'final_partitions', 'final_perceptrons')

end

%% generate sampled partition lists
else 

    partition_list_file = ['Cell Groups/partition_list_size' num2str(set_size) '.mat'];

    % 
    num_partitions_sampled = 5000;
    flip_fractions = [0.05 0.1 0.15 0.2 0.25];
    partitions_per_flip = num_partitions_sampled/length(flip_fractions);
    X_up = load(['Cell Groups/perceptron_list_size' num2str(set_size) '.mat']);

    perceptron_seeds = X_up.final_partitions;

    sampled_partitions = cell(length(flip_fractions), 1);
    starter_partitions = cell(length(flip_fractions), 1);
    for i_flip = 1:length(flip_fractions)
        flip_set_inds = ceil(size(perceptron_seeds, 1)*rand(partitions_per_flip, 1));

        flip_seed = perceptron_seeds(flip_set_inds, :);

        flipped_entries = rand(size(flip_seed)) <= flip_fractions(i_flip);

        flip_seed(flipped_entries) = 1 - flip_seed(flipped_entries);

        % keep zero-word response 0
        flip_seed(:, 1) = 0;
        sampled_partitions{i_flip} = flip_seed;
        starter_partitions{i_flip} = perceptron_seeds(flip_set_inds, :);
    end

    fraction_flipped = cell2mat(cellfun(@(x, y) y*ones(size(x, 1),1), sampled_partitions, ...
        num2cell(flip_fractions'), 'uniformoutput', false));
    sampled_partitions = cell2mat(sampled_partitions);
    sampled_partition_seeds = cell2mat(starter_partitions);

    save(partition_list_file, 'sampled_partitions', 'fraction_flipped', ...
            'sampled_partition_seeds')
end