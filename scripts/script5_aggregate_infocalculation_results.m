%% Script to process the simulation results and save in the aggregated folder. 

% Creating one file for each data source of each set size. 
set_size = 10;
data_sources = {'fishmovieS', 'fishmovieA'};
info_subfolder = 'information calculation results/final_g10/';

%%
for i_ds = 1:length(data_sources)
    %%
    result_list = dir(info_subfolder);
    group_tag = [data_sources{i_ds}([1 end]) num2str(set_size)];

    % look for the sets that match the datasource
    result_matches_datasource = arrayfun(...
        @(x) ~isempty(regexp(x.name, ['set\w*_' group_tag '.mat'], 'match')), ...
        result_list);
%%
    perc_lookup_struct = load([info_subfolder 'sampledPercs_setsOf' num2str(set_size) '.mat']);
    datasource_file_list = result_list(result_matches_datasource);
    % create an empty cell array to hold results
    % Use cellfun in analyses to search for results with particular
    % attributes. 
    info_ds_arr = cell(length(datasource_file_list), 1);    
    for i_dfl = 1:length(datasource_file_list)
        X = load([info_subfolder datasource_file_list(i_dfl).name]);
        
        set_ID = regexp(datasource_file_list(i_dfl).name, 'set(\w*)_', 'tokens');
        X.set_ID = str2double(set_ID{1});
        info_ds_arr{i_dfl} = X;
        
    end
    save(['aggregated results/' group_tag '_info_calcs'], 'info_ds_arr', ...
        'perc_lookup_struct', '-v7.3')
    
    
    


end

%% Create one file with all information calculations for a set of simulation results

% Creating one file for each data source of each set size. 
set_size = 7;
data_sources = {'fishmovieS', 'fishmovieA'};
info_directory = 'information calculation results/';
info_subdirs = {'january17_g7/', 'datetime_20161220/'};

%%
for i_ds = 1:length(data_sources)
    %% load sim results
    group_tag = [data_sources{i_ds}([1 end]) num2str(set_size)];

    X_res = load(['aggregated results/' group_tag '_learning_results.mat']);
    
    set_IDs = unique(cellfun(@(x) x.i_set, X_res.res_ds_arr));
    
    % look for info calculations for each result file. 
    info_res_arr = cell(length(set_IDs), 1);    
    for i_dfl = 1:length(set_IDs)
        
        i_sd = 1;
        while i_sd <= length(info_subdirs)
            file_name = [info_directory info_subdirs{i_sd} 'set' ...
                num2str(set_IDs(i_dfl)) '_' group_tag '.mat'];
            try X = load(file_name);
        
                X.set_ID = set_IDs(i_dfl);
                X.perc_lookup_struct_file = [info_subdirs{i_sd} 'sampledPercParts_setsOf' num2str(set_size) '.mat'];
                info_res_arr{i_dfl} = X;
                i_sd = length(info_subdirs) + 1;
            catch
                i_sd = i_sd + 1;
            end
        end
        
        if isempty(info_res_arr{i_dfl})
            info_res_arr{i_dfl} = 'no info calcs available';
        end
    end
    save(['aggregated results/' group_tag '_info_calcs_for_sims'], 'info_res_arr', ...
        '-v7.3')
%     
    
    


end