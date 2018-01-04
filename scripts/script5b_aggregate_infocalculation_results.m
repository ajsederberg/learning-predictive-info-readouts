%% Script to process the simulation results and save in the aggregated folder for word-word info. 

% Creating one file for each data source of each set size. 
set_size = 10;
data_sources = {'fishmovieS', 'fishmovieA'};
info_subfolder = 'information calculation results/token_word_g10/';

for i_ds = 1:length(data_sources)
    %%
    result_list = dir(info_subfolder);
    group_tag = [data_sources{i_ds}([1 end]) num2str(set_size)];

    % look for the sets that match the datasource
%     result_matches_datasource = arrayfun(...
%         @(x) ~isempty(regexp(x.name, ['wordword_set\w*_' group_tag '.mat'], 'match')), ...
%         result_list);
    result_matches_datasource = arrayfun(...
        @(x) ~isempty(regexp(x.name, ['tokeninfo_set\w*_' group_tag '.mat'], 'match')), ...
        result_list);
    %%
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
    save(['aggregated results/tokenword_' group_tag '_info_calcs'], 'info_ds_arr', ...
        '-v7.3')
    
    
    


end



%% Creating one file for each data source of each set size. 
set_size = 10;
data_sources = {'fishmovieS', 'fishmovieA'};
info_subfolder = 'information calculation results/wordword_sets/';

for i_ds = 1:length(data_sources)
    %%
    result_list = dir(info_subfolder);
    group_tag = [data_sources{i_ds}([1 end]) num2str(set_size)];

    % look for the sets that match the datasource
%     result_matches_datasource = arrayfun(...
%         @(x) ~isempty(regexp(x.name, ['wordword_set\w*_' group_tag '.mat'], 'match')), ...
%         result_list);
    result_matches_datasource = arrayfun(...
        @(x) ~isempty(regexp(x.name, ['wordword_set\w*_' group_tag '.mat'], 'match')), ...
        result_list);
    %%
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
    save(['aggregated results/wordword_' group_tag '_info_calcs'], 'info_ds_arr', ...
        '-v7.3')
    
    
    


end