function script2_run_learning_simulations_script(numL0_list, i_set_list, learning_types, set_w_grid, data_sources)
%Script to run all learning simulations on a single set of cells. 
% Written as a function to run for different set sizes, but can easily go
% back to being a script if numL0_list is set. numL0_list is a subset of
% [3 4 5 6 7 8 9 10]. 
% learning_types = {'vanilla', 'vanilloh', 'triplet', 'triploh'};

if nargin < 5
    data_sources = {'fishmovieS', 'fishmovieA'}; %, 'multimovieS', 'multimovieA'};
end

for numL0 = numL0_list
    for i_set = i_set_list
        for i_lt = 1:length(learning_types)
            learning_type = learning_types{i_lt};
            
            parfor i_ds = 1:length(data_sources)
                
                data_source = data_sources{i_ds};
                
                if set_w_grid
                    perc_start = 0.55*ones(1, numL0);
                    b = 1;
                    num_ics = 2000;
                    max_dw = 0.5;
                    [init_ws, init_ws_Orule] = sampleWSpaceForRule(perc_start, b, num_ics, max_dw);
                    param = setLearningParameters(learning_type, numL0, ...
                        data_source, 'w0_uniform', false, 'w0_values', [init_ws; init_ws_Orule], ...
                        'training_reps', 100, 'saveandprint', true);

                else              
                    param = setLearningParameters(learning_type, numL0, ...
                        data_source, 'saveandprint', true);
                end                
                % get the name of the file that will run
                param_fs = param;
                param_fs.return_file_string_only = true;
                file_string = binary_neuron_learning_simulation(i_set, param_fs);
                file_name = ['learning simulation results/' data_source '_' learning_type '_' file_string '.mat'];

                % if that file doesn't already exist, run the simulation
                if ~exist(file_name, 'file')
                    try binary_neuron_learning_simulation(i_set, param);
                    catch 
                        display([file_name ' failed to complete'])
                    end
                end
            end
        end
        
        close all;
        
        % close the pool to prevent hang-ups (may not be a problem in
        % 2016a, was a problem in 2014b and before.
        delete(gcp('nocreate'))
    end
end