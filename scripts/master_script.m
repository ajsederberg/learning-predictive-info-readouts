%% The general workflow is as follows:

% (script1) : preliminaries. Define the cell groups (already done, see the
% scripts in Cell Groups/Generating Select Groups/ for details.). This is
% done once, before the remaind of analyses. 

% script2 : run the learning simulations
% script3 : aggregate the simulation results and save in a single file
% script4 : these five scripts are used to compute predictive information. 
%     % if all calculations are run at once, script4 will do everything, 
%         but take a very long time. 
%              % script4b is a simplified version that only computes word
%                         information values (not readout information)
%              % script4c is the same as script4, but it only computes
%              information for perceptron readouts
%              % script4d computes "token" information (i.e. information
%              carried by a particular word)
%              % script4e computes word-word info for the checkerboard
%              stimulus 
        % Info Calculations required:
        % computeInternalPredI function:
        %  1) Word-word internal pred I, with lags 1:10
        %  2) Perceptron-word internal pred I, for all perceptrons visited by
        %  learning (1 lag only)
        %  3) Partition-word internal pred I, for sampled partitions
        % computeBarmovieStimulusPredI: 
        %  4) Word-(p,v) stimulus information, with lags -10:1:10 [for fA and fS
        %  sets only)
        %  5) Perceptron-(p,v) stimulus information, for final learned perceptrons,
        %  with lags -9:2:3. [for fA and fS sets only]

% script5 : aggregate the information calculation results and save in a
%               single file

%% run learning simulations. this is broken into chunks. Change "i_set_start" and i_set_list to run
% for smaller subsets. 

% learning parameters are set in the function setLearningParameters, which
% defines default parameters (which can be modified by (Name, Value) input
% pairs). This is called from within "script2"

% run many more learning simulations, with more ICs for basin analysis
learning_types = {'vanilla', 'vanilloh', 'triplet', 'triploh'};
set_w_grid = true;      % if true, this runs for 2000 IC's in a regular grid; this number can be adjusted inside script2
                        % if false, this runs for 200 randomly selected
                        % initial conditions with weights between w0min and
                        % w0max. 
dir_name = 'information calculation results/datetime_20161221/';

for numL0 = 4:3:10
    for i_set_start = 1:20:980
        i_set_list = i_set_start + (0:1:16);
        try script2_run_learning_simulations_script(numL0, i_set_list, learning_types, set_w_grid)
        catch
            display(['something broke for ' num2str(numL0) ', set ' num2str(i_set_list)]);

        end
    
    end

    % aggregate the new information
    try script3_aggregate_simulation_results(numL0, false)
    catch
        display(['aggregation failed for sets of ' num2str(numL0)])
    end
    
    % compute the new information
    for i_set_start = 1:20:180
        i_set_list = i_set_start + (0:1:18);
        try script4_cellset_predictiveinfo_script(numL0, i_set_list, true, dir_name)
        catch
            display(['info failure for ' num2str(numL0) ', set ' num2str(i_set_list)]);
        end
    end

    
end
%% run more pred-I calculations 

% for word info 
for numL0 = 4:3:10

    try script4b_cellset_predictiveinfo_script(numL0, 1:1000)

    catch 
        disp(['oops - fail for sets of ' num2str(numL0)])
    end


end





%% create a file with the frequency counts (rate) for each word, in each group

data_sources = {'fishmovieA', 'fishmovieS'};% , 'multimovieA', 'multimovieS'};

for i_ds = 1:length(data_sources)
    for numL0 = 4:3:10
        param = setLearningParameters('', numL0, ...
                        data_sources{i_ds});

        word_probs = cell(1000, 1);

        parfor i_set = 1:1000
            [~, test_s] = loadDataForLearningSimulation(i_set, param);
            
            [nn, ocnts, words] = words2nnOcnts(test_s');
            
            word_probs{i_set}.words = words;
            word_probs{i_set}.probs = nn/sum(nn);
            
        end
        
        file_tag = [data_sources{i_ds}([1 end]) num2str(numL0)];
        save(['aggregated results/word_probabilites_' file_tag], 'word_probs')
    end
end


%%




%% Running info calcs for missing sets (with learning results)
% This will check for required files based on what has learning results run
% already. If the info calculations were already run, then they can be
% loaded from the 'efficient_g10' folder using this code. This also checks
% that any info calculations are performed for all learned perceptrons
% (even retroactively, for perceptrons "discovered" after the pred-I was
% calculated for an early simulation.)


%%
script4c_cellset_predictiveinfo_script(10, ...
    'information calculation results/final_g10/', ...
    {'information calculation results/old_final_g10/', ...
    'information calculation results/efficient_g10/', ...
    'information calculation results/datetime_20170121/', ...
    'information calculation results/datetime_20170122/'});
%%
script4c_cellset_predictiveinfo_script(7, ...
    'information calculation results/final_g7/', ...
    {'information calculation results/efficient_g7/', ...
    'information calculation results/datetime_20170115/'});
