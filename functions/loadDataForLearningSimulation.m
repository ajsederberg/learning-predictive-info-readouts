function [train_s, test_s, i_set_string, warning_string] = loadDataForLearningSimulation(i_set, param)
% data_source options: 'fishmovieA', 'fishmovieS', 'multimovieA', 'multimovieS' 
% Load training data. Can also be used to load the non-rep bar data, but
% this will not segment non-rep bar data into test/training sets. If
% param.movie_type is 'nrb', the outputs are: 
% [nrb_stim, nrb_data, i_set_string, warning_string]

warning_string = '';

if strcmp(param.movie_type, 'fish')
    %% load the data
    if strcmp(param.data_source, 'fishmovieA') || strcmp(param.data_source, 'fishmovieS')
        X = load('2011-02-25_barfishcheck_experiment_and_calcs_PNAS_paper.mat', ...
            'fishdata', 'dt');
        fishdata = X.fishdata;
        binsize = X.dt;
        clear X
    else
        X = load('2011-11-03_multmovie_binned', 'binned');
        fishdata = squeeze(X.binned(:, :, :, 4));    % fourth movie is fish - need to confirm
            % X.binned is [numReps, numFrames, numCells]
            % fishdata is [numCells, numReps, numFrames]
        fishdata = permute(fishdata, [3 1 2]);
    end
    % load the new groups    
    fG = load(['Cell Groups/' param.data_source '_setsOf' num2str(param.numL0)]);
    
    cell_set = fG.sets_of_N(i_set, :);

    s = fishdata(cell_set, :, :);
    numL0 = length(cell_set);    

    string_tag = param.data_source([1 end]);
    i_set_string = [num2str(i_set) '_' string_tag num2str(numL0)];    % test set (odd trials)
    % take odd trials for test set
    test_s = permute(s(:, 1:2:end, :), [3 2 1]); 
    test_s = reshape(test_s, numel(test_s)/length(cell_set), length(cell_set))';

    % in training set (even trials), will take NT trials at a time and repeat
    % the entire set of trials training_reps number of times
    if param.training_reps == 1
        NT = 1;
        training_reps = param.training_reps;
    else
        training_reps = 4*round(param.training_reps/4);
        NT = training_reps/4;
    end
    
    train_s = permute(s(:, 2:2:end, :), [3 2 1]);
    train_s = repmat(train_s, [1 training_reps 1]);
    
    % present individual 1141-frame "trials" in a random order - this
    % preserves most structure (1140/1141 frames) but doesn't repeat trials
    % in the same order every time
    trial_ord = randperm(size(train_s, 2));
    train_s = train_s(:, trial_ord, :);
    
    [m,n,p] = size(train_s);
 
    
    if rem(n, NT) ~= 0
        error(['NT must be a factor of ' num2str(n)]) 
    end

    train_s = reshape(train_s, [NT*m n/NT p]);
    % this puts train_s in (cells)x(time points)x(reps)
    train_s = permute(train_s, [3 1 2]);
    


elseif strcmp(param.movie_type, 'nrb') && strcmp(param.data_source(1:4), 'fish')
    % only fishmovieA/fishmovieS have nrb data
    
    % load the new groups    
    fG = load(['Cell Groups/' param.data_source '_setsOf' num2str(param.numL0)]);
    
    cell_set = fG.sets_of_N(i_set, :);
    try comp_type = param.comp_type;
    catch
        comp_type = 'CDM_pv';  
    end
    
    [nrb_stim, nrb_data] = barmovie_activity_readout_function(learned_ws, cell_set, comp_type);
    test_s = nrb_data;
    train_s = nrb_stim;
    binsize = 'WARNING THIS IS NOT IN TEST AND TRAIN FORMAT';
    
    string_tag = param.data_source([1 end]);
    i_set_string = [num2str(i_set) '_' string_tag num2str(param.numL0)];    % test set (odd trials)
    
elseif strcmp(param.data_source(1:end-1), 'checkerboard')   
    
    %%
    X = load('2011-02-25_barfishcheck_experiment_and_calcs_PNAS_paper.mat', ...
        'checkerboarddata', 'dt');
    checkerboarddata = X.checkerboarddata;
    binsize = X.dt;
    clear X
    %%
        
    % load the new groups    
    fG = load(['Cell Groups/fishmovie' param.data_source(end) '_setsOf' num2str(param.numL0)]);
    
    cell_set = fG.sets_of_N(i_set, :);

    s = checkerboarddata(cell_set, :);
    numL0 = length(cell_set);    

    string_tag = ['f' param.data_source(end)];
    i_set_string = [num2str(i_set) '_' string_tag num2str(numL0)];    % test set (odd trials)
    
    %% basically copy the fishdata structure, but let trials be about twice as long
    trial_length = 2285;
    num_trials = floor(size(s, 2)/trial_length);
    s_trials = reshape(s(:, 1:(trial_length*num_trials)), numL0, trial_length, num_trials);
    s_trials = permute(s_trials, [2 3 1]);
    %%
     % take odd trials for test set
    test_s = s_trials(:, 1:2:end, :); 
    test_s = reshape(test_s, numel(test_s)/length(cell_set), length(cell_set))';

    % in training set (even trials), will take NT trials at a time and repeat
    % the entire set of trials training_reps number of times
    if param.training_reps == 1
        NT = 1;
        training_reps = param.training_reps;
    else
        training_reps = 4*round(param.training_reps/4);
        NT = training_reps/4;
    end
    
    train_s = s_trials(:, 2:2:end, :);
    train_s = repmat(train_s, [1 training_reps 1]);
    
    % present individual 1141-frame "trials" in a random order - this
    % preserves most structure (1140/1141 frames) but doesn't repeat trials
    % in the same order every time
    trial_ord = randperm(size(train_s, 2));
    train_s = train_s(:, trial_ord, :);
    
    [m,n,p] = size(train_s);
 
    
    if rem(n, NT) ~= 0
        error(['NT must be a factor of ' num2str(n)]) 
    end

    train_s = reshape(train_s, [NT*m n/NT p]);
    % this puts train_s in (cells)x(time points)x(reps)
    train_s = permute(train_s, [3 1 2]);
end