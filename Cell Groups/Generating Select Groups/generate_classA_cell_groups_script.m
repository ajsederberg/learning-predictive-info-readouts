% Script to generate groups in the "All" class

num_groups = 1000;

%% For the multi-movie dataset

% multimovie data exploration : which cells are interesting? which movies
% are interesting? where to start with generating training sets? 

X = load('2011-11-03_multmovie_binned.mat', 'binned');
mmdata = X.binned;

[numReps, numFrames, numCells, numMovies] = size(mmdata);

% look at movie responses for a single cell 

trials_per_movie = [83 80 84 91 85]';   % this is calculated below
div1 = round(trials_per_movie/3);
div2 = round(2*trials_per_movie/3);

max_fr_change = 0.6;    % relative firing rate change permitted over course of experiment 
min_mov_fr = (1/20)/60; % require at least one spike per trial on average
min_cc = zeros(numMovies, numCells);
median_cc = zeros(numMovies, numCells);
active_count = zeros(numMovies, numCells);
last_active = zeros(numMovies, numCells);
first_active = last_active;
trend_test = zeros(numMovies, numCells);
rel_fr_change = zeros(numMovies, numCells);
movie_fr = zeros(numMovies, numCells);

trial_criterion = zeros(numMovies, numCells);
plotit = 1;
for i_cell = 1:numCells
    

    cc_midtrial = zeros(numMovies, numReps);
    if plotit
        h = makeMyFigure(40, 15);
    end
    for i_mov = 1:numMovies
        cellraster = squeeze(mmdata(1:trials_per_movie(i_mov), :, i_cell, i_mov));
        active_trials = find(sum(cellraster, 2) > 0);
        active_count(i_mov, i_cell) = length(active_trials);
        lat = active_trials(end); 
        fat = active_trials(1);
        last_active(i_mov, i_cell) = lat;
        first_active(i_mov, i_cell) = fat;
        
        if mod(length(active_trials),2) == 0
            middle_trial = active_trials((end/2));
        else
            middle_trial = median(active_trials);
        end
        % test for stationarity: is firing rate of the first 1/3 of trials 
        % indistinguishable from the firing rate of the last 1/3 of trials?
        trial_fr = mean(cellraster, 2);
        movie_fr(i_mov, i_cell) = mean(trial_fr);
        trend_test(i_mov, i_cell) = ranksum(trial_fr(1:div1(i_mov)), ...
            trial_fr(div2(i_mov):trials_per_movie(i_mov)));
        pp = polyfit((1:lat)'/lat, trial_fr(1:lat), 1);
        rel_fr_change(i_mov,i_cell) = pp(1)/max(trial_fr(1:lat));  
        
        if plotit
            subplot(2, 6, i_mov)
            imagesc(cellraster, [0 1])
            xyoff
            set(gca, 'fontsize', 12)

            title({[num2str(length(active_trials)) ...
                ' active trials']; ['last active trial ' num2str(lat)]; ...
                ['cell ' num2str(i_cell) ', movie ' num2str(i_mov)]})

            subplot(2, 6, 6 + i_mov)
            plot(60*trial_fr)
            xlabel('trial')
            ylabel('mean fr (hz)')
            axis([0 lat 0 10])
        
        end
        xc = corrcoef(cellraster');
        cc_midtrial(i_mov, 1:trials_per_movie(i_mov)) = xc(middle_trial, :);
        
        min_cc(i_mov, i_cell) = min(cc_midtrial(i_mov, active_trials));
        median_cc(i_mov, i_cell) = median(cc_midtrial(i_mov, active_trials));
    end
    
    % criterion is that the last active trial was past the 90% point
    % thus, low-fr cells are not discarded solely on the basis of firing
    % rate
    trial_criterion(:, i_cell) = ...
        first_active(:, i_cell) < 0.1*trials_per_movie & ...
        last_active(:, i_cell) > 0.9*trials_per_movie & ...
        abs(rel_fr_change(:, i_cell)) < max_fr_change & ...
        movie_fr(:, i_cell) > min_mov_fr;
    
    if plotit
        subplot(2, 6, 6)
        imagesc(cc_midtrial, [0 .8])
        colorbar
        set(gca, 'fontsize', 12)
        if all(trial_criterion(:, i_cell), 1)
            title({'CELL KEPT'; 'corrcoef with middle trial'})
        else
            title({'CELL NOT KEPT'; 'corrcoef with middle trial'})
        end

        print(h, '-dpdf', ['Cell Groups/multimovie cell snapshots/multimovie_cell' num2str(i_cell)])
        pause(.01)
        close(h)
    end
end

% keep cells that passed on 4 or more movies
     % small relative change in fr
     % first responsive trial in first 10% of trials
     % last responsive trial in last 10% of trials
cells_to_use = find(sum(trial_criterion, 1) >= 4);

trials_per_movie = max(last_active, [], 2); 
%%
save('Cell Groups/multimovie cell snapshots/cell_list_20111103', 'cells_to_use');
%% 
makeMyFigure(20, 30)
subplot(4, 1, 1)
imagesc(min_cc(:, cells_to_use), [-.1 .1])
colorbar
title('minimum active trial cross corr')
ylabel('movie')

subplot(4, 1, 2)
imagesc(median_cc(:, cells_to_use), [-.1 0.5])
colorbar
title('median active trial cross corr')
ylabel('movie')

subplot(4, 1, 3)
imagesc(last_active(:, cells_to_use))
colorbar
title('last active trial')
ylabel('movie')

subplot(4, 1, 4)
imagesc(active_count(:, cells_to_use))
colorbar
title('activity tally (#trials per cell per movie)')
xlabel('cell')
ylabel('movie')

%% generate groups of 3 to groups of 10


% draw cell groups, look at average cell_cc
for num_cells_in_group = 3:10
    num_kept_cells = length(cells_to_use);
    cell_sets = zeros(num_groups, num_cells_in_group);
    % for checking uniqueness
    conv_factor = num_kept_cells.^(0:(num_cells_in_group-1));
    group_numbers = conv_factor*cell_sets';

    while length(unique(group_numbers)) < num_groups
        n_uniq = length(unique(group_numbers)); % keep the unique ones
        
        % loop up to num_groups to search for more unique groups
        for i_group = n_uniq:num_groups
            cell_sets(i_group, :) = sort(randperm(num_kept_cells, num_cells_in_group));
        end

        group_numbers = conv_factor*cell_sets';
        
        [group_numbers, i_u] = unique(group_numbers);
        cell_sets(1:length(i_u), :) = cell_sets(i_u, :);

    end
    sets_of_N = cells_to_use(cell_sets);

    
    save(['Cell Groups/multimovieA_setsOf' num2str(num_cells_in_group)], 'sets_of_N');
end

%% For fishmovie: draw from the full set of cells, generate groups of 3 to groups of 10

cells_to_use_fm = 1:53;
% draw cell groups, look at average cell_cc
for num_cells_in_group = 3:10
    num_kept_cells = length(cells_to_use_fm);
    cell_sets = zeros(num_groups, num_cells_in_group);
    
    % this is a special example set from Fig 1
    if num_cells_in_group == 4 
        cell_sets(1, :) = [13 16 27 49];
    end
    
    % for checking uniqueness
    conv_factor = num_kept_cells.^(0:(num_cells_in_group-1));
    group_numbers = conv_factor*cell_sets';

    while length(unique(group_numbers)) < num_groups
        n_uniq = length(unique(group_numbers)); % keep the unique ones
        % loop up to num_groups to search for more unique groups
        for i_group = n_uniq:num_groups
            cell_sets(i_group, :) = sort(randperm(num_kept_cells, num_cells_in_group));
        end

        group_numbers = conv_factor*cell_sets';
        
        [group_numbers, i_u] = unique(group_numbers);
        cell_sets(1:length(i_u), :) = cell_sets(i_u, :);

    end
    sets_of_N = cells_to_use_fm(cell_sets);

    
    save(['Cell Groups/fishmovieA_setsOf' num2str(num_cells_in_group)], 'sets_of_N');
end