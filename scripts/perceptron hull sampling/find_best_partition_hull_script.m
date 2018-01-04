% script to find (1) optimal compression hull and (2) optimal perceptron
% hull for groups of 4 to 16


% this is how to compute pred-I using SEP's code
% [xbins, ~, data_x] = unique(input_wordnum(2:end));
% data_y = output_words(1:end-1) + 1;
% n_x = length(xbins);
% n_y = 2;
% 
% fracs =  [1.0000    0.9500    0.9000    0.8500    0.8000    0.7500    0.7000    0.5000];
% nreps = 50;
% 
% [I,I_err_std,I_err_frac,IR,IR_err_std,IR_err_frac,varargout] = ...
%       calc_info_P_joint(data_x,data_y,n_x,n_y,fracs,nreps);
% output_Ipred2 = [I I_err_std];


%% set up for 10 cells 
num_cells = 10;
sample_dir = 'sampled ten cell partitions/';
testset_prefix = 'test_setoddfishdata_set';
set_list = 1:20;


for i_set = 18:-1:1

X = load([sample_dir testset_prefix num2str(set_list(i_set))]);
fr_part = X.fr_proj;
pI_part = X.pI_proj;

fr_perc = X.fr_perc;
pI_perc = X.pI_perc;
spikingbin_times = X.spikingbin_times;
word_record = X.word_record;
full_word_record = zeros(num_cells, max(spikingbin_times) + 1);
full_word_record(:, spikingbin_times) = word_record;
tau_future = 1;
%% select 9 of 10 cells: would like to choose 9-cell group with highest input 
% efficiency

if ~exist(['g' num2str(num_cells) '_decimation_set' num2str(i_set) '.mat'], 'file')
cell_groups = true(num_cells);

group_Ipred = zeros(num_cells, 2);

[~, inputIpred, inputEntropy] = computePredictiveInfo(full_word_record', ...
                0*full_word_record(1, :)', tau_future);
group_Ipred(1, 1) =  inputIpred(1)/inputEntropy(1);
group_Ipred(1, 2) = inputIpred(1);

%%

figure()

for i_group = 2:(num_cells-1)

    inp_eff = zeros(num_cells, 2);

    for i_cell = 1:num_cells
        % if this cell is still in the partial group, compute the
        % predictive information when it is left out
        if cell_groups(i_group, i_cell)
            tau_future= 1;
            big_group_cells = find(cell_groups(i_group, :));
            cell_inds = setdiff(big_group_cells, i_cell);
            [~, inputIpred, inputEntropy] = computePredictiveInfo(full_word_record(cell_inds, :)', ...
                0*full_word_record(1, :)', tau_future);
            inp_eff(i_cell, 1) = inputIpred(1)/inputEntropy(1);
            inp_eff(i_cell, 2) = inputIpred(1);

        end
    end
    [~, dump_ind] = max(sum(inp_eff,2));
    
    cell_groups(i_group:end, dump_ind) = false;
    
    group_Ipred(i_group, 1) = inp_eff(dump_ind, 1);
    group_Ipred(i_group, 2) = inp_eff(dump_ind, 2);
    
    subplot(1, 3, [1 2])
    imagesc(cell_groups)
    
    subplot(1, 3, 3)
    imagesc(group_Ipred, [0 0.25])
    pause(0.01)
end
    
save(['g' num2str(num_cells) '_decimation_set' num2str(i_set)], 'cell_groups')

else
    load(['g' num2str(num_cells) '_decimation_set' num2str(i_set)], 'cell_groups')
end

%% now loop through all groups
clear effOut
effOut(num_cells-1, 1) = struct; 
%%
for i_group = 1:(num_cells-1)
    group_size = num_cells - i_group + 1;
    effOut(i_group).group_size = group_size;
    effOut(i_group).cell_inds = find(cell_groups(i_group, :));
    effOut(i_group).real_cellIDs = X.cell_set((cell_groups(i_group, :)));
    
    group_word_record = full_word_record(cell_groups(i_group, :), :);
    %% put inputs in number format
    conv_factor = 2.^(0:(group_size-1));
    word_nums = conv_factor*group_word_record;
    %% compute input efficiency
    tau_future= 1;
    [~, inputIpred, inputEntropy] = computePredictiveInfo(group_word_record', 0*group_word_record(1, :)', tau_future);
    input_efficiency = inputIpred(1)/inputEntropy(1);
    
    effOut(i_group).inputIpred = inputIpred;
    effOut(i_group).inputEntropy = inputEntropy;
    effOut(i_group).input_efficiency = input_efficiency;
    

    %% sample w-space
    % initialize - coarse
    num_samples = 1024;
    sample_vec_N = max(2, floor(num_samples^(1/group_size)));
    sample_vec = linspace(0, 1.02, sample_vec_N);
    b = 1;
    %%
    all_w = matrixNDGRID(sample_vec, false(1, group_size), nan(1, group_size));

    [all_Infos, uniqPartitions, uniq_fr, wisUniq] = computePerceptronPI_onMesh(all_w, b, group_word_record);
    old_w.all_w = all_w;
    old_w.all_info = all_Infos;
    old_w.uniqParts = uniqPartitions;
    old_w.uniwFR = uniq_fr;
    old_w.uniwPI = all_Infos(wisUniq == 1);
    old_w.uniW = all_w(wisUniq == 1, :);
    old_w.wisUniq = wisUniq;
    %% now decimate high-pI regions
    new_w = perceptronNeighborhood_infocalc(old_w, b, group_word_record);

    %% save
    effOut(i_group).init_w = old_w;
    effOut(i_group).samp_w = new_w;

    %% make some scatter plots
    makeMyFigure(30, 16)


    % first plot original scatter
    subplot(121)
    x = old_w(1).uniwFR;
    y = old_w(1).all_info(old_w.wisUniq==1)/inputIpred(1);
    plot(x, y, 'ko')

    x_fr = linspace(min(x), max(x), 100);

    % then plot each cluster
    hold on
    cluster_colors = lines(length(new_w));
    for i_k = 1:length(new_w)
        x = new_w(i_k).uniwFR;
        y = new_w(i_k).all_info(new_w(i_k).wisUniq == 1)/inputIpred(1);
        plot(x, y, '.', 'color', cluster_colors(i_k, :), 'markersize', 15);
    end


    hold off
    axis tight
    axis square
    xlabel('firing rate')
    ylabel('frac of Input pred-I')
    title(['groups of ' num2str(group_size)]);
    subplot(122)
    x = old_w(1).uniwFR;
    input_x_ent = input_efficiency*frEntropyBound(x).*x;
    y = old_w(1).all_info(old_w.wisUniq==1)./input_x_ent;
    plot(x, y, 'ko')
    hold on
    x_fr = linspace(min(x), max(x), 100);

    % then plot each cluster
    cluster_colors = lines(length(new_w));
    for i_k = 1:length(new_w)
        x = new_w(i_k).uniwFR;

        input_x_ent = input_efficiency*frEntropyBound(x).*x;

        y = new_w(i_k).all_info(new_w(i_k).wisUniq == 1)./input_x_ent;
        plot(x, y, '.', 'color', cluster_colors(i_k, :), 'markersize', 15);
    end

    plot(x_fr, 1 + 0*x_fr, 'k')
    hold off
    axis tight
    axis square
    xlabel('firing rate')
    ylabel('compression (perceptron efficiency/input efficiency)')

    pause(0.01)
end


save(['perceptron_grouprun_maxg' num2str(num_cells) '_set' num2str(i_set)], 'effOut');
close all


end
%% plotting parts
makeMyFigure(30, 30);

xlim_capacity = [0 0.1];
ylim_capacity = [0 1];
for i_group = 1:(length(effOut)-1)
    old_w = effOut(i_group).init_w;
    new_w = effOut(i_group).samp_w;
    group_size = effOut(i_group).group_size;
    input_efficiency = effOut(i_group).input_efficiency;
    
    if length(old_w.all_info) == length(old_w.uniwFR) && ~isfield(old_w, 'wisUniq')
        old_w.wisUniq = true(size(old_w.uniwFR));
    end
    
    % get the group input Ipred
    group_Ipredin = effOut(i_group).inputIpred(1);
        %% make some scatter plots
    subplot(3, 3, i_group)


    % first plot original scatter
    x = old_w(1).uniwFR;
    y = old_w(1).all_info(old_w.wisUniq==1)/group_Ipredin;
    plot(x, y, 'ko')

    x_fr = linspace(min(x), max(x), 100);

    % then plot each cluster
    hold on
    cluster_colors = lines(length(new_w));
    for i_k = 1:length(new_w)
        x = new_w(i_k).uniwFR;
        y = new_w(i_k).all_info(new_w(i_k).wisUniq == 1)/group_Ipredin;
        plot(x, y, '.', 'color', cluster_colors(i_k, :), 'markersize', 15);
    end


    hold off
    axis([xlim_capacity ylim_capacity])
    axis square
    xlabel('firing rate')
    ylabel('frac of Input pred-I')
    title(['groups of ' num2str(group_size)]);
    
end

%%

xlim_compression = [0 0.1];
ylim_compression = [0.5 2];
makeMyFigure(30, 30);
for i_group = 1:(length(effOut)-1)
    old_w = effOut(i_group).init_w;
    new_w = effOut(i_group).samp_w;
    group_size = effOut(i_group).group_size;
    input_efficiency = effOut(i_group).input_efficiency;
    
    if length(old_w.all_info) == length(old_w.uniwFR) && ~isfield(old_w, 'wisUniq')
        old_w.wisUniq = true(size(old_w.uniwFR));
    end
    
    
    subplot(3, 3, i_group)
    x = old_w(1).uniwFR;
    input_x_ent = input_efficiency*frEntropyBound(x).*x;
    y = old_w(1).all_info(old_w.wisUniq==1)./input_x_ent;
    plot(x, y, 'ko')
    hold on

    % then plot each cluster
    cluster_colors = lines(length(new_w));
    for i_k = 1:length(new_w)
        x = new_w(i_k).uniwFR;

        input_x_ent = input_efficiency*frEntropyBound(x).*x;

        y = new_w(i_k).all_info(new_w(i_k).wisUniq == 1)./input_x_ent;
        plot(x, y, '.', 'color', cluster_colors(i_k, :), 'markersize', 15);
    end

    plot(x_fr, 1 + 0*x_fr, 'k')
    hold off
    axis([xlim_compression ylim_compression])
    axis square
    xlabel('firing rate')
    ylabel('perceptron efficiency/input efficiency')

    pause(0.01)
end

%% plot max efficiency and max recovered pred-I

max_eff_group = ones(length(effOut)+1, 1);
max_rec_predI = ones(length(effOut)+1, 1);
size_of_group = ones(length(effOut)+1, 1);

% haven't computed these two yet
eff_at_maxrec = ones(length(effOut)+1, 1);
rec_at_maxeff = ones(length(effOut)+1, 1);


min_fr = 0.0167;
for i_group = 1:(length(effOut))
    size_of_group(i_group) = effOut(i_group).group_size;
    
    old_w = effOut(i_group).init_w;
    if length(old_w.all_info) == length(old_w.uniwFR) && ~isfield(old_w, 'wisUniq')
        old_w.wisUniq = true(size(old_w.uniwFR));
    end
    
    
    new_w = effOut(i_group).samp_w;
    
    % compute maximum pred-I
    uniq_infos = old_w.all_info(old_w.wisUniq == 1);
    mrg1 = max(uniq_infos(old_w.uniwFR > min_fr));
    %%
    for ii = 1:length(new_w)
        uniq_infos = new_w(ii).all_info(new_w(ii).wisUniq == 1);
        mrg1 = max(mrg1, max(uniq_infos(new_w(ii).uniwFR > min_fr)));
        
    end
    %%
    % normalize by input I pred
    max_rec_predI(i_group) = mrg1/effOut(i_group).inputIpred(1);
    
    % % compute maximum efficiency
    [~, frb] = frEntropyBound(old_w.uniwFR);
    uniq_info_effs = old_w.all_info(old_w.wisUniq == 1)./frb;
    meg1 = max(uniq_info_effs(old_w.uniwFR > min_fr));
    
    for ii = 1:length(new_w)
        [~, frb] = frEntropyBound(new_w(ii).uniwFR);
        uniq_info_effs = new_w(ii).all_info(new_w(ii).wisUniq == 1)./frb;
        meg1 = max(meg1, max(uniq_info_effs(new_w(ii).uniwFR > min_fr)));
        
    end
    % normalize by input I pred
    max_eff_group(i_group) = meg1/effOut(i_group).input_efficiency;
    
    
end

makeMyFigure(20, 10)
subplot(121)
plot(size_of_group, max_rec_predI, 'o-')
xlabel('N (# cells in group)')
ylabel({'maximum recovered pred-I';'(output pred-I/input pred-I)'})
axis square
axis([1 max(size_of_group) 0 1])
set(gca, 'fontsize', 14)

subplot(122)
plot(size_of_group, max_eff_group, 'o-')
xlabel('N (# cells in group)')
ylabel(['maximum compression w/ fr > ' num2str(min_fr*60, '%1.2f') 'Hz'])
axis square 
axis([1 max(size_of_group) 0.5 2])
set(gca, 'fontsize', 14)

% %% fminunc method: does not work on perceptron weights because info gradient 
% % is zero inside any partition 
% % COULD work on probabilities. Need to calculate gradient. 
% b = 1;
% word_record = X.word_record;
% conv_factor = 2.^(0:(num_cells-1));
% word_nums = conv_factor*word_record;
% nz_word_inds = unique(word_nums) + 1;
% 
% fr_lambda = 1;
% pIfrc = @(x) probrulePredI_fr_objfun(x, nz_word_inds, full_word_record, fr_lambda) ;
% 
% p0 = rand(length(nz_word_inds), 1);



% %% set initial condition for optimizing perceptrons
% fr_min = 0.0;
% fr_max = 0.03;
% fr_range = find(fr_perc > fr_min & fr_perc < fr_max);
% 
% [~, ii] = max(pI_perc(fr_range));
% 
% wb_i = X.allowed_wb(fr_range(ii), :);
% b = 1;
% 
% %% compute input pred-I
% 
% [xbins, ~, data_xy] = unique(word_nums(1:end));
% data_x = data_xy(2:end);
% data_y = data_xy(1:(end-1));
% 
% n_x = length(xbins);
% n_y = n_x;
% 
% fracs =  [1.0000    0.950   0.8500    0.7500     0.5000];
% nreps = 50;
% 
% [I,I_err_std,I_err_frac,IR,IR_err_std,IR_err_frac,varargout] = ...
%       calc_info_P_joint(data_x,data_y,n_x,n_y,fracs,nreps);
% input_Ipred = [I I_err_std];
% 
% %% gradient ascent on perceptron weights
% dw = 0.2;
% num_iters = 100;
% wb_0 = wb_i;
% wb_track = zeros(num_cells, num_iters);
% info_track = zeros(1, num_iters);
% for i_iter = 1:num_iters
%     wb_track(:, i_iter) = wb_0;
%     tic
%     [I0, dI] = computeInfoGradient_Perceptron(wb_0, b, dw, full_word_record); 
%     toc
%     wb_0 = wb_0 + dI'/I0;
%     info_track(i_iter) = I0;
% end
