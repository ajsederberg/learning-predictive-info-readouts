% working script while 'find_best_partition_hull_script' is computing

%% set up for 10 cells 
num_cells = 10;
sample_dir = 'sampled ten cell partitions/';
testset_prefix = 'test_setoddfishdata_set';
set_list = 1:20;
%%
for i_set = 18:18

X = load([sample_dir testset_prefix num2str(set_list(i_set))]);
fr_part = X.fr_proj;
pI_part = X.pI_proj;

fr_perc = X.fr_perc;
pI_perc = X.pI_perc;
spikingbin_times = X.spikingbin_times;
word_record = X.word_record;
full_word_record = zeros(num_cells, max(spikingbin_times) + 1);
full_word_record(:, spikingbin_times) = word_record;

%% select 9 of 10 cells: would like to choose 9-cell group with highest input 
% efficiency
tau_future = 1;

x_cg = load(['g' num2str(num_cells) '_decimation_set' num2str(i_set)], 'cell_groups');
cell_groups = x_cg.cell_groups;

%% now loop through all groups


b = 1;
tau_future= 1;
fr_lambda = 200;

X_per = load(['perceptron_grouprun_maxg' num2str(num_cells) '_set' num2str(i_set)], 'effOut');

%%
clear probRuleOut
probRuleOut(num_cells) = struct();
%%
for i_group = 1:(num_cells-1)
    group_size = num_cells - i_group + 1;
    cell_inds = find(cell_groups(i_group, :));
    real_cellIDs = X.cell_set((cell_groups(i_group, :)));
    
    probRuleOut(i_group).group_size = group_size;
    probRuleOut(i_group).cell_inds = cell_inds;
    probRuleOut(i_group).real_cellIDs = real_cellIDs;
    probRuleOut(i_group).setID = i_set;
    
    group_word_record = full_word_record(cell_groups(i_group, :), :);
    %% put inputs in number format
    conv_factor = 2.^(0:(group_size-1));
    word_nums = conv_factor*group_word_record;
    %% compute input efficiency
    [~, inputIpred, inputEntropy] = computePredictiveInfo(group_word_record', 0*group_word_record(1, :)', tau_future);
    input_efficiency = inputIpred(1)/inputEntropy(1);
    
    
    %% find which words are actually observed
    % word number "0" is 0000000...0 and has 'nz_word_inds' of 1
    obs_word_inds = unique(word_nums) + 1;
    obs_word_nums = obs_word_inds - 1;
    % start near a random perceptron
    wb_i = round(2*rand(1, group_size));
    [~, part_i] = perceptronToProjectionRule(wb_i, b);
    obs_proj_probs = part_i(obs_word_inds);
    %% compute input distributions 
    [xbins, ~, data_x] = unique(word_nums);
    data_x1 = data_x(1:end-1);
    data_x2 = data_x(2:end);
    n_x1 = length(xbins);
    n_x2 = length(xbins);
    fracs =  [1.0000    0.950   0.8500    0.7500];
    nreps = 50; 
    [inputIpred2, ~, ~, ~, ~, ~, ~, ~, Pxx1] = calc_info_P_joint(data_x1,data_x2,n_x1,n_x2,fracs,nreps);
    %%
    max_fr = mean(sum(group_word_record) ~= 0);
    fr_bins = linspace(0, max_fr, 11);
    fr_bins = fr_bins(2:end);
    
    clear frSimsOut
    frSimsOut(length(fr_bins)) = struct();
    %%
    tic
    for i_fr = 1:length(fr_bins)
%%
        fr_target = fr_bins(i_fr);
        % define important function for minimization / evaluation 
        prPI_frcon = @(x) probrulePredI_fr_objfun(x, obs_word_nums, word_nums, fr_lambda, fr_target, Pxx1) ;
        prPI_nocon =  @(x) probrulePredI_fr_objfun(x, obs_word_nums, word_nums, 0, fr_target, Pxx1)*-1 ;

        %% run fmincon from an array of initial positions satisfying firing-rate constraints. 
        P_xprior = sum(Pxx1);
        numRules = 50;
        
        per_wsample = X_per.effOut(i_group).init_w;
        % find the indices (in uniq-lists of items) of points close tot he
        % firing rate constraint
        good_startpoints = find(abs(per_wsample.uniwFR - fr_target) < 0.02);
        % want to sort the good_startpoints by their predictive info. to 
        % do this, first get the uniq-PI
        uniq_pI = per_wsample.all_info((per_wsample.wisUniq == 1));
        % now sort the uniq_pI at the good startpoints
        [~, ord] = sort(uniq_pI(good_startpoints), 'descend');
        good_startpoints = good_startpoints(ord);
        % now take up to four points to start the optimization
        partitionsStartPoints = per_wsample.uniqParts(good_startpoints(1:min(length(good_startpoints), 4)), :);
        %%
        
        x0s = generateFRConstrainedProbabilisticRules_fromInitialPoints(P_xprior, fr_target, partitionsStartPoints, numRules) ;
        
        xOut = zeros(size(x0s));
        efOut = zeros(numRules, 1);
        pIOut = zeros(numRules, 1);
        pIconOut = zeros(numRules, 1);
        
        
        parfor i_x0 = 1:numRules
            x0 = x0s(i_x0, :);
            lb = 0*x0;
            ub = 1 + 0*x0;
            ub(1) = 0.0001;

            options = optimoptions('fmincon', 'gradobj', 'on', 'display', 'off');
            [x, fv, exitflag] = fmincon(prPI_frcon,x0,[],[],[],[],lb,ub, [], options);
            
            xOut(i_x0, :) = x;
            efOut(i_x0) = exitflag;
            pIconOut(i_x0) = fv;
            pIOut(i_x0) = prPI_nocon(x);
        end
        % objective function is is: - ( I_wbi - fr_lambda*(w_fr - fr_target)^2 )

        frOut = sqrt((pIconOut + pIOut)/fr_lambda) + fr_target;
        
        
        frSimsOut(i_fr).xOut = xOut;
        frSimsOut(i_fr).efOut = efOut;
        frSimsOut(i_fr).pIconOut = pIconOut;
        frSimsOut(i_fr).pIOut = pIOut;
        frSimsOut(i_fr).frOut = frOut;
        frSimsOut(i_fr).fr_target = fr_target;
        frSimsOut(i_fr).fr_lambda = fr_lambda;
    end
    
    toc
    
    probRuleOut(i_group).frSimsOut = frSimsOut;
end
save(['probrule_grouprun_maxg' num2str(num_cells) '_set' num2str(i_set)], 'probRuleOut');
end

