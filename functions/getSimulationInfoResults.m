function [all_lt_res, all_lt_res_info, processed_results] = getSimulationInfoResults(set_size, group_tag, learning_types, learn_param)
% getSimulationInfoResults looks up the simulation results and info
% calculations for each set with a simulation result and returns a cell
% array with collated results, as well as basic processed results in the
% struct processed_results



%%

% group_tag = ['fA' num2str(set_size)];

X_learn = load(['aggregated results/' group_tag '_learning_results.mat']);
X_learn = X_learn.res_ds_arr;

if strcmp(group_tag(1), 'c')
    i_group_tag = ['f' group_tag(2:end)];
else
    i_group_tag = group_tag;
end
X_info = load(['aggregated results/' i_group_tag '_info_calcs.mat']);
X_dict = X_info.perc_lookup_struct;
X_info = X_info.info_ds_arr;

all_sets_wordprobs = load(['aggregated results/word_probabilites_' i_group_tag '.mat']);
all_set_ids = cellfun(@(x) x.set_ID, X_info);
X_wordprobs = all_sets_wordprobs.word_probs(all_set_ids);

%% For each learning type, compute the fraction of total pred-I learned and
% the similarity to the optimal rule. 
if nargin == 2
    learning_types = {'vanilla', 'triplet', 'vanilloh', 'triploh'};
end
prct_learned = cell(length(learning_types), 1); %  fraction of total pred-I learned
prct_maxFR = cell(length(learning_types), 1);   %  fraction of max FR of learned perc
learned_pI = cell(length(learning_types), 1); 
learned_perc_ID = cell(length(learning_types), 1);
rule_sims = cell(length(learning_types), 1);    %  similarity to optimal rule
all_perc_sims = cell(length(learning_types), 1);    % similarity of all perceptrons to their optimal rule
all_lt_res = cell(length(learning_types), 1);
all_lt_res_info = cell(length(learning_types), 1);
all_lt_perc_fr = cell(length(learning_types), 1);
all_lt_perc_pI = cell(length(learning_types), 1);
all_lt_last10dW = cell(length(learning_types), 1);   % SE of weights over last five timesteps: check. 

% weights of the hull for each set
all_lt_wb_hull = cell(length(learning_types), 1);
all_lt_wb_hull_perc_ID = cell(length(learning_types), 1);

all_lt_res_hull_bin = cell(length(learning_types), 1);
lt_fr_bins = cell(length(learning_types), 1);
 
for i_lt = 1:length(learning_types)
    %%
    % extract all results for specified learning type. 
    if nargin > 3 && ~isempty(learn_param{i_lt})
        %%
        lt_res_sel = cellfun(@(x) strcmp(x.param.learning_type, learning_types{i_lt}), X_learn);

        % extract the field names of the parameter that is specified by
        % learn_param. This is how you select for a particular value of
        % alphaLTD, for example. Other options: 'noise_sig_val', 'tau_y',
        % 'wmax'
        learn_param_field = fieldnames(learn_param{i_lt});
        for i_lpf = 1:length(learn_param_field)
            new_condition = cellfun(@(x) x.(learn_param_field{i_lpf}) == ...
                learn_param{i_lt}.(learn_param_field{i_lpf}), X_learn);
            lt_res_sel = lt_res_sel & new_condition;
        end


    else
        % if no other parameters are specified, then select for
        % learning_type only
        lt_res_sel = cellfun(@(x) strcmp(x.param.learning_type, learning_types{i_lt}), X_learn);
    end

    lt_res = X_learn(lt_res_sel);
    
    % get the set numbers, and re-order results by set number
    lt_setnums = cellfun(@(x) x.set_ID, lt_res);
    [lt_setnums, ord] = sort(lt_setnums);
    lt_res = lt_res(ord);
    
    % get the info results for these sets. Treat X_wordprobs just like
    % X_info
    keep_info_ent = cellfun(@(x) ismember(x.set_ID, lt_setnums), X_info);
    lt_res_info = X_info(keep_info_ent);
    lt_res_wordprobs = X_wordprobs(keep_info_ent);
    
    lt_info_setnums = cellfun(@(x) x.set_ID, lt_res_info);
    [lt_info_setnums, ord] = sort(lt_info_setnums);
    lt_res_info = lt_res_info(ord);
    lt_res_wordprobs = lt_res_wordprobs(ord);
    
    % if any info calculations are missing, drop those results
    keep_res_ent = ismember(lt_setnums, lt_info_setnums);
    lt_res = lt_res(keep_res_ent);
    lt_setnums = lt_setnums(keep_res_ent);
    
    display(['info calcs and learning results are ' ...
        num2str(100*mean(lt_setnums == lt_info_setnums), '%1.0f') '% matched'])
    all_lt_res{i_lt} = lt_res;
    all_lt_res_info{i_lt} = lt_res_info;
    
    % now that we have the sets matched up:
    if set_size < 6
        wb_perc_partitions = X_dict.percX.unique_partitions;
        wb_perc = X_dict.percX.unique_perceptrons;
    else
        wb_perc_partitions = X_dict.percX.final_partitions;
        wb_perc = X_dict.percX.final_perceptrons;
    end       

    
    % get the final weights
    lt_res_finalWB = cellfun(@(x) squeeze(x.trialW(end, :, :)), lt_res, 'uniformoutput', false);
    % % look up the pred-I of the final state
    lt_res_finalPartition = cellfun(@(x) perceptronToProjectionRule(x, 1), ...
        lt_res_finalWB, 'uniformoutput', false);
    [~, lt_final_perc_ind] = cellfun(@(x) checkPartitionUniqueness_RPMethod(wb_perc_partitions, x), ...
        lt_res_finalPartition, 'uniformoutput', false);



    % compute the perceptron hull for each group
    % [pI_low, pI_high, fr_bins, wb_hull] = computePerceptronHull_withWeights(fr_edges, fr_perc, pI_perc, wb_perc)
    fr_perc_all = cellfun(@(x) x.X_pw.fr_perc, lt_res_info, 'uniformoutput', false);
    pI_perc_all = cellfun(@(x) x.X_pw.pI_perc(:, 2), lt_res_info, 'uniformoutput', false);
    
    all_lt_perc_fr{i_lt} = fr_perc_all;
    all_lt_perc_pI{i_lt} = pI_perc_all;
    

       
    
    all_fr_max = max(cellfun(@max, fr_perc_all));
    fr_edges = linspace(0, all_fr_max, 201);
    [~, pI_high, fr_bins_all, wb_hull_all] = cellfun(@(x, y) computePerceptronHull_withWeights(fr_edges, x, y, wb_perc), ...
        fr_perc_all, pI_perc_all, 'uniformoutput', false);
    
    all_lt_wb_hull{i_lt} = wb_hull_all;
    % find the indices of the hull in the full set of sampled perceptrons
    
    wb_hull_Partition = cellfun(@(x) perceptronToProjectionRule(x, 1), ...
        wb_hull_all, 'uniformoutput', false);
    [~, wb_hull_perc_ind] = cellfun(@(x) checkPartitionUniqueness_RPMethod(wb_perc_partitions, x), ...
        wb_hull_Partition, 'uniformoutput', false);
    all_lt_wb_hull_perc_ID{i_lt} = wb_hull_perc_ind;
    
    num_fr_bins = length(fr_bins_all{1});
    lt_fr_bins{i_lt} = fr_edges;
    % compute learned perceptron similarity to nearest point on hull
    % % first find nearest hull point to the final state
    % find the firing rate on the fishmovie (perc) of the learned
    % perceptrons
    try lt_res_finalFR = cellfun(@(x, ind_x) x(ind_x), fr_perc_all, lt_final_perc_ind, ...
        'uniformoutput', false);
    catch
            % if that fails, take the empirical end firing rate
        lt_res_finalFR = cellfun(@(x) x.trialFR(end, :), lt_res, 'uniformoutput', false);
    end
    [lt_res_cts, ~, lt_res_hull_bin] = cellfun(@(x) histcounts(x, fr_edges), ...
        lt_res_finalFR, 'uniformoutput', false);
    %% correct bins by 1
    lt_res_hull_bin = cellfun(@(x) min(x+1, num_fr_bins), lt_res_hull_bin, 'uniformoutput', false);    
    all_lt_res_hull_bin{i_lt} = lt_res_hull_bin;
    
    % compute *all* perceptrons similarity to nearest point on hull 
    % 
    [~, ~, perc_wb_bin] = cellfun(@(x) histcounts(x, fr_edges), fr_perc_all, ...
        'uniformoutput', false);
    perc_wb_bin = cellfun(@(x) min(x+1, num_fr_bins), perc_wb_bin, 'uniformoutput', false);
    
    % calculate how much final weights changed in the last NN steps : in
    % steady state? 
    lt_res_last10Wstd = cellfun(@(x) mean(var(x.trialW(end-9:end, :, :)), 3), ...
        lt_res, 'uniformoutput', false);
    all_lt_last10dW{i_lt} = lt_res_last10Wstd;
    
    % compare the final weights to the appropriate hull weights
    rule_sims{i_lt} = cellfun(@(x,y,z, w) computePerceptronStructuralSimilarity(...
        x(y, :), z, w), wb_hull_all, lt_res_hull_bin, lt_res_finalWB, lt_res_wordprobs,  ... 
        'uniformoutput', false);
    all_perc_sims{i_lt} = cellfun(@(x,y, w) computePerceptronStructuralSimilarity(...
        x(y, :), wb_perc, w), wb_hull_all, perc_wb_bin, lt_res_wordprobs,  ... 
        'uniformoutput', false);
    %%
    % compute percent of hull pred-I learned 

    
    % report total number of not-found partitions
    missing_partitions = sum(cellfun(@(x) sum( x== 0), lt_final_perc_ind));
    display([num2str(missing_partitions) ' were missing from the sampled perceptron set'])
    % if not found, assign to first bin
    lt_final_perc_ind = cellfun(@(x) max(x, 1), lt_final_perc_ind, 'uniformoutput', false);
    learned_perc_ID{i_lt} = lt_final_perc_ind;
    
    
    lt_res_finalPI = cellfun(@(x, y) x(y), pI_perc_all, lt_final_perc_ind, 'uniformoutput', false);
    % % look up the pred-I of the nearest hull point
    lt_res_hull_pI_all = cellfun(@(x, y) x(y)', pI_high, lt_res_hull_bin, ...
        'uniformoutput', false);
    learned_pI{i_lt} = lt_res_finalPI;
    prct_learned{i_lt} = cellfun(@(x, y) 100*x./y, lt_res_finalPI, lt_res_hull_pI_all, ...
        'uniformoutput', false);
    prct_maxFR{i_lt} = cellfun(@(x) 100*x'/max(x), lt_res_finalFR, 'uniformoutput', false);
    
    
end


processed_results.group_tag = group_tag;
processed_results.learning_types = learning_types;
processed_results.learned_pI = learned_pI;
processed_results.learned_perc_ID = learned_perc_ID;
processed_results.all_lt_last10dW = all_lt_last10dW;
processed_results.all_lt_wb_hull = all_lt_wb_hull;
processed_results.all_lt_wb_hull_perc_ID = all_lt_wb_hull_perc_ID;
processed_results.prct_learned = prct_learned;
processed_results.prct_maxFR = prct_maxFR;
processed_results.rule_sims = rule_sims;
processed_results.all_perc_sims = all_perc_sims;
processed_results.pI_high = pI_high;
processed_results.lt_fr_bins = lt_fr_bins;
processed_results.prct_learned = prct_learned;
processed_results.all_lt_res_hull_bin = all_lt_res_hull_bin;
processed_results.wb_hull_all = wb_hull_all;
processed_results.wb_perc_partitions = wb_perc_partitions;
processed_results.all_lt_perc_fr = all_lt_perc_fr;
processed_results.all_lt_perc_pI = all_lt_perc_pI;
