function sim_vals = computePerceptronStructuralSimilarity(perceptron_1, perceptron_2, sim_type)
% this computes similarity based on the partition rule. it counts the
% number of partitions that are identical for perceptron1(1, :) and
% perceptron2(1, :). The sim_type can be used to control which partitions
% count by setting a maximum spikecount for counted words. 

parts_1 = perceptronToProjectionRule(perceptron_1, 1);
[parts_2, word_mat] = perceptronToProjectionRule(perceptron_2, 1);

partition_spikecount = sum(word_mat);

% select which of the rules to compare. Generally, the zero-count word is
% ignored, and anything with 3+ spikes is ignored (since most rules agree
% on these)


if ~isnumeric(sim_type)
    obs_words = double(sim_type.words');
    [~, word_inds] = checkPartitionUniqueness_RPMethod(word_mat', obs_words');
    
    obs_parts_1 = parts_1(:, word_inds);
    obs_parts_2 = parts_2(:, word_inds);
    
    word_probs = sim_type.probs;
    
    word_probs_n0 = word_probs/sum(word_probs(2:end));
    word_probs_n0(1) = 0;
    
    sim_vals = (double(obs_parts_1 == obs_parts_2))*(word_probs_n0);

else
    if sim_type > 0

        compare_parts = partition_spikecount > 0 & partition_spikecount <= sim_type;

    else
        compare_parts = true(size(partition_spikecount));
    end
    sim_vals = mean(parts_1(:, compare_parts) == parts_2(:, compare_parts), 2);
end
