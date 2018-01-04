% generate sets of perceptrons and save
clear;
num_tries = 1000;   % if the guess is correct, this can be very high w/o affecting
                    % runtime much, since all guessed partitions will be
                    % found. If guess is too low, however, process aborts. 
num_guesses = [0 0 19 149 650 3000 10000 15000 15000 15000];

for cells_in_set = 3:10
    
        %% find a selection of allowable and distinct perceptrons
    num_wbs_guess = num_guesses(cells_in_set);    % ~15 for 3 cells; 150 for 4 cells; many for 5 cells

    allowed_wb = zeros(num_wbs_guess, cells_in_set);
    allowed_partitions = zeros(num_wbs_guess, 2^cells_in_set);
    ap_uniq_num = nan(num_wbs_guess, 1);


    % instead of doing something systematic, let's just try random w's
    % and then vary b from just larger than min(w) up to sum(w)
    % sample w's with 0.1-precision
    p_ctr = 0;


    for i_try = 1:num_tries
        if p_ctr < num_wbs_guess
            w = ceil(100*rand(cells_in_set, 1))/100;

            % gradually increase b for this direction of w
            b_vals = (min(w)-.05):0.1:sum(w);
            for i_b = 1:length(b_vals)
                b = b_vals(i_b);

                norm_w = w/b;
                % choose a value just larger than 1 (b = 1 for normalized w).
                norm_w(norm_w > 1.1) = 1.1;
                [~, new_apic] = perceptronToProjectionRule(norm_w, 1);

                % this function checks new_apic against the row
                % vectors of allowed_partitions
                is_new_part = checkPartitionUniqueness(...
                    allowed_partitions(1:p_ctr, :), new_apic);

                if is_new_part && p_ctr < num_wbs_guess
                    p_ctr = p_ctr + 1;
                    allowed_partitions(p_ctr, :) = new_apic;
                    % normalize by b (ie b == 1)
                    allowed_wb(p_ctr, :) = norm_w;
                end
            end
        
            if mod(i_try, 10) == 0
                display([num2str(p_ctr) ' perceptrons found in ' num2str(i_try) ' tries'])
            end
        end
    end
    allowed_partitions = allowed_partitions(1:p_ctr, :);
    allowed_wb = allowed_wb(1:p_ctr, :);
    display([num2str(p_ctr) ' perceptrons found in ' num2str(i_try) ' tries'])
    
    
    save(['postive_perceptrons_dim' num2str(cells_in_set)], 'allowed_wb', 'allowed_partitions');
end