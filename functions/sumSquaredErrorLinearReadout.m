function ss_err = sumSquaredErrorLinearReadout(k_x_vec, s_t, x_t)

C = k_x_vec(1);

k_x = reshape(k_x_vec(2:end), size(s_t, 2), numel(k_x_vec(2:end))/size(s_t, 2));

p_t = linearPositionReadoutFunction(s_t, k_x, C);

ss_err = mean((p_t - x_t).^2) + 5*sum(sum(abs(diff(k_x, 1, 2)).^2));

% grad_ss_err = sum(2*(p_t - x_t)*dpdk);