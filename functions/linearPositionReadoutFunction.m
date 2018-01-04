function p_t = linearPositionReadoutFunction(s_t, k_x, C)


p_it = 0*s_t;
for i_input = 1:size(s_t, 2)
    p_it(:, i_input) = conv(s_t(:, i_input), k_x(i_input, :)', 'same');
end 

p_t = sum(p_it, 2) + C;