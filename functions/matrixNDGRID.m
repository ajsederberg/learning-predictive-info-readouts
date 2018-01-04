function xvn = matrixNDGRID(vec, is_fixed, fixed_vals)
% vec is the vector for creating the grid. is_fixed is a Nx1 logical, with
% 1 for all values of w that are not varied in the grid. fixed_vals are the
% fixed values that should replace 'vec' for those entries of w. 


num_dims = length(is_fixed);
is_fixed = reshape(is_fixed, [1 num_dims]);


x(num_dims) = nan;
if length(fixed_vals) == sum(is_fixed)
    x(is_fixed) = fixed_vals;
else
    x = fixed_vals;
end

% total_num_els = length(vec)^sum(~is_fixed);
repmat_vec = double(repmat(is_fixed, [num_dims 1]));
repmat_vec(repmat_vec == 0) = length(vec);
repmat_vec(repmat_vec == 1) = 1;

in_v = cell(num_dims, 1);
for ii= 1:num_dims
    if is_fixed(ii)
        rep_vec = reshape(x(ii), ones(1, num_dims));
    else
        dim_vec = ones(1, num_dims);
        dim_vec(ii) = length(vec);
        rep_vec = reshape(vec, dim_vec);
    end
    repmat_vec(ii,ii) = 1;
    in_v{ii} = repmat(rep_vec, repmat_vec(ii, :));
end


xvn = zeros(numel(in_v{1}), num_dims);
for ii = 1:num_dims
    xvn(:, ii) = in_v{ii}(:);
end