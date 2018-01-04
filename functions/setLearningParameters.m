function param = setLearningParameters(learning_type, set_size, data_source, varargin)

% learning_type: vanilla, vanilloh, triplet, triploh
param.learning_type = learning_type;

% data_source: this can be 'fishmovieS', 'fishmovieA', 'multimovieS',
% 'multimovieA', 'checkerboardA', or 'checkerboardS'
param.data_source = data_source;

% choose learning type: vanilla, vanilloh, triplet, triploh
% learning_type = 'vanilla';

% this can be set to 3, 4, 5... 10
param.numL0 = set_size;


% restrict to positive weights
param.wmin = 0;

% speed of learning, length of interval at which to record progress (every
% NT trials, weights are recorded), and number of times to repeat the
% training stimulus 
param.alphaLTP = 0.01;
param.training_reps = 200;



% set varargin params -- this overwrites anything already set. 
for i_p = 1:length(varargin)/2
    param.(varargin{2*i_p - 1}) = varargin{2*i_p};
end

% default is 'fish', 'nrb' is also an option
if ~isfield(param, 'movie_type') 
    if strcmp(data_source(1:end-1), 'checkerboard')
        param.movie_type = 'checkerboard';
    else
        param.movie_type = 'fish';
    end
end


% set defaults, if not already set in varargin
if ~isfield(param, 'numL1')
    if ~isfield(param, 'w0_values')
    % this is the number of initial conditions to test
        param.numL1 = 200;
    else
        param.numL1 = size(param.w0_values, 1);
    end
end


if ~isfield(param, 'w0_uniform') 
    % w0_uniform: fix the random number generator and set w0 uniform randomly between w0min and w0max
    % w0_uniform = false: fix the random number generator and generate w0
    % according to b/2 + 0.3*randn, with negative values zeroed. 
    % If w0_uniform is not set and if w0max/w0min are not provided, then set
    % w0_uniform to false
    if ~isfield(param, 'w0max') || ~isfield(param, 'w0min')
        param.w0_uniform = false; 
    else
        param.w0_uniform = true;
    end
elseif param.w0_uniform && ~isfield(param, 'w0max')
    % if w0_uniform is true, but w0max is not provided, set defaults. 
    % this is the initial condition for weights
    param.w0max = 0.9;
    % setting w0min to 0.5 means any simultaneous firing in the inputs will
    % kick off learning 
    param.w0min = 0.5;
end

if strcmp(learning_type, 'vanilloh') || strcmp(learning_type, 'triploh')
    % for vanilloh training, want more variance across initial conditions.
    % i.e. permit the entire range, 0 to 1.1, if uniform distribution used.
    % 
    if param.w0_uniform
       param.w0max = 1.1;
       param.w0min = 0;
    end
end
% this is the maximum weight: either sub or super threshold (<1 or >1)
if strcmp(learning_type, 'vanilloh') || strcmp(learning_type, 'triploh')
    param.wmaxes = round(0.75*set_size*10)/10;   % for 'vanilloh' this is the sum of w, not maximum value of w
else
    param.wmaxes = 1.1; %
                        % edit 4/7/16: sub-threshold is not interesting with vanilla
                        % dropped it.   
end




% set learning type-dependent parameters
if strcmp(learning_type, 'vanilla') || strcmp(learning_type, 'vanilloh')
    if ~isfield(param, 'alphas_ltd')
        param.alphas_ltd = 1.2; 
    end
    
    if ~isfield(param, 'noise_sigmas')
        param.noise_sigmas = 0;
    end
    
elseif strcmp(learning_type, 'triplet') || strcmp(learning_type, 'triploh')
    if ~isfield(param, 'alphas_ltd')
        if set_size == 4
            param.alphas_ltd = 0.7;
        else
            param.alphas_ltd = 0.9;
        end
    end
    if ~isfield(param, 'tau_y')
        param.tau_y = 7;
    end

end

