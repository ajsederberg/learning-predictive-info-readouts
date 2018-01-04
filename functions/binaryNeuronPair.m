function [x, w] = binaryNeuronPair(x,w,inp,syn_sigma, b, learning_type, varargin)
% x is Nx1. w is Nxlength(inp). inp is length(inp)x1. 
% if g_I = 0, no connections between neurons and simply a way to test
% multiple ICs. 
% variable arguments in: depend on learning type

% lateral inhibition? comment out for now
% g_I = 0;
% jij = -g_I*(ones(length(x)) - eye(length(x))); % strength of lateral inhibition
if length(x) ~= size(w, 1)
    keyboard
end
% save previous value of activity 
x0 = x;
% add some random synaptic noise 
syn_noise = randn(size(x))*syn_sigma;
% x = (w*inp + jij*x0 + syn_noise) > b;
% compute updated x
x = (w*inp + syn_noise) > b;
if length(x) ~= size(w, 1)
    keyboard
end
if ~isempty(learning_type)
    if strcmp(learning_type, 'bcmlike')
        alpha_ltp = varargin{1};
        eps_ltd = varargin{2}*alpha_ltp;
        target_fr = varargin{3};
        ave_fr = varargin{4};
        

        % activity-dependent change
        dw = diag(target_fr - ave_fr)*(x*inp'*alpha_ltp - eps_ltd*x0*inp');
        w = w + dw;

        w(w<0) = 0;
    elseif strcmp(learning_type, 'homeostatic')
        
        alpha_ltp = varargin{1};
        eps_ltd = varargin{2}*alpha_ltp;
        wmax = varargin{3};

        if any(w(:) < 0)
            keyboard
        end

        % vanilla 
        dw = (x*inp'*alpha_ltp) - eps_ltd*x0*inp';

        % homeostatically regulate total synaptic weight
        homeo_fac = sign(diag(wmax - sum(w, 2))).*(abs(diag(wmax - sum(w, 2))).^0.1);
        w = w + homeo_fac*(w.^0.1).*(dw);

        if any(imag(w(:)) ~= 0)
            keyboard 
        end
        
    elseif strcmp(learning_type, 'constantLTD')
        alpha_ltp = varargin{1};
        eps_ltd = varargin{2}*alpha_ltp;
        wmax = varargin{3};
        if nargin > 6 + 4   % 6 inputs before varargin, 4 elements of varargin
            wmin = varargin{4};
        else
            wmin = 0;
        end

        % vanilla STDP
        dw = (x*inp'*alpha_ltp) - eps_ltd; 

        % hard bounds
        w = w + dw;
        w(w < wmin) = wmin; % bound by wmin from below
        w(w > wmax) = wmax; % bound by wmax from above

    elseif strcmp(learning_type, 'vanilla') 
        
        alpha_ltp = varargin{1};
        eps_ltd = varargin{2}*alpha_ltp;
        wmax = varargin{3};
        if nargin > 6 + 4   % 6 inputs before varargin, 4 elements of varargin
            wmin = varargin{4};
        else
            wmin = 0;
        end

        % vanilla STDP
        dw = (x*inp'*alpha_ltp) - eps_ltd*x0*inp';
        

        % hard bounds
        w = w + dw;
        w(w < wmin) = wmin; % bound by wmin from below
        w(w > wmax) = wmax; % bound by wmax from above
        
        
    elseif strcmp(learning_type, 'vanilloh')
        % "vanilla with online homeostasis"
        % same input structure as vanilla, but the third input give the
        % weight-sum constraint instead of wmax. 
        
        alpha_ltp = varargin{1};
        eps_ltd = varargin{2}*alpha_ltp;
        fixed_wsum = varargin{3};


        

        % vanilla STDP
        dw = (x*inp'*alpha_ltp) - eps_ltd*x0*inp';
        

        w = w + dw;
      
        % weights must be positive
        wmin = 0;
        w(w < wmin) = wmin; % bound by wmin from below
%         w(w > 1.1*b) = 1.1*b; % bound by 110% of bias (b) from above

        % constrain sum of weights 
        sum_w = sum(w, 2);
        w = diag(fixed_wsum./sum_w)*w;

    elseif strcmp(learning_type, 'antiSTDP') 
        
        alpha_ltp = varargin{1};
        eps_ltd = varargin{2}*alpha_ltp;
        wmax = varargin{3};
        if nargin > 6 + 4   % 6 inputs before varargin, 4 elements of varargin
            wmin = varargin{4};
        else
            wmin = 0;
        end

        % only difference between this and vanilla STDP is the minus sign
        % below. 
        dw = - ( (x*inp'*alpha_ltp) - eps_ltd*x0*inp' );
        

        % hard bounds
        w = w + dw;
        w(w < wmin) = wmin; % bound by wmin from below
        w(w > wmax) = wmax; % bound by wmax from above
    elseif strcmp(learning_type, 'vanilla_homeostasisstep')
        alpha_ltp = varargin{1};
        target_fr = varargin{2};
        ave_fr = varargin{3};

        % activity-dependent change that doesn't depend on presynaptic
        % input
        dw = alpha_ltp*diag(target_fr - ave_fr)*x*ones(size(inp))';
        w = w + dw;

        w(w<0) = 0;
        
    elseif strcmp(learning_type, 'triplet')
        
        alpha_ltp = varargin{1};
        eps_ltd = varargin{2}*alpha_ltp;
        wmax = varargin{3};
        wmin = varargin{4};
        timeSinceLastSpike = varargin{5};
        tau_y = varargin{6};  % this is in time bins
        
        % turn off influence of exp(tSLS) until first spike
        timeSinceLastSpike(isnan(timeSinceLastSpike)) = 0;  
        spk2exp = exp(-timeSinceLastSpike/tau_y);
        
        % vanilla STDP
        dw = ((x.*spk2exp)*inp'*alpha_ltp) - eps_ltd*x0*inp';

        

        % hard bounds
        w = w + dw;
        w(w < wmin) = wmin; % bound by wmin from below
        w(w > wmax) = wmax; % bound by wmax from above
    elseif strcmp(learning_type, 'triploh')
        
        alpha_ltp = varargin{1};
        eps_ltd = varargin{2}*alpha_ltp;
        fixed_wsum = varargin{3};
        wmin = 0; %varargin{4};
        timeSinceLastSpike = varargin{5};
        tau_y = varargin{6};  % this is in time bins
        
        % turn off influence of exp(tSLS) until first spike
        timeSinceLastSpike(isnan(timeSinceLastSpike)) = 0;  
        spk2exp = exp(-timeSinceLastSpike/tau_y);
        
        % vanilla STDP
        dw = ((x.*spk2exp)*inp'*alpha_ltp) - eps_ltd*x0*inp';

        

        w = w + dw;
      
        % weights must be positive
        w(w < wmin) = wmin; % bound by wmin from below
        w(w > 1.1*b) = 1.1*b; % bound by 110% of bias (b) from above

        % constrain sum of weights 
        sum_w = sum(w, 2);
        w = diag(fixed_wsum./sum_w)*w;
    
    end
end

