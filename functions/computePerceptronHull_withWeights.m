function [pI_low, pI_high, fr_bins, wb_hull] = computePerceptronHull_withWeights(fr_edges, varargin)
% [pI_low, pI_high, fr_bins, wb_hull] = computePerceptronHull_withWeights(fr_edges, varargin)
%  options: pass the firing rate bin edges in fr-edges, and a struct with
%  fields 'fr_perc' (firing raes of percetprons) and 'pI_perc' (predictive
%  information of perceptrons) and 'allowed_wb' (perceptron weights) 
%           OR
%  pass fr_perc, pI_perc, wb_perc as inputs 2, 3, and 4. 


if nargin == 2
    fr_perc = varargin{1}.fr_perc;
    pI_perc = varargin{1}.pI_perc;
    wb_perc = varargin{1}.allowed_wb;
    
    % round fr down to nearest 1000th
    fr_perc = floor(fr_perc*1000)/1000;
    % round pI up to nearest 1000th
    pI_perc = ceil(pI_perc*1000)/1000;

    
else
    fr_perc = varargin{1};
    pI_perc = varargin{2}; 
    wb_perc = varargin{3};

end

% fr_edges defines the bin divisions of the firing rate axis
% fr_bins will have the firing rate of the optimal perceptron in each bin
fr_bins = edges2bins(fr_edges);
% pI_low is the minimum pred-I within the bin
pI_low = 0*fr_bins;
% pI_high is the maximum pred-I within the bin
pI_high = 0*fr_bins;
% wb_hull contains the weights of the perceptron with highest pred-I
wb_hull = zeros(length(fr_bins), size(wb_perc, 2));

start_bin = sum(fr_edges <= min(fr_perc)) + 1;
for i_fr = start_bin:length(fr_bins)
    fr_sel_perc = fr_perc < fr_bins(i_fr);
    
    bin_pI = pI_perc(fr_sel_perc);
    bin_fr = fr_perc(fr_sel_perc);
    bin_wb = wb_perc(fr_sel_perc, :);
    
    % maximum at that firing rate or lower
    if  sum(fr_sel_perc)>0
        [pI_high(i_fr), ind] = max(bin_pI);
%         fr_bins(i_fr) = bin_fr(ind);
        wb_hull(i_fr, :) = bin_wb(ind, :);
    elseif i_fr > 1
        pI_high(i_fr) = pI_high(i_fr-1);
        wb_hull(i_fr, :) = wb_hull(i_fr-1, :);
    end
    
    
    % minimum at that firing rate or higher
    fr_sel_perc = fr_perc >= fr_bins(i_fr);
    
    bin_pI = pI_perc(fr_sel_perc);
    if sum(fr_sel_perc)>0
        pI_low(i_fr) = min(bin_pI);
    end
    
end