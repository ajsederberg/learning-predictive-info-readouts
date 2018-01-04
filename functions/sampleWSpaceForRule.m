function [rule_ws, otherrule_ws] = sampleWSpaceForRule(perc_start, b, num_ws, max_dist)


if nargin == 2
    max_dist = 0.1;
    num_ws = 100;
end


%%

perc_rule = perceptronToProjectionRule(perc_start, b);
sat_wvals = false(size(perc_start));  

vec = linspace(-max_dist, max_dist, ceil(num_ws^(1/sum(~sat_wvals))));
dw = matrixNDGRID(vec, sat_wvals, 0);

num_ws = size(dw, 1);

set_of_ws = repmat(perc_start, num_ws, 1) + dw; 
set_of_ws(set_of_ws < 0) = 0;
%%
set_of_rules = perceptronToProjectionRule(set_of_ws, b);

same_rule_ws = all(set_of_rules == repmat(perc_rule, num_ws, 1), 2);

rule_ws = set_of_ws(same_rule_ws, :);
otherrule_ws = set_of_ws(~same_rule_ws, :);