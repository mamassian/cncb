% log-odds transform
%
% based on:
% Zhang, H., & Maloney, L. T. (2012). Ubiquitous log odds: a common
%   representation of probability and frequency distortion in perception,
%   action, and cognition. Frontiers in Neuroscience, 6, 1. 
%   https://doi.org/10.3389/fnins.2012.00001 
%
%   The parameter 'param_p0' in the log-odds function is the "fixed point" 
% of the transformation, the value of p which is mapped to itself. 
%   The parameter 'param_gamma', is the slope of the linear transformation
% on log odds scales, and on linear scales, is the slope of the curve at
% the crossover point 'param_p0'.
%   Identity mapping is obtained for 'param_gamma = 1'.
%
% USAGE:
%   new_prob = cncb_log_odds(0.4, 0.75, 0.4);
%   prob_lst = cncb_log_odds(0:0.01:1, 0.75, 0.4);
%
% NOTE:
%   The inverse function is obtained from:
%   old_prob = cncb_log_odds(new_prob, 1/param_gamma, param_p0)
%
% 19-SEP-2021 - pascal mamassian

function new_prob = cncb_log_odds(old_prob, param_gamma, param_p0)

    % -> log-odds function of a probability
    fun_lo = @(pp) log(pp ./ (1 - pp));

    % -> inverse log-odds function
    fun_lo_inv = @(lo) 1 ./ (1 + exp(-lo));

    yy =  param_gamma * fun_lo(old_prob) + (1 - param_gamma) * fun_lo(param_p0);
    
    new_prob = fun_lo_inv(yy);
    
end
