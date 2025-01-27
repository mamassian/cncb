% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.2
%
% cncb_get_bnds
%   This function is useful to find confidence boundaries so that the
%   number of confidence judgments has a particular distribution across
%   all levels. For now, this distribution can be uniform (default), or
%   triangular with a linear increase (or decrease) for consective
%   confidence levels. The parameter 'conf_bias' describes the ratio of 
%   confidence judgments between the top and bottom levels. If it is larger
%   than 1, there will be a dominance of high confidence values.
%
%
% INPUT:
%   'simul_params': structure of simulation parameters:
%       'sens_intens': difficulty levels
%       'nb_trials'  : number of confidence judgments to simulate
%
%   'model_params': model values parameters, as a structure
%       'sens_noise' : sensory (Type 1) sdtev of noise (0 = perfectly sensitive)
%       'sens_crit'  : sensory (Type 1) criterion
%       'conf_noise' : confidence (Type 2) sdtev of noise (0 = ideal)
%       'conf_boost' : fraction super-ideal (0 = ideal, 1 = super-ideal)
%       'conf_crit'  : confidence (Type 2) criterion
%       'conf_bnds'  : confidence boundaries
%
%   'conf_nbbins': confidence nb of levels
%
%   OPTIONAL:
%   'conf_bias'  : confidence bias (ratio of conf. choices b/w top & bottom levels)
%
% OUTPUT: 
%   'conf_bnds'      : confidence boundaries
%
%
% 02-JAN-2025 - pascal mamassian


function conf_bnds = cncb_get_bnds(simul_params, model_params, ...
    conf_nbbins, conf_bias)

    
    if ~exist('conf_bias', 'var')
        % -> by default, uniform distribution across all confidence levels
        conf_bias = 1.0;
    end


    model_orig_params2 = model_params;
    model_orig_params2.conf_continuous = true;
    
    % -> do 1 run just to get the distribution of confidence values
    crs_orig_data = cncb_simul(simul_params, model_orig_params2);
    
    % -> ratio of confidence choices between top and bottom confidence levels
    %    frac_bin_top / frac_bin_bottom = conf_bias

    quant_a = 2 * (conf_bias - 1) / ...
        ((1 + conf_bias) * conf_nbbins * (conf_nbbins + 1));
    quant_b = 2 * (conf_nbbins + 1 - conf_bias) / ...
        ((1 + conf_bias) * conf_nbbins * (conf_nbbins + 1));
    quant_vals = NaN(1, conf_nbbins);
    for qq = 1:conf_nbbins
        quant_vals(qq) = quant_a * qq*(qq + 1)/2 + quant_b * qq;
    end

    % quant_vals = linspace(0, 1, rating_bnds_nb+2);
    % quant_vals = quant_vals(2:(end-1));
    conf_bnds = NaN(1, conf_nbbins);    % redefine bounds
    

    conf_vals = crs_orig_data.raw_data(:,3);

    for kk = 1:conf_nbbins
        q_val = quant_vals(kk);
        c_kk = quantile(conf_vals, q_val);
        conf_bnds(kk) = c_kk;
    end


end


% -> THE END <-
