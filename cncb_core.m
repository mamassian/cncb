% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.2
%
% cncb_core
%   This function generates predictions of the model for confidence ratings
%   for a Type 1 discrimination task. 
%
%
% INPUT:
%   'stim_lst'    : vector of stimuli (ordered)
%   'model_params': model values parameters, as a structure
%       'sens_noise' : sensory (Type 1) sdtev of noise (0 = perfectly sensitive)
%       'sens_crit'  : sensory (Type 1) criterion
%       'conf_noise' : confidence (Type 2) sdtev of noise (0 = ideal)
%       'conf_boost' : fraction super-ideal (0 = ideal, 1 = super-ideal)
%       'conf_crit'  : confidence (Type 2) criterion
%       'conf_continuous'  : are confidence ratings on continuous scale?
%       'conf_half_scale'  : if continuous, is it on a half-scale?
%       'conf_bnds'  : if continuous, confidence boundaries
%       'conf_cont_range'  : if continuous, range
%       'conf_cont_nb_levels'  : if continuous, resolution to fit
%
% OUTPUT: cncb_rating_mat 
%    matrix where each line is [stim, resp, conf, prob(resp, conf | stim)]
%
%
% 19-AUG-2021 - pascal mamassian
% 16-FEB-2022 - pm: added continuous ratings
% 21-OCT-2023 - pm: fixed continuous ratings
% 08-JAN-2024 - pm: change 1st input argument to 'stim_lst'
% 06-MAR-2024 - pm: break continuous ratings by quantiles


function cncb_rating_mat = cncb_core(stim_lst, model_params)

stim_nb = length(stim_lst);


% -> PARAMETERS <-

% -> Type 1 parameters
sens_noise = model_params.sens_noise;  % sensory noise
sens_crit  = model_params.sens_crit;   % sensory criterion

% -> Type 2 parameters (compulsory)
conf_noise = model_params.conf_noise;  % confidence noise
conf_boost = model_params.conf_boost;  % confidence boost

% -> default values
min_conf_noise = 1e-6;
default_conf_crit = 0.0;
default_continuous = 1;
default_gamma = 1.0;
default_p0 = 0.5;
default_conf_range = [0.0, 1.0];
default_conf_half_scale = false;
default_conf_levels_nb = 10;


% -> avoid singular covariance matrix 'lcl_distA_type1n2_var'
conf_noise = max(min_conf_noise, conf_noise);

% -> Type 2 parameters (optional)
% -> confidence criterion
if (any(strcmp(fieldnames(model_params), 'conf_crit')))
    conf_crit  = model_params.conf_crit;   
else
    conf_crit = default_conf_crit;
end


% -> use discrete or continuous confidence judgments
if (any(strcmp(fieldnames(model_params), 'conf_continuous')))
    conf_continuous = model_params.conf_continuous;
else
    % -> if nothing specified, assume continuous confidence on [0, 1]
    conf_continuous = default_continuous;
end

if (any(strcmp(fieldnames(model_params), 'llo_gamma')))
    llo_gamma = model_params.llo_gamma;
else
    llo_gamma = default_gamma;
end

if (any(strcmp(fieldnames(model_params), 'llo_p0')))
    llo_p0 = model_params.llo_p0;
else
    llo_p0 = default_p0;
end

% -> use discrete or continuous confidence judgments
if (conf_continuous)
    % -> confidence range
    if (any(strcmp(fieldnames(model_params), 'conf_cont_range')))
        conf_range = model_params.conf_cont_range;
    else
        conf_range = default_conf_range;
    end
    
    if (any(strcmp(fieldnames(model_params), 'conf_half_scale')))
        conf_half_scale = model_params.conf_half_scale;
    else
        conf_half_scale = default_conf_half_scale;
    end


    % -> approximate continuous with discrete (useful to do faster fit)
    if (any(strcmp(fieldnames(model_params), 'conf_bnds')))
        conf_bnd_lst = model_params.conf_bnds;
        
        % -> reorder confidence boundaries (sometimes out-of-order because
        %    of fminsearch)
        conf_bnd_labels = sort(conf_bnd_lst);
        conf_levels_nb = length(conf_bnd_labels);
        conf_bnd_edges = [conf_range(1), conf_bnd_labels];
    else
        % -> we don't have any other choice than doing a uniform
        %    discretization (not optimal, because some cells can be empty)
        if (any(strcmp(fieldnames(model_params), 'conf_cont_nb_levels')))
            conf_levels_nb = model_params.conf_cont_nb_levels;
        else
            conf_levels_nb = default_conf_levels_nb;   % default resolution
        end
        
        conf_bnd_edges = linspace(conf_range(1), conf_range(2), conf_levels_nb+1);
        conf_bnd_labels = conf_bnd_edges(2:end);
    end

else
    % -> discrete confidence boundaries
    if (any(strcmp(fieldnames(model_params), 'conf_bnds')))
        conf_bnd_lst = model_params.conf_bnds;
        
        % -> reorder confidence boundaries (sometimes out-of-order because
        %    of fminsearch)
        conf_bnd_lst = sort(conf_bnd_lst);
        conf_levels_nb = length(conf_bnd_lst) + 1;

        % -> check confidence boundaries are probabilities
        if ((conf_bnd_lst(end) > 1.0) || (conf_bnd_lst(1) < 0.0))
            fprintf('Error: confidence boundaries are not probabilities\n');
            return;
        end
        conf_bnd_edges = [0, conf_bnd_lst, 1];
    else
        % -> in theory, this should not happen, but if the user "forgot" to
        %    pass the confidence boundaries, spread confidence in 10 linear bins
        conf_levels_nb = default_conf_levels_nb;
        conf_bnd_edges = linspace(0, 1, conf_levels_nb+1);
    end
    conf_bnd_labels = conf_bnd_edges(2:end);
    
end


if (conf_continuous)
    % -> transformation rules
    conf_bnd_min = conf_range(1);
    conf_bnd_max = conf_range(end);

    % -> 1. log-odds <-> confidence ratings
    if (~conf_half_scale)
        conf2lo = @(x_conf) (x_conf - conf_bnd_min) ./ ...
            (conf_bnd_max - conf_bnd_min);
    else
        % -> this is a bizarre formula to fall back on an equivalent full scale
        conf2lo = @(x_conf) (x_conf - 2*conf_bnd_min + conf_bnd_max) ./ ...
            (2*conf_bnd_max - 2*conf_bnd_min);
    end

    % -> 2. confidence prob <-> log-odds
    lo2prob = @(x_lo) cncb_log_odds(x_lo, 1/llo_gamma, llo_p0);
    
    % -> 3. confidence evidence <-> confidence prob
    prob2evid = @(x_prob) norminv(x_prob);
    
    % -> 4. whole transformation: confidence evidence <-> confidence ratings
    conf2evid = @(x_conf) prob2evid(lo2prob(conf2lo(x_conf)));
    
    evid_bnd_edges = conf2evid(conf_bnd_edges);
else
    % -> transformation: confidence evidence <-> confidence prob
    prob2evid = @(x_prob) norminv(x_prob);

    evid_bnd_edges = prob2evid(conf_bnd_edges);
end

% -> adjust boundaries so that we make sure to cover the whole space
evid_bnd_edges(1) = -Inf;
evid_bnd_edges(end) = +Inf;


sens_var_noise = sens_noise^2;
conf_var_noise = conf_noise^2;


rtg_resp0 = NaN(conf_levels_nb, 4);  % stim, resp=0, conf, prob(conf & resp | stim)
rtg_resp1 = NaN(conf_levels_nb, 4);  % stim, resp=1, conf, prob(conf & resp | stim)

nb_conds = stim_nb * 2 * conf_levels_nb;
nb_conds2 = 2 * conf_levels_nb;
cncb_rating_mat = NaN(nb_conds, 4);

for ss = 1:stim_nb
    sens_mean = stim_lst(ss);

    lcl_distA_type1n2_mean = [sens_mean, ...
        (sens_mean - sens_crit - conf_crit)/sens_noise];
    lcl_distA_type1n2_covar = (1 - conf_boost) * sens_noise;
    lcl_distA_type1n2_var2 = (1 - conf_boost)^2 + conf_var_noise;
    lcl_distA_type1n2_var = [sens_var_noise,    lcl_distA_type1n2_covar; ...
                        lcl_distA_type1n2_covar, lcl_distA_type1n2_var2];

    % -> Make sure Sigma is a valid covariance matrix
    % [~, err] = cholcov(lcl_distA_type1n2_var, 0);
    % if err ~= 0
    %     fprintf('Erreur: singular covariance matrix\n');
    %     lcl_distA_type1n2_var(2,2) = lcl_distA_type1n2_var(2,2) + 0.001;
    % end


    for rr = 1:conf_levels_nb
        % -> lower bound
        evid_lower_val = evid_bnd_edges(rr);
        
        % -> upper bound
        evid_upper_val = evid_bnd_edges(rr+1);

        pp_L_val = mvncdf([-Inf, -evid_upper_val], [sens_crit, -evid_lower_val], ...
            lcl_distA_type1n2_mean, lcl_distA_type1n2_var);
        pp_R_val = mvncdf([sens_crit, evid_lower_val], [+Inf, evid_upper_val], ...
            lcl_distA_type1n2_mean, lcl_distA_type1n2_var);

        rtg_resp0(rr, :) = [sens_mean, 0, conf_bnd_labels(rr), pp_L_val];
        rtg_resp1(rr, :) = [sens_mean, 1, conf_bnd_labels(rr), pp_R_val];
    end

    rng = (ss - 1)*nb_conds2 + (1:nb_conds2);
    cncb_rating_mat(rng, :) = [rtg_resp0; rtg_resp1];
end


end

% -> THE END <-
