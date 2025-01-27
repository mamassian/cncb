% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.2
%
% cncb_fit
%   This is the wrap function for confidence model fit.
%
% INPUT:
%   'data_src': grouped data per (stim, resp, conf_id, conf_prob):
%       1st col: stimulus intensity
%       2nd col: perceptual decision
%       3rd col: confidence level (index 1, 2, ...)
%       4th col: confidence probability (for that conf_id | stim, resp)
% note: the 4th column can also just be the count of trials for this condition
%
%
% OPTIONAL PARAMETERS:
%
%   'model_parameters': model free parameters, as a structure
%                       increment numbers for desired parameters
%                       (use NaN for padding):
%       'sens_noise' : sensory (Type 1) sdtev of noise (0 = perfectly sensitive)
%       'sens_crit'  : sensory (Type 1) criterion
%       'conf_noise' : confidence (Type 2) sdtev of noise (0 = ideal)
%       'conf_boost' : fraction super-ideal (0 = ideal, 1 = super-ideal)
%       'conf_crit'  : confidence (Type 2) criterion
%
%   'model_fixed_values': model fixed values parameters, as a structure
%       'sens_noise' : sensory (Type 1) sdtev of noise (0 = perfectly sensitive)
%       'sens_crit'  : sensory (Type 1) criterion
%       'conf_noise' : confidence (Type 2) sdtev of noise (0 = ideal)
%       'conf_boost' : fraction super-ideal (0 = ideal, 1 = super-ideal)
%       'conf_crit'  : confidence (Type 2) criterion
%
%   'boost_init': redo the fit from multiple confidence boost starting points
%       (provide list here as a vector) to avoid local minima
%
%   'noise_init': redo the fit from multiple confidence noise starting points
%       (provide list here as a vector) to avoid local minima
%
%   'skip_efficiency': skip computing efficiency (for faster computations)
%
%   'only_efficiency': compute only efficiency, i.e. skip full model
%                       (useful when number of trials is too small)
%
%   'confidence_half_scale': half scale has its lowest rating for chance,
%                       full scale had its lowest rating for sure error
%                       (default is 'false', i.e. full scale)
%
%   'confidence_is_continuous': are confidence judgments on a continuous scale?
%                (true / false). Default is 'false' (i.e. discrete ratings).
%
%   'confidence_cont_range': if continuous, range of confidence ratings
%                [min, max]. Default is [0, 1].
%
%   'confidence_cont_nb_levels': if continuous, the number of points used to
%                estimate the confidence function. Default is 100.
%
%   'verbose': verbose flag:
%       0: remove all online reports
%       1: print parameter estimates
%       2: print parameter estimates and fitting progress
%
%
%
% OUTPUT: cncb_fit_struct = struct
%   cncb_fit_struct.sens_noise        sdtev of sensory noise (0 = perfect sensitivity)
%   cncb_fit_struct.sens_crit         sensory criterion
%   cncb_fit_struct.conf_noise        sdtev of confidence noise (0 = ideal)
%   cncb_fit_struct.conf_boost        confidence boost (fraction super-ideal)
%   cncb_fit_struct.conf_crit         confidence criterion
%   cncb_fit_struct.intrvl_bias       intrvl_bias (bias in favour of interval 1)
%   cncb_fit_struct.conf_bias         confidence bias (relative to a task set to 1.0)
%   cncb_fit_struct.efficiency        efficiency
%   cncb_fit_struct.choice_prob_ideal            ideal choice probabilities
%   cncb_fit_struct.choice_prob_super_ideal      super-ideal choice probabilities
%   cncb_fit_struct.choice_prob_eff              efficiency choice probabilities
%   cncb_fit_struct.choice_prob_model            full-model choice probabilities
%   cncb_fit_struct.loglike                      log-likelihood of best fit
%
% 
% EXAMPLES OF USE:
%   cncb_fit_struct = cncb_fit(data_src)
%   cncb_fit_struct = cncb_fit(data_src, 'boost2_init', [0.2, 0.5, 0.8])
% 
%
% 23-AUG-2021 - pascal mamassian

% 16-FEB-2022 - pm: added continuous ratings
% 03-NOV-2022 - pm: define ideal observer to match high conf frequency
% 19-OCT-2023 - pm: updated continuous ratings fit
% 11-FEB-2024 - pm: fit cumulative distributions instead of distributions
% 10-MAR-2024 - pm: cleaned up
% 27-APR-2024 - pm: fixed 'cellpadding' & cropped 'conf_cumul'
% 16-JUN-2024 - pm: added confidence boundaries as parameters for continuous
% 17-JUL-2024 - pm: added Deviance statistics
% 22-JUL-2024 - pm: corrected nb trials in likelihood computation
% 22-DEC-2024 - pm: removed fit of cumulative distributions


function cncb_fit_struct = cncb_fit(cncb_data, varargin)

% -> default optional arguments
dflt_model            = struct; % pass on specific model
dflt_fixed            = struct; % pass on fixed values
dflt_boost2_init      = [];     % initial value (or list of values) for boost2
dflt_noise2_init      = [];     % initial value (or list of values) for noise2
dflt_skip_efficiency  = false;  % skip computing efficiency
dflt_only_efficiency  = false;  % compute only efficiency
dflt_is_continuous    = false;  % continuous (T) or discrete (F) rating scale
dflt_conf_cont_half_scale  = false;
dflt_conf_cont_range       = [0, 1];
dflt_conf_cont_nb_levels   = 100;
dflt_verbose          = 0;      % verbose flag


% -> parse all arguments
ip = inputParser;
ip.StructExpand = false;
addRequired(ip, 'data_src', @isnumeric);
addParameter(ip, 'model_parameters', dflt_model, @isstruct);
addParameter(ip, 'model_fixed_values', dflt_fixed, @isstruct);
addParameter(ip, 'boost_init', dflt_boost2_init, @isnumeric);
addParameter(ip, 'noise_init', dflt_noise2_init, @isnumeric);
addParameter(ip, 'skip_efficiency', dflt_skip_efficiency, @islogical);
addParameter(ip, 'only_efficiency', dflt_only_efficiency, @islogical);
addParameter(ip, 'confidence_is_continuous', dflt_is_continuous, @islogical);
addParameter(ip, 'confidence_half_scale', dflt_conf_cont_half_scale, @islogical);
addParameter(ip, 'confidence_cont_range', dflt_conf_cont_range, @isnumeric);
addParameter(ip, 'confidence_cont_nb_levels', dflt_conf_cont_nb_levels, @isnumeric);
addParameter(ip, 'verbose', dflt_verbose, @isnumeric);


parse(ip, cncb_data, varargin{:});
my_model_params = ip.Results.model_parameters;
my_fixed_params = ip.Results.model_fixed_values;
boost_init_lst  = ip.Results.boost_init;
noise_init_lst  = ip.Results.noise_init;
skip_efficiency = ip.Results.skip_efficiency;
only_efficiency = ip.Results.only_efficiency;
conf_continuous = ip.Results.confidence_is_continuous;
conf_half_scale = ip.Results.confidence_half_scale;
conf_cont_range = ip.Results.confidence_cont_range;
conf_cont_nb_levels = ip.Results.confidence_cont_nb_levels;
verbose_flag    = ip.Results.verbose;



options1 = optimset;
% options2 = optimset('TolFun',1e-3, 'TolX',1e-3);
options2 = optimset;
if (verbose_flag >= 2)
    options3 = optimset('Display','iter');
    tstart = tic;   % start timer
else
    options3 = optimset;
end



% -> 0. extract proportion of choices for each confidence level
%    for the original human data

% -> clean up input data
if (conf_continuous)
    src_data_count = cncb_group(cncb_data, ...
        'confidence_is_continuous', conf_continuous, ...
        'confidence_cont_range', conf_cont_range, ...
        'confidence_half_scale', conf_half_scale, ...
        'confidence_cont_nb_levels', conf_cont_nb_levels);
else
    src_data_count = cncb_group(cncb_data, ...
        'confidence_is_continuous', conf_continuous);
end

col_stim = 1;
col_resp = 2;
col_conf_levl = 3;
col_conf_prob = 4;

[stim_lst, ~, stim_ic] = unique(cncb_data(:, col_stim));
stim_nb = length(stim_lst);

% -> count how many confidence judgments were made for each stimulus
%    (useful to compute likelihood)
conf_count_stim = NaN(stim_nb, 1);
for ss = 1:stim_nb
    cond_inds = (stim_ic == ss);
    conf_count_stim(ss) = sum(cncb_data(cond_inds, col_conf_prob));
end


% -> count how many confidence judgments were made for each confidence level
[conf_lst, ~, conf_ic] = unique(cncb_data(:, col_conf_levl));
conf_nb = length(conf_lst);
conf_count_conf = count_conf(cncb_data, conf_nb, conf_ic);


% -> extended stimulus list, with number of trials per stimulus
stim_lst_extended = [stim_lst, conf_count_stim];

resp_lst = unique(cncb_data(:, col_resp));
resp_nb = length(resp_lst);

if (size(cncb_data,1) == (stim_nb*resp_nb))
    % -> it looks like all the data have the same confidence
    fprintf('CNCB fit error: all the data have the same confidence\n');

    cncb_fit_struct = -1;   % return error
    return;
end

src_data_norm = src_data_count;
data_SRC_norm = NaN(size(src_data_norm,1),1);
data_SRC_count = NaN(size(src_data_norm,1),1);

rating_inds = unique(src_data_norm(:, col_conf_levl));
conf_levels_nb = length(rating_inds);
rating_bnds_nb = conf_levels_nb - 1;

% -> convert counts of confidence ratings into probabilities 
%    (conditional on stimulus)
stim_prob = NaN(1, stim_nb);    % prob. occurence of each stimulus
for ss = 1:stim_nb
    sens_mean = stim_lst(ss);
    cond_inds = (src_data_norm(:, col_stim) == sens_mean);
    stim_prob(ss) = sum(src_data_norm(cond_inds, col_conf_prob));
    data_SRC_norm(cond_inds) = ...
        src_data_norm(cond_inds, col_conf_prob) ./ stim_prob(ss);
end
src_data_norm(:, col_conf_prob) = data_SRC_norm;

% -> use cell padding technique advised in 'type2_SDT_MLE'
% cellpadding = 1 / (10*conf_levels_nb);
cellpadding = 1 / (100000*conf_levels_nb);
if any(src_data_norm(:,col_conf_prob) < cellpadding)
    inds = find(src_data_norm(:,col_conf_prob) < cellpadding);
    src_data_norm(inds, col_conf_prob) = src_data_norm(inds, col_conf_prob) + ...
        cellpadding;
end

% -> fix the normalization so that cumulative prob is 1
stim_prob = NaN(1, stim_nb);    % prob. occurence of each stimulus
for ss = 1:stim_nb
    sens_mean = stim_lst(ss);
    cond_inds = (src_data_norm(:, col_stim) == sens_mean);
    stim_prob(ss) = sum(src_data_norm(cond_inds, col_conf_prob));
    data_SRC_norm(cond_inds) = ...
        src_data_norm(cond_inds, col_conf_prob) ./ stim_prob(ss);
end
src_data_norm(:, col_conf_prob) = data_SRC_norm;

if (conf_continuous)
    conf_bnd_labels = rating_inds';
    conf_levels_nb = length(conf_bnd_labels) + 1;
    rating_bnds_nb = conf_levels_nb - 1;
end

% -> normalize to get sum(prob | stim)=1
for ss = 1:stim_nb
    sens_mean = stim_lst(ss);
    cond_inds = (src_data_norm(:, col_stim) == sens_mean);
    data_SRC_norm(cond_inds) = ...
        src_data_norm(cond_inds, col_conf_prob) ./ ...
        sum(src_data_norm(cond_inds, col_conf_prob));
    data_SRC_count(cond_inds) = data_SRC_norm(cond_inds) .* conf_count_stim(ss);
end
src_data_count(:, col_conf_prob) = data_SRC_count;


% -> compute confidence sums over (stim)
conf_sum_stim = NaN(1, stim_nb);
for cnd = 1:stim_nb
    conf_inds5 = (stim_ic == cnd);
    cum_tmp5 = sum(src_data_count(conf_inds5, 4));
    conf_sum_stim(cnd) = cum_tmp5;
end

% -> human data to fit
human_tofit = lastcol(src_data_count);


% -> default values of the parameters
default_sens_noise  = 1.0;
default_sens_crit   = 0.0;
default_conf_noise  = 0.001;    % do not set to 0 to avoid singularity
default_conf_boost  = 0.0;
default_conf_crit   = 0.0;
if (~conf_continuous)
    bnds = linspace(0.5, 1.0, rating_bnds_nb+2);
    default_conf_bnds = bnds(2:(end-1));
else
    default_llo_gamma   = 1.0;
    default_llo_p0      = 0.5;

    bnds = linspace(0.5, 1.0, rating_bnds_nb+2) * conf_cont_range(2);
    default_conf_bnds = bnds(2:(end-1));
end

% -> initial values of the parameters
initial_sens_noise  = 1.2;
initial_sens_crit   = 0.1;
initial_conf_noise  = 0.5;
initial_conf_boost  = 0.5;
initial_conf_crit   = 0.4;
if (~conf_continuous)
    % -> avoid boundaries below 0.5 because ideally there should be no data
    initial_conf_bnds   = linspace(0.4, 0.9, rating_bnds_nb);

else
    initial_llo_gamma   = 0.8;
    initial_llo_p0      = 0.4;

    initial_conf_bnds   = linspace(0.4, 0.9, rating_bnds_nb) * conf_cont_range(2);
end


% -> lower and upper bound values of the parameters
lo_bnd_sens_noise   = 0.0;      hi_bnd_sens_noise   = Inf;
lo_bnd_sens_crit    = -Inf;     hi_bnd_sens_crit    = Inf;
lo_bnd_conf_noise   = 0.0;      hi_bnd_conf_noise   = Inf;
lo_bnd_conf_boost   = 0.0;      hi_bnd_conf_boost   = 1.0;
lo_bnd_conf_crit    = -Inf;     hi_bnd_conf_crit    = Inf;
if (~conf_continuous)
    lo_bnd_conf_bnds = zeros(1, rating_bnds_nb);
    hi_bnd_conf_bnds =  ones(1, rating_bnds_nb);
else
    lo_bnd_llo_gamma    = 0.0;     hi_bnd_llo_gamma     = Inf;
    lo_bnd_llo_p0       = 1e-6;     hi_bnd_llo_p0       = 1-1e-6;

    lo_bnd_conf_bnds = zeros(1, rating_bnds_nb);
    hi_bnd_conf_bnds =  ones(1, rating_bnds_nb) * conf_cont_range(2);
end


% -> initial values of the parameters
initial_params = struct;
initial_params.sens_noise = initial_sens_noise;
initial_params.sens_crit  = initial_sens_crit;
initial_params.conf_noise = initial_conf_noise;
initial_params.conf_boost = initial_conf_boost;
initial_params.conf_crit  = initial_conf_crit;
if (~conf_continuous)
    initial_params.conf_bnds = initial_conf_bnds;
else
    initial_params.llo_gamma  = initial_llo_gamma;
    initial_params.llo_p0     = initial_llo_p0;
    initial_params.conf_bnds = initial_conf_bnds;
end

% -> lower and upper bound values of the parameters
lo_bnd_params = struct;
lo_bnd_params.sens_noise = lo_bnd_sens_noise;
lo_bnd_params.sens_crit  = lo_bnd_sens_crit;
lo_bnd_params.conf_noise = lo_bnd_conf_noise;
lo_bnd_params.conf_boost = lo_bnd_conf_boost;
lo_bnd_params.conf_crit  = lo_bnd_conf_crit;
if (~conf_continuous)
    lo_bnd_params.conf_bnds = lo_bnd_conf_bnds;
else
    lo_bnd_params.llo_gamma  = lo_bnd_llo_gamma;
    lo_bnd_params.llo_p0     = lo_bnd_llo_p0;
    lo_bnd_params.conf_bnds = lo_bnd_conf_bnds;
end

hi_bnd_params = struct;
hi_bnd_params.sens_noise = hi_bnd_sens_noise;
hi_bnd_params.sens_crit  = hi_bnd_sens_crit;
hi_bnd_params.conf_noise = hi_bnd_conf_noise;
hi_bnd_params.conf_boost = hi_bnd_conf_boost;
hi_bnd_params.conf_crit  = hi_bnd_conf_crit;
if (~conf_continuous)
    hi_bnd_params.conf_bnds = hi_bnd_conf_bnds;
else
    hi_bnd_params.llo_gamma  = hi_bnd_llo_gamma;
    hi_bnd_params.llo_p0     = hi_bnd_llo_p0;
    hi_bnd_params.conf_bnds = hi_bnd_conf_bnds;
end


% -> 'params_set' is a struct that defines the free parameters, and 
%    their order. Put '0' if the variable is not a free parameter.
default_params_set = struct;
default_params_set.sens_noise = 0;
default_params_set.sens_crit  = 0;
default_params_set.conf_noise = 0;
default_params_set.conf_boost = 0;
default_params_set.conf_crit  = 0;  % by default, we don't fit this parameter
if (~conf_continuous)
    default_params_set.conf_bnds = 0;
else
    default_params_set.llo_gamma  = 0;
    default_params_set.llo_p0     = 0;
    default_params_set.conf_bnds = 0;
end

% -> default values of the parameters
default_params_val = struct;
default_params_val.sens_noise = default_sens_noise;
default_params_val.sens_crit  = default_sens_crit;
default_params_val.conf_noise = default_conf_noise;
default_params_val.conf_boost = default_conf_boost;
default_params_val.conf_crit  = default_conf_crit;
if (~conf_continuous)
    default_params_val.conf_bnds = default_conf_bnds;
else
    default_params_val.llo_gamma  = default_llo_gamma;
    default_params_val.llo_p0     = default_llo_p0;
    default_params_val.conf_bnds = default_conf_bnds;
end


% -> check how many parameters we should fit
my_params_cell = struct2cell(my_model_params);
my_params_mat = cell2mat(my_params_cell');
param_free_nb = max(my_params_mat);

% -> check which parameters we should fit
params_set = default_params_set;    % reset params_set
params_set1 = default_params_set;   % params_set for only Type 1
params_set2 = default_params_set;   % params_set for only Type 2
param_free_nb1 = 0;     % nb free parameters for Type 1
param_free_nb2 = 0;     % nb free parameters for Type 2

my_fld_nms = fieldnames(my_model_params);
for my_kk = 1:param_free_nb
    my_fld_ind = cellfun(@(xx) find(xx == my_kk), my_params_cell, 'UniformOutput', false);
    my_uu = find(~cellfun(@isempty, my_fld_ind));

    % -> allow for multiple use of a variable across model parameters
    for my_pp = 1:length(my_uu)
        my_vv = my_fld_ind{my_uu(my_pp)};

        % -> allow for multiple use of a variable within model parameters
        if (strcmp('sens_noise', my_fld_nms{my_uu(my_pp)}) || ...
                    strcmp('sens_crit', my_fld_nms{my_uu(my_pp)}))
            param_free_nb1 = param_free_nb1 + 1;
        else
            param_free_nb2 = param_free_nb2 + 1;
        end

        for my_qq = 1:length(my_vv)
            params_set.(my_fld_nms{my_uu(my_pp)})(my_vv(my_qq)) = my_kk;
            if (strcmp('sens_noise', my_fld_nms{my_uu(my_pp)}) || ...
                    strcmp('sens_crit', my_fld_nms{my_uu(my_pp)}))
                params_set1.(my_fld_nms{my_uu(my_pp)})(my_vv(my_qq)) = ...
                    my_kk;
            else
                params_set2.(my_fld_nms{my_uu(my_pp)})(my_vv(my_qq)) = ...
                    my_kk - param_free_nb1;
            end
        end
    end
end


% -> if (part of) a model was set, pick the corresponding parameters
% -> check which values are fixed
fixed_set = default_params_set;

% -> use the best fit of Type 1 to estimate Type 2 parameters
fixed_vals = default_params_val;

my_params_cell = struct2cell(my_fixed_params);

my_fld_nms = fieldnames(my_fixed_params);
my_fld_ind = cellfun(@(xx) ~isnan(xx), my_params_cell, 'UniformOutput', false);
my_uu = find(~cellfun(@isempty, my_fld_ind));

% -> allow for multiple use of a variable across model parameters
for my_pp = 1:length(my_uu)
    my_ww = find(my_fld_ind{my_uu(my_pp)});

    % -> allow for multiple use of a variable within model parameters
    for my_qq = my_ww

        % -> remember that this variable was set by user
        fixed_set.(my_fld_nms{my_uu(my_pp)})(my_qq) = 1;

        % -> remember the value of this variable
        fixed_vals.(my_fld_nms{my_uu(my_pp)})(my_qq) = ...
            my_fixed_params.(my_fld_nms{my_uu(my_pp)})(my_qq);
    end
end


% -> params for ideal confidence observer
ideal_conf_boost = 0.0;
ideal_conf_noise = 0.001;

% -> params for super-ideal confidence observer
superideal_conf_boost = 1.0;


% ------------------------------------------------------------------------
% -> 1. 'type1': fit Type 1 model (cumulative Gaussian)
%       to extract the Type 1 parameters, i.e. (sens_noise) and (sens_crit)
resp_vec = NaN(stim_nb, 1);
for ss = 1:stim_nb
    stim_val = stim_lst(ss);
    lcl_righ_inds = ((src_data_count(:,col_stim) == stim_val) & ...
                     (src_data_count(:,col_resp) == resp_lst(2)));
    lcl_all_inds = (src_data_count(:,col_stim) == stim_val);
    resp_vec(ss) = sum(src_data_count(lcl_righ_inds, col_conf_prob)) / ...
                   sum(src_data_count(lcl_all_inds, col_conf_prob));
end
nn1_resp_list = [resp_vec, 1-resp_vec]; 
nn1_resp_list = nn1_resp_list .* conf_sum_stim';

% -> check that Type 1 is indeed fitted
params_set_type1 = struct;
nb_params_to_fit = 1;
if (isfield(my_model_params, 'sens_noise'))
    params_set_type1.sens_noise = my_model_params.sens_noise;
else
    params_set_type1.sens_noise = nb_params_to_fit;
end
if (params_set_type1.sens_noise == 1)
    nb_params_to_fit = nb_params_to_fit + 1;
end
if (isfield(my_model_params, 'sens_crit'))
    params_set_type1.sens_crit = my_model_params.sens_crit;
else
    params_set_type1.sens_crit = nb_params_to_fit;
end

% -> fit cumulative Gaussian
params_set_cumul = params_set_type1;
params0_cumul  = extract_params_from_struct(params_set_cumul, initial_params);
paramsLB_cumul = extract_params_from_struct(params_set_cumul, lo_bnd_params);
paramsUB_cumul = extract_params_from_struct(params_set_cumul, hi_bnd_params);

my_fun_0 = @(pp) cfc_type1(stim_lst, pp, params_set_type1, fixed_vals);

param_type1 = fitnllbin(my_fun_0, nn1_resp_list, ...
    params0_cumul, paramsLB_cumul, paramsUB_cumul, options1);


% -> store fitted params for step 2
model_params_type2 = pack_params_in_struct(param_type1, params_set_type1, fixed_vals);

idl_sens_noise = model_params_type2.sens_noise;
idl_sens_crit  = model_params_type2.sens_crit;


% ------------------------------------------------------------------------
% -> 2. 'ideal': build the ideal observer that is constrained by Type 1 
%       performance (inferred from data), 
%       that has (conf_noise = 0) & (conf_boost = 0)
%       and allow its confidence ratings to best match the human observer

if (~skip_efficiency)
    
model_idl_params = struct;
model_idl_params.sens_noise = idl_sens_noise;
model_idl_params.sens_crit  = idl_sens_crit;
model_idl_params.conf_noise = ideal_conf_noise;
model_idl_params.conf_boost = ideal_conf_boost;
model_idl_params.conf_crit  = default_conf_crit;
model_idl_params.conf_continuous  = conf_continuous;

if (conf_continuous)
    
    % -> compute the confidence ratings for the ideal confidence observer
    model_idl_params.conf_cont_range  = conf_cont_range;
    model_idl_params.conf_half_scale  = conf_half_scale;
    model_idl_params.conf_cont_nb_levels = conf_cont_nb_levels;

    % -> tries to match (gamma, p0) to confidence data
    fixed_vals_idl = struct;
    fixed_vals_idl.sens_noise = idl_sens_noise;
    fixed_vals_idl.sens_crit  = idl_sens_crit;
    fixed_vals_idl.conf_noise = ideal_conf_noise;
    fixed_vals_idl.conf_boost = ideal_conf_boost;
    fixed_vals_idl.conf_crit  = default_conf_crit;
    fixed_vals_idl.conf_continuous = conf_continuous;
    fixed_vals_idl.conf_cont_range = conf_cont_range;
    fixed_vals_idl.conf_half_scale = conf_half_scale;
    fixed_vals_idl.conf_cont_nb_levels = conf_cont_nb_levels;

    params_set_idl = default_params_set;

    nb_params_to_fit = 1;
    if (isfield(my_model_params, 'llo_gamma'))
        params_set_idl.llo_gamma = params_set2.llo_gamma;
    else
        params_set_idl.llo_gamma = nb_params_to_fit;
    end
    if (params_set_idl.llo_gamma == nb_params_to_fit)
        nb_params_to_fit = nb_params_to_fit + 1;
    else
        fixed_vals_idl.llo_gamma = fixed_vals.llo_gamma;
    end

    if (isfield(my_model_params, 'llo_p0'))
        params_set_idl.llo_p0 = params_set2.llo_p0;
    else
        params_set_idl.llo_p0 = nb_params_to_fit;
    end
    if (params_set_idl.llo_p0 == nb_params_to_fit)
        nb_params_to_fit = nb_params_to_fit + 1; 
    else
        fixed_vals_idl.llo_p0 = fixed_vals.llo_p0;
    end

    % -> try to adjust the confidence boundaries because here we are 
    %    interested in fitting a curve to the Type 2 ROC, not the
    %    exact points where the human observer set her boundaries.
    params_set_idl.conf_bnds = nb_params_to_fit - 1 + (1:rating_bnds_nb);

    params0_idl  = extract_params_from_struct(params_set_idl, initial_params);
    paramsLB_idl = extract_params_from_struct(params_set_idl, lo_bnd_params);
    paramsUB_idl = extract_params_from_struct(params_set_idl, hi_bnd_params);

    % -> extract (gamma, p0) for the ideal confidence observer
    %    by approximating the confidence frequency in each level
    my_fun_ideal3 = @(pp) count_conf(cncb_core_count_wrap(stim_lst_extended, ...
        pp, params_set_idl, fixed_vals_idl), conf_nb, conf_ic);

    [best_params_ideal, ~] = ...
        fitnll(my_fun_ideal3, conf_count_conf, ...
        params0_idl, paramsLB_idl, paramsUB_idl, options1);


    paramBest_struct = pack_params_in_struct(best_params_ideal, params_set_idl, fixed_vals_idl);
    model_idl_params.llo_gamma = paramBest_struct.llo_gamma;
    model_idl_params.llo_p0 = paramBest_struct.llo_p0;

    model_idl_params.conf_bnds = sort(paramBest_struct.conf_bnds);

else
    params_set_idl = default_params_set;
    params_set_idl.conf_bnds  = (1:rating_bnds_nb);

    fixed_vals_idl = struct;
    fixed_vals_idl.sens_noise = idl_sens_noise;
    fixed_vals_idl.sens_crit  = idl_sens_crit;
    fixed_vals_idl.conf_noise = ideal_conf_noise;
    fixed_vals_idl.conf_boost = ideal_conf_boost;
    fixed_vals_idl.conf_crit  = default_conf_crit;
    fixed_vals_idl.conf_continuous = conf_continuous;

    params0_idl = initial_conf_bnds;
    paramsLB_idl = lo_bnd_conf_bnds;
    paramsUB_idl = hi_bnd_conf_bnds;
    

    % -> extract the confidence boundaries for the ideal confidence observer
    %    by adjusting the confidence frequency in each level
    my_fun_ideal3 = @(pp) count_conf(cncb_core_count_wrap(stim_lst_extended, ...
        pp, params_set_idl, fixed_vals_idl), conf_nb, conf_ic); 
    
    [best_params_ideal, ~] = ...
        fitnll(my_fun_ideal3, conf_count_conf, ...
        params0_idl, paramsLB_idl, paramsUB_idl, options1);
    
    % -> boundaries are sometimes in random order
    bnd_lst_ideal = sort(best_params_ideal); 
    model_idl_params.conf_bnds  = bnd_lst_ideal;

end

% -> compute the confidence ratings for the ideal confidence observer
ideal_rating_SRC = cncb_core(stim_lst, model_idl_params);

ideal_rating_count = cncb_core_count(stim_lst_extended, model_idl_params);
ideal_tofit = lastcol(ideal_rating_count);


% ------------------------------------------------------------------------
% -> 3. 'super-human': build the super-ideal confidence observer
%    that has (conf_boost = 1) and free (conf_noise)
%    and extract equivalent confidence noise for this 'super-human'
params_set_super = default_params_set;
params_set_super.conf_noise = 1;

fixed_vals_super = struct;
fixed_vals_super.sens_noise = idl_sens_noise;
fixed_vals_super.sens_crit  = idl_sens_crit;
fixed_vals_super.conf_crit  = default_conf_crit;
fixed_vals_super.conf_boost = superideal_conf_boost;
fixed_vals_super.conf_continuous  = conf_continuous;

params0_super = initial_conf_noise;
paramsLB_super = lo_bnd_conf_noise;
paramsUB_super = hi_bnd_conf_noise;


% -> do fit
if (conf_continuous)
    % -> extract the log-odds parameters for the 'super-human'
    
    fixed_vals_super.conf_cont_range  = conf_cont_range;
    fixed_vals_super.conf_half_scale  = conf_half_scale;
    fixed_vals_super.conf_cont_nb_levels = conf_cont_nb_levels;
    
    params_set_super.llo_gamma = 2;
    params_set_super.llo_p0 = 3;
    params_set_super.conf_bnds = 3 + (1:rating_bnds_nb);

    params0_super = [params0_super, initial_llo_gamma, initial_llo_p0, initial_conf_bnds];
    paramsLB_super = [paramsLB_super, lo_bnd_llo_gamma, lo_bnd_llo_p0, lo_bnd_conf_bnds];
    paramsUB_super = [paramsUB_super, hi_bnd_llo_gamma, hi_bnd_llo_p0, hi_bnd_conf_bnds];

    my_fun_super = @(pp) lastcol(cncb_core_wrap(stim_lst, ...
        pp, params_set_super, fixed_vals_super)); 

    [super_human_params, super_human_loglike] = fitnll(my_fun_super, human_tofit, ...
        params0_super, paramsLB_super, paramsUB_super, options2);

else
    % -> extract the confidence boundaries for the 'super-human'
    params_set_super.conf_bnds = 1 + (1:rating_bnds_nb);

    params0_super = [params0_super, initial_conf_bnds];
    paramsLB_super = [paramsLB_super, lo_bnd_conf_bnds];
    paramsUB_super = [paramsUB_super, hi_bnd_conf_bnds];
            
    my_fun_super = @(pp) lastcol(cncb_core_wrap(stim_lst, ...
        pp, params_set_super, fixed_vals_super)); 

    [super_human_params, super_human_loglike] = fitnll(my_fun_super, human_tofit, ...
        params0_super, paramsLB_super, paramsUB_super, options2);
        
end

super_human_df = length(super_human_params);

equiv_noise_hum = super_human_params(1);


% -> compute the confidence ratings for the 'super-human'
params_super_human = struct;
params_super_human.sens_noise = idl_sens_noise;
params_super_human.sens_crit  = idl_sens_crit;
params_super_human.conf_noise = equiv_noise_hum;
params_super_human.conf_boost = superideal_conf_boost;
params_super_human.conf_crit  = default_conf_crit;
params_super_human.conf_continuous  = conf_continuous;

if (conf_continuous)
    llo_gamma_super_human = super_human_params(2);
    llo_p0_super_human = super_human_params(3);
    
    params_super_human.conf_cont_range  = conf_cont_range;
    params_super_human.conf_half_scale  = conf_half_scale;
    params_super_human.conf_cont_nb_levels = conf_cont_nb_levels;
    params_super_human.llo_gamma = llo_gamma_super_human;
    params_super_human.llo_p0 = llo_p0_super_human;
    params_super_human.conf_bnds = sort(super_human_params(4:end));

else
    % -> sort confidence boundaries in increasing order
    super_human_rtng_bnd_lst = sort(super_human_params(2:end));

    % -> compute the confidence ratings for the 'super-human'
    params_super_human.conf_bnds  = super_human_rtng_bnd_lst;

end

super_human_rating_SRC = cncb_core(stim_lst, params_super_human);


% ------------------------------------------------------------------------
% -> 4. 'super-ideal': get the equivalent noise for the 'ideal' super-ideal 
%    confidence observer

[super_ideal_params, super_ideal_loglike] = fitnll(my_fun_super, ...
    ideal_tofit, params0_super, paramsLB_super, paramsUB_super, ...
    options2);   

super_ideal_df = length(super_ideal_params);

equiv_noise_idl = super_ideal_params(1);


% -> compute the confidence ratings for the 'super-ideal'
params_super_ideal = struct;
params_super_ideal.sens_noise = idl_sens_noise;
params_super_ideal.sens_crit  = idl_sens_crit;
params_super_ideal.conf_noise = equiv_noise_idl;
params_super_ideal.conf_boost = superideal_conf_boost;
params_super_ideal.conf_crit  = default_conf_crit;
params_super_ideal.conf_continuous  = conf_continuous;

if (conf_continuous)
    params_super_ideal.conf_cont_range  = conf_cont_range;
    params_super_ideal.conf_half_scale  = conf_half_scale;
    params_super_ideal.conf_cont_nb_levels = conf_cont_nb_levels;
    params_super_ideal.conf_bnds = sort(super_ideal_params(4:end));

    params_super_ideal.llo_gamma = super_ideal_params(2);
    params_super_ideal.llo_p0 = super_ideal_params(3);
else
    super_ideal_rtng_bnd_lst = sort(super_ideal_params(2:end));
    params_super_ideal.conf_bnds  = super_ideal_rtng_bnd_lst;
end

super_ideal_rating_SRC = cncb_core(stim_lst, params_super_ideal);


% -> ********
% -> use ratio of variance (std dev squared) for definition of efficiency
efficiency = (equiv_noise_idl / equiv_noise_hum)^2;

% -> compute Deviance for efficiency
super_human_fun_sat = @(pp) lastcol(src_data_norm); 
super_human_loglike_sat = - loglikefcn([], super_human_fun_sat, human_tofit);
super_human_df_sat = size(human_tofit, 1);

super_ideal_fun_sat = @(pp) lastcol(ideal_rating_SRC); 
super_ideal_loglike_sat = - loglikefcn([], super_ideal_fun_sat, ideal_tofit);
super_ideal_df_sat = size(ideal_tofit, 1);

eff_model_loglike = super_human_loglike + super_ideal_loglike;
eff_sat_loglike = super_human_loglike_sat + super_ideal_loglike_sat;

eff_G2 = -2 * (eff_model_loglike - eff_sat_loglike);

eff_model_df = super_human_df + super_ideal_df;
eff_sat_df = super_human_df_sat + super_ideal_df_sat;

eff_df_chi2 = eff_sat_df - eff_model_df;
eff_p_val = chi2cdf(eff_G2, eff_model_df, 'upper');   % if p>0.01, good fit


end     % ~skip_efficiency


% ------------------------------------------------------------------------
% -> 5. 'full-model'

if (~only_efficiency)

params_set_full = struct;

% -> keep the Type 1 parameters from the ideal fit
fixed_vals_full = struct;
fixed_vals_full.sens_noise = idl_sens_noise;
fixed_vals_full.sens_crit  = idl_sens_crit;
fixed_vals_full.conf_continuous  = conf_continuous;


% -> check that Type 2 parameters are indeed fitted
nb_params_to_fit = 1;
if (isfield(my_model_params, 'conf_noise'))
    params_set_full.conf_noise = params_set2.conf_noise;
else
    params_set_full.conf_noise = nb_params_to_fit;
end
if (params_set_full.conf_noise == nb_params_to_fit)
    param_noise_ind = nb_params_to_fit;
    nb_params_to_fit = nb_params_to_fit + 1;
else
    fixed_vals_full.conf_noise = fixed_vals.conf_noise;
end

if (isfield(my_model_params, 'conf_boost'))
    params_set_full.conf_boost = params_set2.conf_boost;
else
    params_set_full.conf_boost = nb_params_to_fit;
end
if (params_set_full.conf_boost == nb_params_to_fit)
    param_boost_ind = nb_params_to_fit;
    nb_params_to_fit = nb_params_to_fit + 1;
else
    fixed_vals_full.conf_boost = fixed_vals.conf_boost;
end

if (isfield(my_model_params, 'conf_crit'))
    params_set_full.conf_crit = params_set2.conf_crit;
else
    params_set_full.conf_crit = 0;  % by default, do not fit this param
end
if (params_set_full.conf_crit == nb_params_to_fit)
    nb_params_to_fit = nb_params_to_fit + 1;
else
    fixed_vals_full.conf_crit = fixed_vals.conf_crit;
end


if (conf_continuous)
    
    fixed_vals_full.conf_cont_range  = conf_cont_range;
    fixed_vals_full.conf_half_scale  = conf_half_scale;
    fixed_vals_full.conf_cont_nb_levels = conf_cont_nb_levels;

    % -> do not fit the confidence boundaries since there seems
    %    to be a trade-off with confidence noise and boost
    fixed_vals_full.conf_bnds = conf_bnd_labels;

    if (isfield(my_model_params, 'llo_gamma'))
        params_set_full.llo_gamma = params_set2.llo_gamma;
    else
        params_set_full.llo_gamma = nb_params_to_fit;
    end
    if (params_set_full.llo_gamma == nb_params_to_fit)
        nb_params_to_fit = nb_params_to_fit + 1;
    else
        fixed_vals_full.llo_gamma = fixed_vals.llo_gamma;
    end

    if (isfield(my_model_params, 'llo_p0'))
        params_set_full.llo_p0 = params_set2.llo_p0;
    else
        params_set_full.llo_p0 = nb_params_to_fit;
    end
    if (params_set_full.llo_p0 == nb_params_to_fit)
        nb_params_to_fit = nb_params_to_fit + 1;
    else
        fixed_vals_full.llo_p0 = fixed_vals.llo_p0;
    end


else

    if (isfield(my_model_params, 'conf_bnds'))
        params_set_full.conf_bnds = params_set2.conf_bnds;
    else
        param_bnds_ind = nb_params_to_fit:(rating_bnds_nb + nb_params_to_fit - 1);
        params_set_full.conf_bnds  = param_bnds_ind;
    end
    if (params_set_full.conf_bnds(1) == nb_params_to_fit)
        nb_params_to_fit = nb_params_to_fit + rating_bnds_nb;
    else
        fixed_vals_full.conf_bnds = fixed_vals.conf_bnds;
    end

end
param_free_nb2 = nb_params_to_fit - 1;

params0_full  = extract_params_from_struct(params_set_full, initial_params);
paramsLB_full = extract_params_from_struct(params_set_full, lo_bnd_params);
paramsUB_full = extract_params_from_struct(params_set_full, hi_bnd_params);


% -> redo the fit from multiple boost2 starting points to avoid local minima
if (params_set_full.conf_boost)
    if (isempty(boost_init_lst))
        boost_init_lst = initial_conf_boost;
    end
    boost2_init_nb = length(boost_init_lst);
else
    boost2_init_nb = 1; % do one round of fit
end

% -> redo the fit from multiple noise2 starting points to avoid local minima
if (params_set_full.conf_noise)
    if (isempty(noise_init_lst))
        noise_init_lst = initial_conf_noise;
    end
    noise2_init_nb = length(noise_init_lst);
else
    noise2_init_nb = 1; % do one round of fit
end

paramBest_mat = NaN(boost2_init_nb, noise2_init_nb, param_free_nb2);
loglike_lst = NaN(boost2_init_nb, noise2_init_nb, 1);


% -> do the fit for multiple starting values of 'boost2'
for bb = 1:boost2_init_nb
for nn = 1:noise2_init_nb
    
    if (params_set_full.conf_boost)
        boost2_init_val = boost_init_lst(bb);
        params0_full(param_boost_ind) = boost2_init_val;
    end
    if (params_set_full.conf_noise)
        noise2_init_val = noise_init_lst(nn);
        params0_full(param_noise_ind) = noise2_init_val;
    end
    

    my_fun_full = @(pp) lastcol(cncb_core_wrap(stim_lst, ...
        pp, params_set_full, fixed_vals_full));

    [paramBest, loglike] = fitnll(my_fun_full, human_tofit, ...
        params0_full, paramsLB_full, paramsUB_full, options3);


    paramBest_mat(bb, nn, :) = paramBest;
    loglike_lst(bb, nn) = loglike;

    if (verbose_flag >= 1)
        if ((boost2_init_nb > 1) && (noise2_init_nb > 1))
            fprintf(['Init (noise2, boost2) = (%5.3f, %5.3f) > ', ...
                '(noise2, boost2) = (%5.3f, %5.3f), loglike = %7.4f\n'], ...
                noise2_init_val, boost2_init_val, ...
                paramBest(1), paramBest(2), loglike);
        elseif (noise2_init_nb > 1)
            fprintf('Init noise2 = %5.3f > noise2 = %5.3f, loglike = %7.4f\n', ...
                noise2_init_val, paramBest(1), loglike);
        elseif (boost2_init_nb > 1)
            fprintf('Init boost2 = %5.3f > boost2 = %5.3f, loglike = %7.4f\n', ...
                boost2_init_val, paramBest(1), loglike);
        end
    end
end
end

% -> pick the initial boost value that led to maximum likelihood
if ((boost2_init_nb > 1) && (noise2_init_nb > 1))
    [jj1, ii1] = max(loglike_lst);
    ii1 = ii1(1);
    [~, ii2] = max(jj1);
elseif (noise2_init_nb > 1)
    [~, ii2] = max(loglike_lst);
    ii1 = 1;
elseif (boost2_init_nb > 1)
    [~, ii1] = max(loglike_lst);
    ii2 = 1;
else
    ii1 = 1;
    ii2 = 1;
end
model_full_params = pack_params_in_struct(paramBest_mat(ii1, ii2, :), params_set_full, fixed_vals_full);
model_full_loglike = loglike_lst(ii1, ii2);
if (verbose_flag >= 1)
    if ((boost2_init_nb > 1) && (noise2_init_nb > 1))
        fprintf(['> choosing fit #(%d, %d): conf_noise = %5.3f, ', ...
            'conf_boost = %5.3f, loglike = %7.3f\n'], ... 
            ii1, ii2, model_full_params.conf_noise, ...
            model_full_params.conf_boost, model_full_loglike);
    elseif (noise2_init_nb > 1)
        fprintf('> choosing fit #%d: conf_noise = %5.3f, loglike = %7.3f\n', ... 
            ii2, model_full_params.conf_noise, model_full_loglike);
    elseif (boost2_init_nb > 1)
        fprintf('> choosing fit #%d: conf_boost = %5.3f, loglike = %7.3f\n', ... 
            ii1, model_full_params.conf_boost, model_full_loglike);
   end
end

% -> extract the confidence parameters for the 'full-model'
full_conf_noise = model_full_params.conf_noise;
full_conf_boost = model_full_params.conf_boost;
full_conf_crit  = model_full_params.conf_crit;

full_rating_SRC = cncb_core(stim_lst, model_full_params);


% -> compute Deviance
my_fun_sat = @(pp) lastcol(src_data_norm); 
loglike_sat = - loglikefcn([], my_fun_sat, human_tofit);
G2 = -2 * (model_full_loglike - loglike_sat);
df_sat = size(human_tofit, 1);
df_full = length(paramBest);
df_chi2 = df_sat - df_full;
p_val = chi2cdf(G2, df_chi2, 'upper');    % if p>0.01, good fit


end

% ------------------------------------------------------------------------
% -> fill in struct
cncb_fit_struct = struct;
cncb_fit_struct.sens_noise = idl_sens_noise;     % sensory noise
cncb_fit_struct.sens_crit  = idl_sens_crit;      % sensory criterion
if (~only_efficiency)
    cncb_fit_struct.conf_noise = full_conf_noise;    % sdtev of noise for Type 2 decision (0 = ideal)
    cncb_fit_struct.conf_boost = full_conf_boost;    % fraction super-ideal (1 - fraction ideal)

    if (default_params_set.conf_crit > 0)
        cncb_fit_struct.conf_crit  = full_conf_crit;     % criterion for Type 2 decision (criterion)
    end
    
    full_struct = struct;
end


if (~skip_efficiency)
    cncb_fit_struct.efficiency = efficiency;                   % efficiency

    eff_struct = struct;
    eff_struct.equiv_conf_noise_ideal = equiv_noise_idl;  % equivalent noise 2 for ideal perf
    eff_struct.equiv_conf_noise_human = equiv_noise_hum;  % equivalent noise 2 for human perf

end

if (conf_continuous)
    if (~only_efficiency)
        eff_struct.llo_gamma_ideal       = model_idl_params.llo_gamma;
        eff_struct.llo_p0_ideal          = model_idl_params.llo_p0;
        eff_struct.llo_gamma_super_ideal = params_super_ideal.llo_gamma;
        eff_struct.llo_p0_super_ideal    = params_super_ideal.llo_p0;
        eff_struct.llo_gamma_super_human = params_super_human.llo_gamma;
        eff_struct.llo_p0_super_human    = params_super_human.llo_p0;

        full_struct.llo_gamma_full        = model_full_params.llo_gamma;
        full_struct.llo_p0_full           = model_full_params.llo_p0;
    end
else
    if (~skip_efficiency)
        eff_struct.conf_bnd_ideal       = bnd_lst_ideal;
        eff_struct.conf_bnd_super_ideal = super_ideal_rtng_bnd_lst;
        eff_struct.conf_bnd_super_human = super_human_rtng_bnd_lst;        
    end
    if (~only_efficiency)
        full_struct.conf_bnd_full        = sort(model_full_params.conf_bnds);
    end
end

if (~skip_efficiency)
    eff_struct.conf_rating_ideal_SRC       = ideal_rating_SRC;
    eff_struct.conf_rating_super_ideal_SRC = super_ideal_rating_SRC;
    eff_struct.conf_rating_super_human_SRC = super_human_rating_SRC;

    eff_struct.loglike_model = eff_model_loglike;
    eff_struct.loglike_saturated = eff_sat_loglike;
    eff_struct.df_model = eff_model_df;
    eff_struct.df_saturated = eff_sat_df;
    eff_struct.chi2_G2 = eff_G2;
    eff_struct.chi2_df = eff_df_chi2;
    eff_struct.chi2_p = eff_p_val;     % if p>0.01, good fit
    
    cncb_fit_struct.eff_struct = eff_struct;
end

if (~only_efficiency)
    full_struct.conf_rating_full_SRC    = full_rating_SRC;

    full_struct.loglike_model = model_full_loglike;  % log-likelihood of best fit
    full_struct.loglike_saturated = loglike_sat;
    full_struct.df_model = df_full;
    full_struct.df_saturated = df_sat;
    full_struct.chi2_G2 = G2;
    full_struct.chi2_df = df_chi2;
    full_struct.chi2_p = p_val;     % if p>0.01, good fit

    cncb_fit_struct.full_struct = full_struct;

end

if (verbose_flag >= 2)
    toc(tstart);
end

% -> EXIT <-


% ------------------------------------------------------------------------
% -> local functions


% -> This is the Type 1 function for confidence model fit
function yy_vals = cfc_type1(xx_vals, variable_prms, params_set, fixed_vals)

    noise_set_lst = params_set.sens_noise;
    crit_set_lst = params_set.sens_crit;
    noise_fix_lst = fixed_vals.sens_noise;
    crit_fix_lst = fixed_vals.sens_crit;

    var_ii = 1;

    if (noise_set_lst)
        noise = variable_prms(var_ii);
        var_ii = var_ii + 1;
    else
        noise = noise_fix_lst;
    end

    if (crit_set_lst)
        crit = variable_prms(var_ii);
    else
        crit = crit_fix_lst;
    end

    yy_vals = basic_normcdf(xx_vals, crit, noise);
end


% -> extract the last column of a matrix
function my_vec = lastcol(my_mat)
    my_vec = my_mat(:,end);
end


% -> extract counts of confidence judgments for each level
function my_vec = count_conf(my_mat, my_conf_nb, my_conf_ic)
    my_vec = NaN(my_conf_nb, 1);
    for my_cc = 1:my_conf_nb
        my_cond_inds = (my_conf_ic == my_cc);
        my_vec(my_cc) = sum(my_mat(my_cond_inds, 4));
    end
end


% -> wrapper around cfc_core function to flexibly add parameters
function pred_rating_mat = cncb_core_wrap(stims, variable_prms, params_set, fixed_vals)

    model_params = pack_params_in_struct(variable_prms, params_set, fixed_vals);
    pred_rating_mat = cncb_core(stims, model_params);
end


% -> wrapper around cfc_core function to flexibly add parameters
    function pred_rating_mat = cncb_core_count_wrap(stim_ext, variable_prms, params_set, fixed_vals)

    model_params = pack_params_in_struct(variable_prms, params_set, fixed_vals);
    pred_rating_mat = cncb_core_count(stim_ext, model_params); 
end


% -> get confidence counts rather than confidence probability 
function cncb_rating_mat = cncb_core_count(stim_ext, model_params)
    cncb_rating_mat = cncb_core(stim_ext(:,1), model_params);

    % -> scale confidence probabilities by counts for each stimulus
    for ii = 1:size(cncb_rating_mat, 1)
        my_stim_val = cncb_rating_mat(ii, 1);
        my_count_line = (stim_ext(:, 1) == my_stim_val);
        my_count_val = stim_ext(my_count_line, 2);
        cncb_rating_mat(ii, 4) = cncb_rating_mat(ii, 4) * my_count_val;
    end
end


% -> format parameters in a 'struct' variable
function params_struct = pack_params_in_struct(variable_prms, params_set, fixed_vals)

    params_struct = fixed_vals;
    params_cell = struct2cell(params_set);
    fld_nms = fieldnames(params_set);

    variable_nb = length(variable_prms);

    for kk = 1:variable_nb
        fld_ind = cellfun(@(xx) find(xx == kk), params_cell, 'UniformOutput', false);
        uu = find(~cellfun(@isempty, fld_ind));

        % -> allow for multiple use of a variable across model parameters
        for pp = 1:length(uu)
            vv = fld_ind{uu(pp)};

            % -> allow for multiple use of a variable within model parameters
            for qq = 1:length(vv)
                params_struct.(fld_nms{uu(pp)})(vv(qq)) = variable_prms(kk);
            end
        end
    end
end


% -> transform parameter structure into a vector of free parameters
function params_vals = extract_params_from_struct(params_set, params_struct)

    params_cell = struct2cell(params_set);
    fld_nms = fieldnames(params_set);

    params_mat = cell2mat(params_cell');
    variable_nb = max(params_mat);

    params_vals = NaN(1, variable_nb);

    for kk = 1:variable_nb
        fld_ind = cellfun(@(xx) find(xx == kk), params_cell, 'UniformOutput', false);
        uu = find(~cellfun(@isempty, fld_ind));
        vv = fld_ind{uu(1)};
        params_vals(kk) = params_struct.(fld_nms{uu(1)})(vv(1));
    end
end


% -> fit by maximizing the log-likelihood
function [params_best, loglike] = fitnll(fit_fcn, nn1_lst, ...
        params_0, params_LB, params_UB, fit_options)

    fun = @(xx) loglikefcn(xx, fit_fcn, nn1_lst);
    [params_best, nll_best] = fminsearchbnd(fun, ...
        params_0, params_LB, params_UB, fit_options);
    loglike = - nll_best;
end


% -> fit by maximizing the log-likelihood for binomials
function [params_best, loglike] = fitnllbin(fit_fcn, nn1_lst, ...
        params_0, params_LB, params_UB, fit_options)

    fun = @(xx) loglikebinfcn(xx, fit_fcn, nn1_lst);
    [params_best, nll_best] = fminsearchbnd(fun, ...
        params_0, params_LB, params_UB, fit_options);
    loglike = - nll_best;
end


% -> negative summed log-likelihood
function nglglk = loglikefcn(pp, ff, nn)

    ypred = ff(pp);
    ypred(ypred < 1e-6) = 1e-6;
    ypred(ypred > (1 - 1e-6)) = 1 - 1e-6;

    ll = log(ypred);

    ll_vect = nn .* ll;  % vector of log-likelihoods
    nglglk = - sum(ll_vect);      % minimize (-log) likelihood
end


% -> negative summed log-likelihood for binomials
%    this is the log-likelihood of the binomial functions, ignoring the
%    binomial coefficients that usually cancel out in models comparison
function nglglk = loglikebinfcn(pp, ff, nn)

    ypred = ff(pp);
    ypred(ypred < 1e-6) = 1e-6;
    ypred(ypred > (1 - 1e-6)) = 1 - 1e-6;

    ll1 = log(ypred);
    ll0 = log(1.0 - ypred);

    ll_vect = nn(:,1) .* ll1 + nn(:,2) .* ll0;  % vector of log-likelihoods
    nglglk = - sum(ll_vect);      % minimize (-log) likelihood
end


% ------------------------------------------------------------------------

end

% -> THE END <-
