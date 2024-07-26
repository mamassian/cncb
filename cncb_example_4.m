% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.1
%
% cncb_example_4
%   Confidence example for simulations and fit.
%   Scenario with continuous confidence judgments
%
% 19-AUG-2021 - pascal mamassian
% 14-FEB-2022 - pm: added continuous confidence resp
% 15-OCT-2022 - pm: cleaned up
% 10-MAR-2024 - pm: cleaned up
% 17-JUL-2024 - pm: added Deviance statistics


% ----------------------
% -> prepare to simulate an experiment
% ----------------------

% clear all;
close all;

% ----------------------
% -> sensory parameters
% ----------------------

% -> stimulus strengths
sens_strengths = [-1.0, 1.0];  % <<
nb_strengths = length(sens_strengths);

% -> noise of measurement (transduction): stdev
sens_noise = 1.0;  % <<

% -> sensory criterion
sens_crit = 0.0;    % <<


% ----------------------
% -> confidence parameters
% ----------------------

% -> added noise to combined (original + new sample)
conf_noise = 1.0;  

% -> boost towards super-ideal: 0 = ideal;  1 = super-ideal
conf_boost = 0.2;    % <<

% -> confidence criterion
default_conf_crit = 0.0;
conf_crit = 0.0;

% -> use discrete or continuous confidence judgments
conf_continuous = true;


% -> log-odds parameters
    % llo_gamma = 1.0;        % log-odds parameters 'gamma' (default)
    llo_gamma = 0.5;        % log-odds parameters 'gamma' <<

    % llo_p0 = 0.5;           % log-odds parameters 'p0' (default)
    llo_p0 = 0.4;           % log-odds parameters 'p0' <<

% -> full range
    conf_range = [0, 100];    % min and max confidence values
    conf_half_scale = false;
    
% -> half range
    % conf_range = [50, 100];    % min and max confidence values
    % conf_half_scale = true;
    
    conf_nb_levels = 12;      % resolution (number of points in interval)
    % conf_nb_levels = 100;      % resolution (number of points in interval)



% ----------------------
% -> simulation parameters
% ----------------------

simul_orig_params = struct;

% -> list of stimuli with different difficulty levels
simul_orig_params.sens_intens = sens_strengths;

% -> number of confidence judgments (associates to perceptual decisions)
simul_orig_params.nb_trials = 10000;    % <<



% ----------------------
% -> simulate the experiment
% ----------------------

% -> create structure of parameters for simulated experiment
model_orig_params = struct;
model_orig_params.sens_noise = sens_noise;
model_orig_params.sens_crit  = sens_crit;
model_orig_params.conf_noise = conf_noise;
model_orig_params.conf_boost = conf_boost;
model_orig_params.conf_crit  = conf_crit;


model_orig_params.conf_continuous  = conf_continuous;
model_orig_params.conf_range  = conf_range;

if (exist('conf_half_scale', 'var'))
    model_orig_params.conf_half_scale  = conf_half_scale;
end

model_orig_params.llo_gamma = llo_gamma;
model_orig_params.llo_p0 = llo_p0;



% -> simulate the experiment and store data in 'raw_data' matrix
simul_orig_data = cncb_simul(simul_orig_params, model_orig_params);
raw_data = simul_orig_data.raw_data;


cncb_data_grouped = cncb_group(raw_data, ...
    'confidence_is_continuous', true, ...
    'confidence_cont_range', conf_range, ...
    'confidence_half_scale', conf_half_scale, ...
    'confidence_cont_nb_levels', conf_nb_levels);


% ----------------------
% -> fit model to data
% ----------------------


tic;
    data_tofit = cncb_data_grouped;

    if (~exist('conf_half_scale', 'var'))
        cncb_fit_struct = cncb_fit(data_tofit, ...
            'confidence_is_continuous', conf_continuous, ...
            'confidence_cont_range', conf_range, ...
            'confidence_cont_nb_levels', conf_nb_levels);
    else
        cncb_fit_struct = cncb_fit(data_tofit, ...
            'confidence_is_continuous', conf_continuous, ...
            'confidence_cont_range', conf_range, ...
            'confidence_half_scale', conf_half_scale, ...
            'confidence_cont_nb_levels', conf_nb_levels);
    end
    
toc;


% ------------------------------------------------------------------------
% -> EXTRACT some stuff from the fit

% ----------------------
% -> efficiency fit
% ----------------------
efficiency               = cncb_fit_struct.efficiency;


% ----------------------
% -> full model fit
% ----------------------
best_fit_conf_noise     = cncb_fit_struct.conf_noise;
best_fit_conf_boost     = cncb_fit_struct.conf_boost;

best_fit_loglike        = cncb_fit_struct.full_struct.loglike_model;
best_fit_G2             = cncb_fit_struct.full_struct.chi2_G2;
best_fit_df             = cncb_fit_struct.full_struct.chi2_df;
best_fit_p              = cncb_fit_struct.full_struct.chi2_p;

best_fit_llo_gamma      = cncb_fit_struct.full_struct.llo_gamma_full;
best_fit_llo_p0         = cncb_fit_struct.full_struct.llo_p0_full;


% -> create structure of parameters for best fit
best_fit_params = struct;
best_fit_params.sens_noise = cncb_fit_struct.sens_noise;
best_fit_params.sens_crit  = cncb_fit_struct.sens_crit;
best_fit_params.conf_noise = cncb_fit_struct.conf_noise;
best_fit_params.conf_boost = cncb_fit_struct.conf_boost;
best_fit_params.conf_continuous  = conf_continuous;
best_fit_params.conf_cont_range  = conf_range;
if (exist('conf_half_scale', 'var'))
    best_fit_params.conf_half_scale  = conf_half_scale;
end
best_fit_params.llo_gamma = cncb_fit_struct.full_struct.llo_gamma_full;
best_fit_params.llo_p0    = cncb_fit_struct.full_struct.llo_p0_full;

% model_best_params.conf_cont_nb_levels = 20;
best_fit_params.conf_cont_nb_levels = 100;

model_SRC = cncb_core(sens_strengths, best_fit_params);


% ------------------------------------------------------------------------
% -> PRINT

% -> print stuff from efficiency
fprintf('\nCNCB Efficiency:\n');
fprintf('Confidence Efficiency:   %7.3f\n', efficiency);

% -> print stuff from full model
fprintf('\nCNCB Full model:\n');
fprintf('Estimated confidence noise: %7.3f\n', best_fit_conf_noise);
fprintf('Estimated confidence boost: %7.3f\n', best_fit_conf_boost);

fprintf('Estimated log-odds gamma:   %7.3f\n', best_fit_llo_gamma);
fprintf('Estimated log-odds p0:      %7.3f\n', best_fit_llo_p0);

fprintf(['Goodness of fit:  log-likelihood = %7.3f, ', ...
    'chi2(%d) = %7.3f, p = %5.3f\n\n'], best_fit_loglike, ...
    best_fit_df, best_fit_G2, best_fit_p);



% ------------------------------------------------------------------------
% -> PLOTS

% ********************************
%   Type 1 psychometric function
% ********************************
cncb_plot(cncb_data_grouped, 'type1_psychometric', true);


% ********************************
%   Type 2 ratings
% ********************************

% -> plot choices from simulated observer against best fit of model
model2_SRC = cncb_fit_struct.full_struct.conf_rating_full_SRC;
best_SRC = cncb_group(model2_SRC, ...
    'confidence_is_continuous', true, 'confidence_cont_range', conf_range, ...
    'confidence_cont_nb_levels', conf_nb_levels);
cncb_plot(cncb_data_grouped, 'human_model', best_SRC);
axis([0, 0.1, 0, 0.1]);


% -> plot Type 2 ratings together with best fit
%    warning: for continuous ratings, this plot rescales the confidence
%    probability so that the area under the curve for any stimulus is 1
cncb_plot(cncb_data_grouped, 'type2_ratings', model_SRC, ...
    'confidence_cont_range', conf_range, ...
    'confidence_is_continuous', conf_continuous);


% -> plot Type 2 ROC for the original data with best fit
cncb_plot(cncb_data_grouped, 'type2_roc', model_SRC);


% -> THE END
