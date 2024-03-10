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
    % llo_gamma = 0.8;        % log-odds parameters 'gamma' (default)
    llo_gamma = 0.5;        % log-odds parameters 'gamma' <<
    % llo_gamma = 2.0;        % log-odds parameters 'gamma'

    % llo_p0 = 0.5;           % log-odds parameters 'p0' (default)
    llo_p0 = 0.4;           % log-odds parameters 'p0' <<
    % llo_p0 = 0.8;           % log-odds parameters 'p0'

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

% -> extract stuff from fit
equiv_noise_hum          = cncb_fit_struct.equiv_conf_noise_human;
equiv_noise_idl          = cncb_fit_struct.equiv_conf_noise_ideal;
efficiency               = cncb_fit_struct.efficiency;

ideal_rating_mat         = cncb_fit_struct.conf_rating_ideal_SRC;
super_human_rating_mat   = cncb_fit_struct.conf_rating_super_human_SRC;
super_ideal_rating_mat   = cncb_fit_struct.conf_rating_super_ideal_SRC;

% -> print some stuff
fprintf('Human: Equivalent noise: %7.3f\n', equiv_noise_hum);
fprintf('Ideal: Equivalent noise: %7.3f\n', equiv_noise_idl);
fprintf('Confidence Efficiency:   %7.3f\n\n', efficiency);

best_fit_loglike         = cncb_fit_struct.loglike;

llo_gamma_full           = cncb_fit_struct.llo_gamma_full;
llo_p0_full              = cncb_fit_struct.llo_p0_full;


% -> create structure of parameters for best fit
model_best_params = struct;
model_best_params.sens_noise = cncb_fit_struct.sens_noise;
model_best_params.sens_crit  = cncb_fit_struct.sens_crit;
model_best_params.conf_noise = cncb_fit_struct.conf_noise;
model_best_params.conf_boost = cncb_fit_struct.conf_boost;
model_best_params.conf_crit  = cncb_fit_struct.conf_crit;
model_best_params.conf_continuous  = conf_continuous;
model_best_params.conf_cont_range  = conf_range;
if (exist('conf_half_scale', 'var'))
    model_best_params.conf_half_scale  = conf_half_scale;
end
model_best_params.llo_gamma = cncb_fit_struct.llo_gamma_full;
model_best_params.llo_p0 = cncb_fit_struct.llo_p0_full;

% model_best_params.conf_cont_nb_levels = 20;
model_best_params.conf_cont_nb_levels = 100;

model_SRC = cncb_core(sens_strengths, model_best_params);



fprintf('Goodness of fit:        %7.3f\n\n', best_fit_loglike);



% ------------------------------------------------------------------------
% -> PLOTS

% ********************************
%   Type 1 psychometric function
% ********************************

cncb_plot(cncb_data_grouped, 'type1_psychometric', true);


% ********************************
% -> plot choices from simulated observer against best fit of model
model2_SRC = cncb_fit_struct.conf_rating_full_SRC;
best_SRC = cncb_group(model2_SRC, ...
    'confidence_is_continuous', true, 'confidence_cont_range', conf_range, ...
    'confidence_cont_nb_levels', conf_nb_levels);
cncb_plot(cncb_data_grouped, 'human_model', best_SRC);


% ********************************
%   Type 2 ratings
% ********************************

% -> plot Type 2 ratings together with best fit
cncb_plot(cncb_data_grouped, 'type2_ratings', model_SRC, ...
    'confidence_cont_range', conf_range, ...
    'confidence_is_continuous', conf_continuous);
axis([0.0, 0.2, 0.0, 0.2]);


% ********************************
%   Type 2 ROC 
% ********************************

% -> plot Type 2 ROC for the original data with best fit
cncb_plot(cncb_data_grouped, 'type2_roc', model_SRC);


% -> THE END
