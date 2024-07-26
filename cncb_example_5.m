% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.1
%
% cncb_example_5
%   Confidence example for simulations and fit.
%   Scenario with 2 stimulus strengths and 2 confidence levels
%   Examples of advances features
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
% sens_noise = 2.0;

% -> sensory criterion
sens_crit = 0.0;    % <<
% sens_crit = 0.25;


% ----------------------
% -> confidence parameters
% ----------------------

% -> added noise to combined (original + new sample)
% conf_noise = 2.0;
conf_noise = 1.0;  
% conf_noise = 0.5;    % <<
% conf_noise = 0.001;   % avoid 0.0

% -> boost towards super-ideal: 0 = ideal;  1 = super-ideal
% conf_boost = 0.8;
conf_boost = 0.2;    % <<
% conf_boost = 0.0;

% -> confidence criterion
default_conf_crit = 0.0;
conf_crit = 0.0;

% -> use discrete or continuous confidence judgments
conf_continuous = false;


% -> confidence boundaries (in [0..1])
conf_bnds = normcdf(1.05);  % optimal for 2

rating_bnds_nb = length(conf_bnds);
rating_labels = (1:(rating_bnds_nb+1));


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
model_orig_params.conf_bnds  = conf_bnds;



% -> simulate the experiment and store data in 'raw_data' matrix
simul_orig_data = cncb_simul(simul_orig_params, model_orig_params);
raw_data = simul_orig_data.raw_data;


cncb_data_grouped = cncb_group(raw_data, ...
    'confidence_is_continuous', conf_continuous, ...
    'confidence_disc_labels', rating_labels);


% ----------------------
% -> fit model to data
% ----------------------

% -> create structure of chosen parameters to fit for model
%    start with '1' and increment
params_set = struct;
params_set.sens_noise = 0;
params_set.sens_crit  = 0;
params_set.conf_noise = 1;
params_set.conf_boost = 2;

params_vals = struct;
params_vals.sens_noise = 1.00;
params_vals.sens_crit  = 0.00;

data_tofit = simul_orig_data.data_SRC;

tic;
cncb_fit_struct = cncb_fit(data_tofit, ...
    'model_parameters', params_set, ...
    'model_fixed_values', params_vals, ...
    'boost_init', [0.2, 0.5, 0.8]);
toc;


% ------------------------------------------------------------------------
% -> EXTRACT some stuff from the fit

% ----------------------
% -> efficiency fit
% ----------------------
equiv_noise_hum          = cncb_fit_struct.eff_struct.equiv_conf_noise_human;
equiv_noise_idl          = cncb_fit_struct.eff_struct.equiv_conf_noise_ideal;
efficiency               = cncb_fit_struct.efficiency;

ideal_rating_mat         = cncb_fit_struct.eff_struct.conf_rating_ideal_SRC;
super_human_rating_mat   = cncb_fit_struct.eff_struct.conf_rating_super_human_SRC;
super_ideal_rating_mat   = cncb_fit_struct.eff_struct.conf_rating_super_ideal_SRC;


% ----------------------
% -> full model fit
% ----------------------
best_fit_conf_noise     = cncb_fit_struct.conf_noise;
best_fit_conf_boost     = cncb_fit_struct.conf_boost;

best_fit_loglike        = cncb_fit_struct.full_struct.loglike_model;
best_fit_G2             = cncb_fit_struct.full_struct.chi2_G2;
best_fit_df             = cncb_fit_struct.full_struct.chi2_df;
best_fit_p              = cncb_fit_struct.full_struct.chi2_p;

bnd_lst_full            = cncb_fit_struct.full_struct.conf_bnd_full;

% -> create structure of parameters for best fit
model_SRC = cncb_fit_struct.full_struct.conf_rating_full_SRC;


% ------------------------------------------------------------------------
% -> PRINT

% -> print some stuff from efficiency fit
fprintf('\nCNCB Efficiency:\n');
fprintf('Equivalent noise (Human): %7.3f\n', equiv_noise_hum);
fprintf('Equivalent noise (Ideal): %7.3f\n', equiv_noise_idl);
fprintf('Confidence Efficiency:    %7.3f\n', efficiency);

% -> print stuff from full model
fprintf('\nCNCB Full model:\n');
fprintf('Estimated confidence noise: %7.3f\n', best_fit_conf_noise);
fprintf('Estimated confidence boost: %7.3f\n', best_fit_conf_boost);

% -> print the confidence boundaries
fprintf('Simulated confidence boundaries:  [');
fprintf('%7.3f  ', model_orig_params.conf_bnds);
fprintf(']\n');
fprintf('Estimated confidence boundaries:  [');
fprintf('%7.3f  ', bnd_lst_full);
fprintf(']\n');

% -> print the goodness of fit
fprintf(['Goodness of fit:  log-likelihood = %7.3f, ', ...
    'chi2(%d) = %7.3f, p = %5.3f\n'], best_fit_loglike, ...
    best_fit_df, best_fit_G2, best_fit_p);


% ------------------------------------------------------------------------
% -> PLOTS

% ********************************
%   Type 1 psychometric function
% ********************************
cncb_plot(cncb_data_grouped, 'type1_psychometric', true);


% ********************************
%   Type 2 choices against model
% ********************************

% -> plot choices from simulated observer against best fit of model
cncb_plot(cncb_data_grouped, 'human_model', model_SRC);


% -> plot choices from simulated observer against ideal goodness:
cncb_plot(cncb_data_grouped, 'human_model', ideal_rating_mat);
xlabel('Ideal Choice');


% ********************************
%   Type 2 ratings
% ********************************

% -> 1. plot Type 2 ratings just with the original data
cncb_plot(cncb_data_grouped, 'type2_ratings', true);
for sb = 1:nb_strengths
    subplot(1,2,sb);
    set(gca, 'XTick', [-0.4, -0.2, 0.0, 0.2, 0.4]);
    set(gca, 'XTickLabel', {'0.4', '0.2', '0', '0.2', '0.4'});
    set(gca, 'YTick', [1.0, 2.0]);
    set(gca, 'YTickLabel', {'low', 'high'});
end

% -> 2. plot Type 2 ratings together with best fit
cncb_plot(cncb_data_grouped, 'type2_ratings', model_SRC, ...
    'confidence_is_continuous', conf_continuous);


% ********************************
%   Type 2 ROC 
% ********************************

% -> 1. plot Type 2 ROC just with the original data
cncb_plot(cncb_data_grouped, 'type2_roc', true);
text(0.6, 0.1, 'no fit', 'FontSize', 20);

% -> 2. example of Type 2 ROC call with a model: e.g. the ideal observer
model_params = model_orig_params;
model_params.conf_noise = 0.001;    % singularity if 0
model_params.conf_boost = 0.0;
model_params.conf_bnds  = cncb_fit_struct.eff_struct.conf_bnd_ideal;
cncb_plot(cncb_data_grouped, 'type2_roc', model_params);
text(0.6, 0.1, 'ideal obs.', 'FontSize', 20);

% -> 3. plot Type 2 ROC for the original data with best fit
cncb_plot(cncb_data_grouped, 'type2_roc', model_SRC);
text(0.6, 0.1, 'best fit', 'FontSize', 20);


% -> THE END
