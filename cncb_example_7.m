% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.2
%
% cncb_example_7
%   Confidence example for simulations and fit.
%   Stimuli drawn from a range of strengths.
%   Use only 2 confidence levels.
%
% 19-AUG-2021 - pascal mamassian
% 14-FEB-2022 - pm: added continuous confidence resp
% 15-OCT-2022 - pm: cleaned up
% 10-MAR-2024 - pm: cleaned up
% 17-JUL-2024 - pm: added Deviance statistics
% 05-APR-2025 - pm: added 'sensory_bins' function to 'cncb_group'


% ----------------------
% -> prepare to simulate an experiment
% ----------------------

% clear all;
close all;

% ----------------------
% -> sensory parameters
% ----------------------

% -> stimulus strengths
sens_method = 0;     % continuous stimulus sampling
sens_intens_min = -2.0;
sens_intens_max =  2.0;

% -> for the fit, group data in bins
sens_bins = 4;

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
conf_continuous = false;


% -> confidence boundaries (in [0..1])
conf_bnds = normcdf(1.05);  % optimal for 2

rating_bnds_nb = length(conf_bnds);
rating_labels = (1:(rating_bnds_nb+1));


% ----------------------
% -> simulation parameters
% ----------------------

simul_orig_params = struct;

% -> draw stimuli from a range of difficulty levels
simul_orig_params.sens_method = sens_method;
simul_orig_params.sens_intens_min = sens_intens_min;
simul_orig_params.sens_intens_max = sens_intens_max;

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


% -> group the simulations by binning the sensory stimuli for the fit
cncb_data_grouped = cncb_group(raw_data, ...
    'sensory_bins', sens_bins, ...
    'confidence_is_continuous', conf_continuous, ...
    'confidence_disc_labels', rating_labels);


% ----------------------
% -> fit model to data
% ----------------------


data_tofit = cncb_data_grouped;

tic;
cncb_fit_struct = cncb_fit(data_tofit);
toc;


% ------------------------------------------------------------------------
% -> EXTRACT some stuff from the fit

% ----------------------
% -> efficiency fit
% ----------------------
equiv_noise_hum          = cncb_fit_struct.eff_struct.equiv_conf_noise_human;
equiv_noise_idl          = cncb_fit_struct.eff_struct.equiv_conf_noise_ideal;
efficiency               = cncb_fit_struct.efficiency;

best_fit_loglike         = cncb_fit_struct.eff_struct.loglike_model;
best_fit_G2              = cncb_fit_struct.eff_struct.chi2_G2;
best_fit_df              = cncb_fit_struct.eff_struct.chi2_df;
best_fit_p               = cncb_fit_struct.eff_struct.chi2_p;


% ----------------------
% -> full model fit
% ----------------------
% -> location of confidence boundaries
bnd_lst_full             = cncb_fit_struct.full_struct.conf_bnd_full;

% -> create structure of parameters for best fit
model_SRC = cncb_fit_struct.full_struct.conf_rating_full_SRC;


% ------------------------------------------------------------------------
% -> PRINT

% -> print some stuff from efficiency fit
fprintf('\nCNCB Efficiency:\n');
fprintf('Equivalent noise (Human): %7.3f\n', equiv_noise_hum);
fprintf('Equivalent noise (Ideal): %7.3f\n', equiv_noise_idl);
fprintf('Confidence Efficiency:    %7.3f\n', efficiency);

% -> print the goodness of fit of the efficiency computation
fprintf(['Goodness of fit:  log-likelihood = %7.3f, ', ...
    'chi2(%d) = %7.3f, p = %5.3f\n'], best_fit_loglike, ...
    best_fit_df, best_fit_G2, best_fit_p);

% -> print the confidence boundaries from the full model
fprintf('\nCNCB Full model:\n');
fprintf('Simulated confidence boundaries:  [');
fprintf('%7.3f  ', model_orig_params.conf_bnds);
fprintf(']\n');
fprintf('Estimated confidence boundaries:  [');
fprintf('%7.3f  ', bnd_lst_full);
fprintf(']\n');



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
cncb_plot(cncb_data_grouped, 'human_model', model_SRC);


% -> plot Type 2 ratings together with best fit
cncb_plot(cncb_data_grouped, 'type2_ratings', model_SRC, ...
    'confidence_is_continuous', conf_continuous);

% -> plot Type 2 ROC for the original data with best fit
cncb_plot(cncb_data_grouped, 'type2_roc', model_SRC);


% -> THE END
