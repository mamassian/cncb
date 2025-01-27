% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.2
%
% cncb_example_6
%   Confidence example for simulations and fit.
%   Scenario with 2 stimulus strengths and 4 confidence levels
%     with a confidence bias.
%
% 02-JAN-2025 - pascal mamassian


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
conf_boost = 0.2;

% -> confidence criterion
default_conf_crit = 0.0;
conf_crit = 0.0;

% -> use discrete or continuous confidence judgments
conf_continuous = false;


% -> number of confidence boundaries
rating_bnds_nb = 3;
rating_labels = (1:(rating_bnds_nb+1));


% -> confidence bias (ratio of conf. choices b/w top & bottom levels)
conf_bias = 2.0;


do_mratio_fit = 1;


% ----------------------
% -> simulation parameters
% ----------------------

simul_orig_params = struct;

% -> list of stimuli with different difficulty levels
simul_orig_params.sens_intens = sens_strengths;

% -> number of confidence judgments (associates to perceptual decisions)
simul_orig_params.nb_trials = 10000; 



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


% -> location of the confidence boundaries corresponding to confidence bias
conf_bnds = cncb_get_bnds(simul_orig_params, model_orig_params, ...
    rating_bnds_nb, conf_bias);

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



data_tofit = simul_orig_data.data_SRC;

tic;
cncb_fit_struct = cncb_fit(data_tofit);
toc;


% ------------------------------------------------------------------------
% -> EXTRACT some stuff from the fit

% ----------------------
% -> efficiency fit
% ----------------------
efficiency              = cncb_fit_struct.efficiency;


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

% -> print some stuff
fprintf('\nCNCB Efficiency:\n');
fprintf('Confidence Efficiency:   %7.3f\n', efficiency);

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
% -> meta-d' analysis (for comparison)
if (do_mratio_fit)
    [nR_S1, nR_S2] = CNCBtoML12(cncb_data_grouped);
    myfit = fit_meta_d_MLE(nR_S1, nR_S2);
    mratio = myfit.M_ratio;
    m_bnds1 = myfit.t2ca_rS1;
    m_bnds2 = myfit.t2ca_rS2;
    fprintf('\nMeta-d'':\n');
    fprintf('m-ratio:  %7.3f\n', mratio);
    fprintf('Estimated confidence boundaries:  [');
    fprintf('%7.3f  ', normcdf(m_bnds2));
    fprintf(']\n\n');
end


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
