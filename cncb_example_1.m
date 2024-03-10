% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.1
%
% cncb_example_1
%   Confidence example for simulations and fit.
%   Simple scenario with 2 stimulus strengths and 2 confidence levels
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
conf_continuous = false;


% -> confidence boundaries (in [0..1])
conf_bnds = normcdf(1.05);  % optimal for 2

rating_bnds_nb = length(conf_bnds);
rating_labels = (1:(rating_bnds_nb+1));

do_mratio_fit = 1;


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


data_tofit = simul_orig_data.data_SRC;

tic;
cncb_fit_struct = cncb_fit(data_tofit);
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

bnd_lst_full             = cncb_fit_struct.conf_bnd_full;


% -> create structure of parameters for best fit
model_SRC = cncb_fit_struct.conf_rating_full_SRC;




% -> print the confidence boundaries
fprintf('Simulated confidence boundaries:  [');
fprintf('%7.3f  ', model_orig_params.conf_bnds);
fprintf(']\n');
fprintf('Estimated confidence boundaries:  [');
fprintf('%7.3f  ', bnd_lst_full);
fprintf(']\n\n');

fprintf('Goodness of fit:        %7.3f\n\n', best_fit_loglike);


if (do_mratio_fit)
    [nR_S1, nR_S2] = CNCBtoML12(cncb_data_grouped);
    myfit = fit_meta_d_MLE(nR_S1, nR_S2);
    mratio = myfit.M_ratio;
    m_bnds1 = myfit.t2ca_rS1;
    m_bnds2 = myfit.t2ca_rS2;
    fprintf('m-ratio:  %7.3f\n', mratio);
    fprintf('meta-d'' confidence boundaries:  [');
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
% -> plot choices from simulated observer against best fit of model
best_SRC = cncb_fit_struct.conf_rating_full_SRC;
cncb_plot(cncb_data_grouped, 'human_model', best_SRC);



% ********************************
%   Type 2 ratings
% ********************************


% -> plot Type 2 ratings together with best fit
cncb_plot(cncb_data_grouped, 'type2_ratings', model_SRC, ...
    'confidence_is_continuous', conf_continuous);


% ********************************
%   Type 2 ROC 
% ********************************

% -> plot Type 2 ROC for the original data with best fit
cncb_plot(cncb_data_grouped, 'type2_roc', model_SRC);


% -> THE END
