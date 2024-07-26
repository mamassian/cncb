% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.1
%
% cncb_simul
%   This function generates simulation of data of confidence rating
%   for a Type 1 discrimination task. 
%   Its output data are both in a neutral format (stim, resp, conf) and
%   in the format used by Maniscalco & Lau (2012).
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
% OUTPUT: cncb_data = struct
%   'raw_data'       : raw data
%   'data_SRC'       : data in format (stimulus, response, confidence bin, nb)
%                      where 'nb' is the count per (stim, resp, conf)
%   'data_SRC_prob'  : data in format (stimulus, response, confidence bin, prob)
%                      where 'prob' is the probability of (resp, conf | stim)
%   'nR_S1'          : data in format of Maniscalco & Lau for S1
%   'nR_S2'          : data in format of Maniscalco & Lau for S2
%
% both nR_S1 and nR_S2 are vectors containing the total number of responses
% in each response category, conditional on presentation of S1 and S2.
%
% e.g. if nR_S1 = [100 50 20 10 5 1], then when stimulus S1 was
% presented, the subject had the following response counts:
% responded S1, rating=3 : 100 times
% responded S1, rating=2 : 50 times
% responded S1, rating=1 : 20 times
% responded S2, rating=1 : 10 times
% responded S2, rating=2 : 5 times
% responded S2, rating=3 : 1 time
%
% The ordering of response / rating counts for S2 should be the same as it
% is for S1. e.g. if nR_S2 = [3 7 8 12 27 89], then when stimulus S2 was
% presented, the subject had the following response counts:
% responded S1, rating=3 : 3 times
% responded S1, rating=2 : 7 times
% responded S1, rating=1 : 8 times
% responded S2, rating=1 : 12 times
% responded S2, rating=2 : 27 times
% responded S2, rating=3 : 89 times
%
%
% 19-AUG-2021 - pascal mamassian
% 19-OCT-2023 - pm: do not discretize continuous simulations
% 10-MAR-2024 - pm: cleaned up


function cncb_data = cncb_simul(simul_params, model_params)


% *************************************************************************
% -> STIMULI <-

sens_intens = simul_params.sens_intens;   % range of stimuli
nb_sens_intens = length(sens_intens);     % nb of stimuli

nb_trials = simul_params.nb_trials;     % number of trials


% -> MODEL PARAMETERS <-

% -> Type 1 parameters
sens_noise = model_params.sens_noise;  % sensory noise
sens_crit  = model_params.sens_crit;   % sensory criterion

% -> Type 2 parameters (compulsory)
conf_noise = model_params.conf_noise;  % confidence noise
conf_boost = model_params.conf_boost;  % confidence boost

% -> Type 2 parameters (optional)
% -> confidence criterion
if (any(strcmp(fieldnames(model_params), 'conf_crit')))
    conf_crit_task  = model_params.conf_crit;   
else
    conf_crit_task = 0.0;
end



% -> use discrete or continuous confidence judgments
if (any(strcmp(fieldnames(model_params), 'conf_continuous')))
    conf_continuous = model_params.conf_continuous;
else
    % -> if nothing specified, assume continuous confidence on [0, 1]
    conf_continuous = 1;
end

% -> log-odds parameters for continuous confidence judgments
if (any(strcmp(fieldnames(model_params), 'llo_gamma')))
    llo_gamma = model_params.llo_gamma;
else
    llo_gamma = 1.0;
end

if (any(strcmp(fieldnames(model_params), 'llo_p0')))
    llo_p0 = model_params.llo_p0;
else
    llo_p0 = 0.5;
end

% -> confidence boundaries
if (any(strcmp(fieldnames(model_params), 'conf_bnds')))
    conf_bnd_lst = model_params.conf_bnds;
    conf_levels_nb = length(conf_bnd_lst) + 1;

    % -> check confidence boundaries are probabilities
    if ((max(conf_bnd_lst) > 1.0) || (min(conf_bnd_lst) < 0.0))
        fprintf('Error: confidence boundaries are not probabilities\n');
        return;
    end
    conf_bnd_edges = [0, conf_bnd_lst, 1];

    % -> if confidence bounds specified, enforce discrete
    conf_continuous = 0;
end

conf_half_scale = false;    % by default, assume full range

if (conf_continuous)
    % -> confidence range
    if (any(strcmp(fieldnames(model_params), 'conf_range')))
        conf_range = model_params.conf_range;
    else
        conf_range = [0.0, 1.0];
    end
    
    if (any(strcmp(fieldnames(model_params), 'conf_half_scale')))
        conf_half_scale = model_params.conf_half_scale;
    end
end



if (conf_continuous)
    conf_bnd_min = conf_range(1);
    conf_bnd_max = conf_range(end);

    if (conf_half_scale)
        conf_bnd_min2 = conf_bnd_min - (conf_bnd_max - conf_bnd_min);
    else
        conf_bnd_min2 = conf_bnd_min;
    end

    lo2conf = @(x_lo) max(conf_bnd_min, conf_bnd_min2 + ...
        (conf_bnd_max - conf_bnd_min2) * x_lo);

end


if (conf_continuous)
    data_SRC = NaN;
    data_SRC_norm = NaN;
    nR_S1 = NaN;
    nR_S2 = NaN;
    
else
    nb_conds = nb_sens_intens * 2 * conf_levels_nb;
    
    % -> data in format: stimulus_intensity, response (0, 1), confidence_index
    %    the 4th column is: prob(conf_ind | stim)
    data_SRC = zeros(nb_conds, 4);
    data_SRC_norm = zeros(nb_conds, 4);


    for ss = 1:nb_sens_intens
        for rr = 1:2
            for cc = 1:conf_levels_nb
                cnd_ind = (ss - 1)*2*conf_levels_nb + (rr - 1)*conf_levels_nb + cc;
                data_SRC(cnd_ind, 1) = sens_intens(ss);
                data_SRC(cnd_ind, 2) = rr - 1;
                data_SRC(cnd_ind, 3) = cc;
            end
        end
    end

    nR_S1_R1 = zeros(1, conf_levels_nb);
    nR_S1_R2 = zeros(1, conf_levels_nb);
    nR_S2_R1 = zeros(1, conf_levels_nb);
    nR_S2_R2 = zeros(1, conf_levels_nb);
end



% *************************************************************************

stims_nn = NaN(nb_trials, 1);   % stimulus intensity for each trial
percs_nn = NaN(nb_trials, 1);   % percept for each trial
conf_ind_nn = NaN(nb_trials, 1);   % confidence rating index
conf_mag_nn = NaN(nb_trials, 1);   % confidence rating magnitude

for tt = 1:nb_trials

    % -> choose randomly a stimulus in the set of possible stimuli
    ind_intens = randi(nb_sens_intens);   % index of stimulus for interval 'intrv'
    stims_nn(tt) = sens_intens(ind_intens);  % (noise-less) stimulus intensity

    % -> take independent noisy samples of the stimuli
    sens_smpl = stims_nn(tt) + randn * sens_noise;

    % -> sensory decision based on side of sensory criterion
    resp1 = (sens_smpl > sens_crit);
    percs_nn(tt) = resp1;

    % -> add fraction super-ideal
    boosted_sens_smpl = (1 - conf_boost)*sens_smpl + ...
        conf_boost*stims_nn(tt);

    % -> confidence sample
    conf_smpl = boosted_sens_smpl - sens_crit - conf_crit_task;

    % -> scale to sensory sensitivity
    scaled_smpl = conf_smpl / sens_noise;

    % -> confidence evidence: corrupt by Type 2 noise
    conf_evd = scaled_smpl + randn * conf_noise;

    % -> confidence magnitude
    conf_mag = abs(conf_evd);

    % -> pseudo-perceptual decision based on confidence sample
    resp2 = (conf_evd > 0);    

    % -> if pseudo-perceptual decision is inconsistent with
    %    perceptual decision, that perceptual decision might well be wrong
    if (resp2 ~= resp1)
        conf_mag = - conf_mag;
        
        % -> if half-scale, crop confidence distribution
        if (conf_half_scale)
            conf_mag = 0.0;
        end
    end

    % -> transform into confidence probability
    conf_prob = normcdf(conf_mag);
        
    if (conf_continuous)
        % -> log-odds transform
        conf_llo = cncb_log_odds(conf_prob, llo_gamma, llo_p0);
        
        conf_mag_nn(tt) = lo2conf(conf_llo);
    else
        conf_mag_nn(tt) = conf_prob;
        conf_ind = discretize(conf_prob, conf_bnd_edges);
        conf_ind_nn(tt) = conf_ind;

        % -> fill in the format stim-resp-conf
        cnd_ind = (ind_intens - 1)*2*conf_levels_nb + resp1*conf_levels_nb + conf_ind;
        data_SRC(cnd_ind, 4) = data_SRC(cnd_ind, 4) + 1;

        % -> fill in the format from Maniscalco&Lau
        if (nb_sens_intens == 2)
            if (ind_intens == 1)
                if (~resp1)
                    nR_S1_R1(conf_ind) = nR_S1_R1(conf_ind) + 1;
                else
                    nR_S1_R2(conf_ind) = nR_S1_R2(conf_ind) + 1;
                end
            else
                if (~resp1)
                    nR_S2_R1(conf_ind) = nR_S2_R1(conf_ind) + 1;
                else
                    nR_S2_R2(conf_ind) = nR_S2_R2(conf_ind) + 1;
                end
            end
        end
    
    end

end


% -> compute conf_prob in last columns of 'data_SRC_norm'
if (~conf_continuous)
    for ss = 1:nb_sens_intens
        cond_i_min = (ss - 1)*2*conf_levels_nb + 1;
        cond_i_max = ss*2*conf_levels_nb;
        data_SRC_norm(cond_i_min:cond_i_max, 4) = ...
            data_SRC(cond_i_min:cond_i_max, 4) ./ sum(data_SRC(cond_i_min:cond_i_max, 4));

        for rr = 1:2
            for cc = 1:conf_levels_nb
                cnd_ind = (ss - 1)*2*conf_levels_nb + (rr - 1)*conf_levels_nb + cc;
                data_SRC(cnd_ind, 1) = sens_intens(ss);
                data_SRC(cnd_ind, 2) = rr - 1;
                if (conf_continuous)
                    data_SRC(cnd_ind, 3) = conf_bnd_labels(cc);
                else
                    data_SRC(cnd_ind, 3) = cc;
                end
            end
        end
        
    end
    data_SRC_norm(:, 1:3) = data_SRC(:, 1:3);
end

% -> format data for output
if (conf_continuous)
    raw_data = [stims_nn, percs_nn, conf_mag_nn];
else
    raw_data = [stims_nn, percs_nn, conf_ind_nn];
    if (nb_sens_intens == 2)
        nR_S1 = [fliplr(nR_S1_R1), nR_S1_R2];
        nR_S2 = [fliplr(nR_S2_R1), nR_S2_R2];
    end
end


% -> write output
cncb_data = struct;

cncb_data.raw_data = raw_data;
cncb_data.data_SRC = data_SRC;
cncb_data.data_SRC_prob = data_SRC_norm;
if (nb_sens_intens == 2)
    cncb_data.nR_S1 = nR_S1;
    cncb_data.nR_S2 = nR_S2;
end

end

% -> THE END <-
