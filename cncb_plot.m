% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.1
%
% cncb_plot
%   This is a collection of plot functions for confidence rating scale.
%
% INPUT:
%   'cncb_data': grouped data per (stim, resp, conf_id, conf_prob):
%       1st col: stimulus intensity
%       2nd col: perceptual decision
%       3rd col: confidence level (index 1, 2, ...)
%       4th col: confidence probability (for that conf_id | stim, resp)
% note: the 4th column can also just be the count of trials for this condition
%
% OUTPUT:
%   'cncb_plot_struct' = struct
%     cncb_plot_struct.sensory_strength    list of stimuli
%     cncb_plot_struct.t1_response_prob    Type1 responses (unsorted)
%     cncb_plot_struct.t1_response_count   nb. of Type1 responses (unsorted)
%     cncb_plot_struct.human_type2_roc     human Type 2 ROC [fa2; hit2]
%     cncb_plot_struct.model_type2_roc     model Type 2 ROC [fa2; hit2]
%   
%
% POSSIBLE PARAMETERS:
%   'human_model': plot confidence probability for human against a model 
%           provided as argument
%
%   'type1_psychometric': plot psychometric functions for percepts
%           if there is an argument 'arg', plot model data in 'arg'
%           format (n1, n0) for each 'knd'
%
%   'type2_roc': plot Type II ROC
%
%   'type2_ratings': plot confidence ratings
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
%
%
% 19-AUG-2021 - pascal mamassian
% 20-FEB-2022 - pm: added continuous ratings
% 15-OCT-2022 - pm: added 'cncb_plot_struct' output
% 25-OCT-2023 - pm: added continuous confidence plot
% 10-MAR-2024 - pm: cleaned up
% 01-APR-2024 - pm: fixed 'xmax' in continuous case
% 01-JUN-2024 - pm: reverse the 'x' and 'y' axes for confidence plots


function cncb_plot_struct = cncb_plot(cncb_data, varargin)

% -> default optional arguments
dflt_human_model         = [];	   % human against model choice
dflt_t1_psychometric     = false;  % type-1 psychometric functions
dflt_t2_ratings          = false;  % type-2 ratings
dflt_t2_roc              = false;  % type-2 ROC: logical, model_params, or rating_mat
dflt_t2_zroc             = false;  % type-2 zROC: logical, model_params, or rating_mat
dflt_is_continuous       = false;  % continuous (T) or discrete (F) rating scale
dflt_conf_cont_half_scale  = false;
dflt_conf_cont_range     = [0, 1];
dflt_conf_cont_nb_levels = 100;

% -> parse all arguments
ip = inputParser;
addRequired(ip, 'cncb_data', @isnumeric);
addParameter(ip, 'human_model', dflt_human_model, @isnumeric);
addParameter(ip, 'type1_psychometric', dflt_t1_psychometric);
addParameter(ip, 'type2_roc', dflt_t2_roc);
addParameter(ip, 'type2_zroc', dflt_t2_zroc);
addParameter(ip, 'type2_ratings', dflt_t2_ratings);
addParameter(ip, 'confidence_is_continuous', dflt_is_continuous, @islogical);
addParameter(ip, 'confidence_half_scale', dflt_conf_cont_half_scale, @islogical);
addParameter(ip, 'confidence_cont_range', dflt_conf_cont_range, @isnumeric);
addParameter(ip, 'confidence_cont_nb_levels', dflt_conf_cont_nb_levels, @isnumeric);

parse(ip, cncb_data, varargin{:});
plot_human_model = ip.Results.human_model;
plot_t1_psychometric = ip.Results.type1_psychometric;
plot_t2_ratings = ip.Results.type2_ratings;
plot_t2_roc = ip.Results.type2_roc;
plot_t2_zroc = ip.Results.type2_zroc;
conf_continuous = ip.Results.confidence_is_continuous;
conf_half_scale = ip.Results.confidence_half_scale;
conf_cont_range = ip.Results.confidence_cont_range;
conf_cont_nb_levels = ip.Results.confidence_cont_nb_levels;

% -> prepare output variables
cncb_plot_struct = struct;


col_stim = 1;
col_resp = 2;
col_conf_levl = 3;
col_conf_prob = 4;

% -> parameters for 'type2_ratings':
%    set thresholds on number of confidence levels to display different stuff
many_confs1 = 6;
many_confs2 = 12;

% -> random stuff
mtlb_colors = [...
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880;
    0.3010, 0.7450, 0.9330;
    0.6350, 0.0780, 0.1840];

% -> ********************************************** <-

% -> if raw data have only 3 columns, add a 4th one
if (size(cncb_data, 2) < col_conf_prob)
    nb_lines = size(cncb_data, 1);
    cncb_data = [cncb_data, ones(nb_lines, 1)];
end

% -> clean up input data
if (conf_continuous)
    cncb_data = cncb_group(cncb_data, ...
        'confidence_is_continuous', conf_continuous, ...
        'confidence_cont_range', conf_cont_range, ...
        'confidence_half_scale', conf_half_scale, ...
        'confidence_cont_nb_levels', conf_cont_nb_levels);
else
    cncb_data = cncb_group(cncb_data, ...
        'confidence_is_continuous', conf_continuous);
end

[stim_lst, ~, stim_ic] = unique(cncb_data(:, col_stim));
stim_nb = length(stim_lst);

[resp_lst, ~, resp_ic] = unique(cncb_data(:, col_resp));
resp_nb = length(resp_lst);

[conf_lst, ~, conf_ic] = unique(cncb_data(:, col_conf_levl));
conf_nb = length(conf_lst);

% -> nb kinds of trials (stim, resp, conf)
nb_knds = size(cncb_data, 1);

% -> un-wrap data matrix
resps2 = cncb_data(:, col_resp);
choic2 = cncb_data(:, col_conf_prob);  % nb trials per (stim, resp)


% -> internal function variables for 'plot_t1_psychometric'
resp_prob_lst     = NaN(1, stim_nb);
resp_count_lst    = NaN(1, stim_nb);

for ww = 1:stim_nb
    inds = (stim_ic == ww);      % indices that have the same stim

    nb_resps = sum(choic2(inds));       % nb. of responses for this particular stim
    resp_vals = resps2(inds) .* choic2(inds);
    rsp_prob = sum(resp_vals) / nb_resps;   % unsorted responses
    resp_prob_lst(ww) = rsp_prob;
    resp_count_lst(ww) = nb_resps;
end


cncb_plot_struct.sensory_strength = stim_lst;
cncb_plot_struct.t1_response_prob = resp_prob_lst;
cncb_plot_struct.t1_response_count = resp_count_lst;


if (conf_continuous)
    conf_bnd_min = conf_cont_range(1);
    conf_bnd_max = conf_cont_range(end);

    % -> edges of confidence ratings
    conf_labels_step = (conf_bnd_max - conf_bnd_min) / conf_cont_nb_levels;
    conf_bnd_edges = conf_bnd_min:conf_labels_step:conf_bnd_max;
    conf_bnd_labels = conf_bnd_edges(2:end);
else
    [conf_bnd_labels, ~, ~] = unique(cncb_data(:, col_conf_levl));
end
conf_levels_nb = length(conf_bnd_labels);
rating_bnds_nb = conf_levels_nb - 1;


% -> normalize to get prob(conf | stim)
cncb_data_norm = cncb_data;
data_SRC_norm = NaN(size(cncb_data_norm,1),1);
for ss = 1:stim_nb
    sens_mean = stim_lst(ss);
    cond_inds = (cncb_data(:, col_stim) == sens_mean);

    data_SRC_norm(cond_inds) = ...
        cncb_data(cond_inds, col_conf_prob) ./ ...
        sum(cncb_data(cond_inds, col_conf_prob));
end
cncb_data_norm(:, col_conf_prob) = data_SRC_norm;


% -> ******************************************************************* <-
% -> plot psychometric function for perceptual decisions
goaheadandplot = 0;
if (islogical(plot_t1_psychometric))
    if (plot_t1_psychometric)
        goaheadandplot = 1;
    end
elseif (isstruct(plot_t1_psychometric))
    goaheadandplot = 2;
end

if (goaheadandplot)
    fig = figure;
    fig_pos = get(fig, 'Position');
    set(fig, 'Position', [fig_pos(1) fig_pos(2) 430 340]);

    hold on;
    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);

    xvals = stim_lst;
    xrange = xvals(end) - xvals(1);
    xmin = xvals(1) - xrange/10;
    xmax = xvals(end) + xrange/10;

    line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    line([xmin xmax], [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

    max_size = max(resp_count_lst) * 0.4;
    col = [0.5, 0.5, 0.5];
    sizes = 100 * resp_count_lst ./ max_size + 0.01;
    sizes(sizes==0) = 0.1;
    scatter(stim_lst, resp_prob_lst, sizes, col, 'filled', 'LineWidth', 2);

    if (goaheadandplot == 2)
        sens_noise = plot_t1_psychometric.sens_noise;  % sensory noise
        sens_crit  = plot_t1_psychometric.sens_crit;   % sensory criterion
        xvals = linspace(xmin, xmax, 100);
        yvals = basic_normcdf(xvals, sens_crit, sens_noise);
        plot(xvals, yvals, '-', 'Color', col, 'LineWidth', 3);

        cncb_plot_struct.t1_psychometric = [xvals; yvals];
    else
        plot(stim_lst, resp_prob_lst, '-', 'Color', col, 'LineWidth', 3);
    end

    xlim([xmin, xmax]);
    ylim([0, 1]);
    xlabel('Stimulus Strength');
    ylabel('Proportion ''Right''');
end


% -> ******************************************************************* <-
% -> plot of fraction chosen for human against model
if (size(plot_human_model, 1) ~= 0)

    if (length(stim_ic) ~= size(plot_human_model, 1))
        fprintf('Error: model does not have the same size as data\n');
        return;
    end

    % -> fraction confidence for each trial kind
    human_chosen = NaN(1, nb_knds);
    model_chosen = NaN(1, nb_knds);
    total1 = NaN(1, nb_knds);

    [~, ~, mdl_resp_ic] = unique(plot_human_model(:, col_resp));
    [~, ~, mdl_conf_ic] = unique(plot_human_model(:, col_conf_levl));
    mdl_choic2 = plot_human_model(:, col_conf_prob);  % nb trials per (stim, resp)

    for ww = 1:stim_nb
        logical_filter_ss = (stim_ic == ww);

        % -> nb. of responses for this particular stim
        hum_nb_resps = sum(choic2(logical_filter_ss));
        mdl_nb_resps = sum(mdl_choic2(logical_filter_ss));

        for rr = 1:resp_nb
            hum_logical_filter_rr = (logical_filter_ss & (resp_ic == rr));
            mdl_logical_filter_rr = (logical_filter_ss & (mdl_resp_ic == rr));

            total_resp = sum(choic2(hum_logical_filter_rr));

            for cc = 1:conf_nb
                conf_ind = (ww - 1)*(resp_nb*conf_nb) + (rr - 1)*conf_nb + cc;
                total1(conf_ind) = total_resp;

                % -> 'human_chosen' should be the same as 'cncb_data_norm(:,4)',
                %    but recompute it just to be consistent with 'model_chosen'
                hum_logical_filter = (hum_logical_filter_rr & (conf_ic == cc));
                conf_val = sum(cncb_data(hum_logical_filter, col_conf_prob));
                conf_val = conf_val / hum_nb_resps;
                human_chosen(conf_ind) = conf_val;

                mdl_logical_filter = (mdl_logical_filter_rr & (mdl_conf_ic == cc));
                conf_val = sum(plot_human_model(mdl_logical_filter, col_conf_prob));
                conf_val = conf_val / mdl_nb_resps;
                model_chosen(conf_ind) = conf_val;
            end
        end
    end

    % -> take care of the rare event when a pair (stim, resp) is not happening
    total1(total1==0) = NaN;

    fig = figure;
    fig_pos = get(fig, 'Position');
    set(fig, 'Position', [fig_pos(1) fig_pos(2) 430 340]);
    hold on;
    set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);
    line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

    sizes = 100 * total1 ./ max(total1);
    col = [0.8, 0.0, 0.5];
    scatter(model_chosen, human_chosen, sizes, col, 'filled', 'LineWidth', 2);

    axis([0, 1, 0, 1]);
    axis('square');
    xlabel('Model Choice');
    ylabel('Human Choice');
end



% -> ******************************************************************* <-
% -> plot Type 2 ratings
if (islogical(plot_t2_ratings))
    if (plot_t2_ratings)
        plot_t2_ratings_cond = 1;
    else
        plot_t2_ratings_cond = 0;
    end
elseif (ismatrix(plot_t2_ratings))
    % -> a struct is also an array, but of size 1
    if (size(plot_t2_ratings, 1) == 1)
        plot_t2_ratings_cond = 3;
        model_params = plot_t2_ratings;
        rating_mat = cncb_core(stim_lst, model_params);
    else    
        plot_t2_ratings_cond = 2;
        rating_mat = plot_t2_ratings;
    end
    [~, ~, mod_stim_ic] = unique(rating_mat(:, col_stim));
    [~, ~, mod_resp_ic] = unique(rating_mat(:, col_resp));
end

if (plot_t2_ratings_cond)
    fig = figure('Name', 'Type 2 ratings');
    fig_pos = get(fig, 'Position');
    if (stim_nb == 2)
        set(fig, 'Position', [fig_pos(1), fig_pos(2), 770, 420]);
    else
        set(fig, 'Position', [fig_pos(1), fig_pos(2), 850, 660]);
    end

    leg_label = strings(1, stim_nb);
    xmax = 0.0;

    for ss = 1:stim_nb

        if (stim_nb == 2)
            subplot(1, 2, ss);
        else
            max_panels = ceil(stim_nb/2);
            if (ss <= max_panels)
                ind_panel = ss;
            else
                ind_panel = stim_nb - ss + max_panels + 1;
            end
            subplot(2, max_panels, ind_panel);
        end
        
        set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);
        hold on;
        
        stim_val = stim_lst(ss);
        leg_label(ss) = sprintf('stim: %5.2f', stim_val);
        title(leg_label(ss));
        
        for rr = 1:resp_nb
            
            if (rr == 1)
                my_col = [1, 0.5, 0];
            elseif (rr == 2)
                my_col = [0, 0.5, 1];
            else
            my_col = mtlb_colors(rr, :);
            end
            
            conf_labels = cncb_data_norm((stim_ic==ss) & (resp_ic==rr), ...
                col_conf_levl);
            conf_counts = cncb_data_norm((stim_ic==ss) & (resp_ic==rr), ...
                col_conf_prob);
            hist_nb_bins = length(conf_labels);
            inter_label = conf_labels(2) - conf_labels(1);
            
            % -> adjust label to represent end of bin
            too_many_confs = (hist_nb_bins > many_confs1);
            
            if (rr == 1)
                conf_counts = - conf_counts;
            end
            
            if ((too_many_confs) || (conf_continuous))
                conf_width = inter_label;
            else
                hh = bar(conf_labels, conf_counts, 'FaceColor', [1, 1, 1], ...
                    'EdgeColor', my_col, 'LineWidth', 3);
                hh.EdgeAlpha = 0.5;

                conf_width = hh.BarWidth;
            end
            
            for conf_ii = 1:hist_nb_bins
                conf_xx = conf_counts(conf_ii);

                if (conf_continuous)
                    if (conf_ii ~= 1)
                        conf_y1 = conf_labels(conf_ii-1);
                    else
                        if (~ismember('confidence_cont_range', ip.UsingDefaults))
                            conf_y1 = conf_cont_range(1);
                        else
                            conf_y1 = 2*conf_labels(conf_ii) - conf_labels(conf_ii+1);
                        end
                    end
                    conf_y2 = conf_labels(conf_ii);

                    % -> normalize to get probabiity
                    conf_xx = conf_xx / (conf_y2 - conf_y1);
                else
                    conf_y1 = conf_labels(conf_ii) - conf_width/2;
                    conf_y2 = conf_labels(conf_ii) + conf_width/2;
                end

                xmax = max(xmax, abs(conf_xx));

                conf_xvals = [conf_xx, 0.0, 0.0, conf_xx];
                conf_yvals = [conf_y1, conf_y1, conf_y2, conf_y2];
                
                alpha_val = 0.2 + 0.6 * (conf_ii / hist_nb_bins);
                fill(conf_yvals, conf_xvals, my_col, 'FaceAlpha', alpha_val, ...
                    'EdgeColor', 'none');
            end
            

            % -> draw a black base for the bars
            %    (there is a bug to be fixed for continuous distributions)
            conf_labels = unique(conf_labels);
            bar(conf_labels, zeros(size(conf_labels)), 1.0, 'FaceColor', [1, 1, 1], ...
                'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 3);

            % -> plot model
            if (plot_t2_ratings_cond >= 2)
                conf_model = rating_mat((mod_stim_ic==ss) & (mod_resp_ic==rr), ...
                    col_conf_prob);
                mod_conf_labels = rating_mat((mod_stim_ic==ss) & (mod_resp_ic==rr), ...
                    col_conf_levl);

                xmax = max(xmax, max(abs(conf_model)));

                for conf_ii = 1:length(mod_conf_labels)
                    if (conf_ii ~= 1)
                        conf_y1 = mod_conf_labels(conf_ii-1);
                    else
                        if (~ismember('confidence_cont_range', ip.UsingDefaults))
                            conf_y1 = conf_cont_range(1);
                        else
                            conf_y1 = 2*mod_conf_labels(conf_ii) - mod_conf_labels(conf_ii+1);
                        end
                    end
                    conf_y2 = mod_conf_labels(conf_ii);
    
                    if (conf_continuous)
                        % -> normalize to get probabiity
                        conf_model(conf_ii) = conf_model(conf_ii) / (conf_y2 - conf_y1);
                    end
                end


                if (~conf_continuous)
                    % -> replace confidence boundaries by labels
                    mod_conf_labels = (1:length(mod_conf_labels));
                end


                if (rr == 1)
                    conf_model = - conf_model;
                end

                if (length(mod_conf_labels) <= many_confs2)
                    plot(mod_conf_labels, conf_model, 'o-', ...
                        'MarkerSize', 14, 'MarkerFaceColor', [1 1 1], ...
                        'Color', my_col, 'LineWidth', 4);
                else
                    plot(mod_conf_labels, conf_model, '-', ...
                        'Color', my_col, 'LineWidth', 4);
                end
            end
            
        end

        if (~ismember('confidence_cont_range', ip.UsingDefaults))
            xlim(conf_cont_range);
        else
            xlim([conf_labels(1)-0.5*inter_label, ...
                conf_labels(end)+0.5*inter_label]);
            if (hist_nb_bins <= 6)
                set(gca, 'XTick', conf_labels);
            end
        end
        
        if (stim_nb == 2)
            ylabel('Confidence Prob.', 'FontSize', 26);
            if (ss == 1)
                xlabel('Confidence Level', 'FontSize', 26);
            end
        elseif (ss == stim_nb)
            ylabel('Confidence Prob.', 'FontSize', 26);
            xlabel('Confidence Level', 'FontSize', 26);
        end

    end
    
    % -> adjust y-axes
    xmax2 = floor(log10(xmax));
    xmax = ceil(xmax / (10^xmax2)) * (10^xmax2);
    ytic = linspace(-xmax, xmax, 5);
    ylbl = num2cell(abs(ytic));
    
    for ss = 1:stim_nb
        if (stim_nb == 2)
            subplot(1, 2, ss);
        else
            if (ss <= max_panels)
                ind_panel = ss;
            else
                ind_panel = stim_nb - ss + max_panels + 1;
            end
            subplot(2, max_panels, ind_panel);
        end
        ylim([-xmax, xmax]);
        set(gca, 'YTick', ytic, 'YTickLabel', ylbl);
    end
    
end



% -> ******************************************************************* <-
% -> plot Type 2 ROC
if ((ismatrix(plot_t2_roc)) && ~(islogical(plot_t2_roc)))
    % -> a struct is also an array, but of size 1
    if (size(plot_t2_roc, 1) == 1)
        plot_t2_roc_cond = 3;
        model_params = plot_t2_roc;
        rating_mat = cncb_core(stim_lst, model_params);
    else    
        plot_t2_roc_cond = 2;
        rating_mat = plot_t2_roc;
    end
    plot_zscore = 0;

elseif ((ismatrix(plot_t2_zroc)) && ~(islogical(plot_t2_zroc)))
    if (size(plot_t2_zroc, 1) == 1)
        plot_t2_roc_cond = 3;
        model_params = plot_t2_zroc;
        rating_mat = cncb_core(stim_lst, model_params);
    else    
        plot_t2_roc_cond = 2;
        rating_mat = plot_t2_zroc;
    end
    plot_zscore = 1;

elseif ((islogical(plot_t2_roc)) || (islogical(plot_t2_zroc)))
    if ((plot_t2_roc) || (plot_t2_zroc))
        plot_t2_roc_cond = 1;
        if (plot_t2_zroc)
            plot_zscore = 1;
        else
            plot_zscore = 0;
        end
    else
        plot_t2_roc_cond = 0;
    end

else
    plot_t2_roc_cond = 0;
end


if (plot_t2_roc_cond)
    hum_roc_lst = src_to_roc(cncb_data);

    if (plot_zscore)
        fig = figure('Name', 'Type 2 zROC');
    else
        fig = figure('Name', 'Type 2 ROC');
    end
    fig_pos = get(fig, 'Position');
    
    if (stim_nb == 2)
        set(fig, 'Position', [fig_pos(1), fig_pos(2), 770, 420]);
    else
        set(fig, 'Position', [fig_pos(1), fig_pos(2), 850, 660]);
    end

    leg_label = strings(1, stim_nb);
    
    for mm = 1:stim_nb
        if (stim_nb == 2)
            subplot(1, 2, mm);
        else
            max_panels = ceil(stim_nb/2);
            if (mm <= max_panels)
                ind_panel = mm;
            else
                ind_panel = stim_nb - mm + max_panels + 1;
            end
            subplot(2, max_panels, ind_panel);
        end
        
        set(gca, 'FontName', 'Arial'); set(gca, 'FontSize', 24);
        hold on;

        if (plot_zscore)
            line([-3, 3], [-3, 3], 'LineStyle', '-', ...
                'Color', [0.5 0.5 0.5], 'LineWidth', 2);
        else
            line([0, 1], [0, 1], 'LineStyle', '-', ...
                'Color', [0.5 0.5 0.5], 'LineWidth', 2);
        end

        stim_val = stim_lst(mm);
        leg_label(mm) = sprintf('stim: %5.2f', stim_val);
        title(leg_label(mm));

        my_col1 = [0, 0, 0];
        my_col2 = [0.5, 0.5, 0.5];
        
        human_rng_inds = (mm-1)*rating_bnds_nb + (1:rating_bnds_nb);

        human_hit2_lst = hum_roc_lst(human_rng_inds, 3);
        human_fa2_lst  = hum_roc_lst(human_rng_inds, 4);

        human_hit2_lst2 = [1.0; human_hit2_lst; 0.0];
        human_fa2_lst2 = [1.0; human_fa2_lst; 0.0];

        if (plot_zscore)
            human_hit2_lst = norminv(human_hit2_lst);
            human_fa2_lst  = norminv(human_fa2_lst);
        end
        
        if (conf_levels_nb <= many_confs2)
            plot(human_fa2_lst, human_hit2_lst, 'o', ...
                'MarkerSize', 16, 'MarkerFaceColor', my_col2, ...
                'MarkerEdgeColor', 'none');
        else
            plot(human_fa2_lst, human_hit2_lst, 'o', ...
                'MarkerSize', 12, 'MarkerFaceColor', my_col2, ...
                'MarkerEdgeColor', 'none');
        end
        
        if (plot_t2_roc_cond == 1)
            if (plot_zscore)
                human_fa2_lst2 = norminv(human_fa2_lst2);
                human_hit2_lst2 = norminv(human_hit2_lst2);
            end
            plot(human_fa2_lst2, human_hit2_lst2, '-', ...
                'Color', my_col1, 'LineWidth', 3);

        elseif ((plot_t2_roc_cond == 2) || (plot_t2_roc_cond == 3))
            mdl_conf_levels_nb = size(rating_mat,1) / stim_nb / resp_nb;
            mdl_rating_bnds_nb = mdl_conf_levels_nb - 1;
            model_rng_inds = (mm-1)*mdl_rating_bnds_nb + (1:mdl_rating_bnds_nb);
            
            mdl_roc_lst = src_to_roc(rating_mat);
            model_hit2_lst = mdl_roc_lst(model_rng_inds, 3);
            model_fa2_lst  = mdl_roc_lst(model_rng_inds, 4);

            model_hit2_lst2 = [1.0; model_hit2_lst; 0.0];
            model_fa2_lst2 = [1.0; model_fa2_lst; 0.0];
            
            if (plot_zscore)
                model_hit2_lst2 = norminv(model_hit2_lst);
                model_fa2_lst2  = norminv(model_fa2_lst);
            end

            plot(model_fa2_lst2, model_hit2_lst2, '-', ...
                'Color', my_col1, 'LineWidth', 3);
            
            cncb_plot_struct.model_type2_roc = [model_fa2_lst; model_hit2_lst];
        end
        
        cncb_plot_struct.human_type2_roc = [human_fa2_lst; human_hit2_lst];

        axis('square');
        if (plot_zscore)
            axis([-3, 3, -3, 3]);
            xticks(-3:1:3);
            yticks(-3:1:3);
        else
            axis([0, 1, 0, 1]);
            xticks(0.0:0.2:1.0);
            yticks(0.0:0.2:1.0);
        end
        
        if (stim_nb == 2)
            xlabel('p(high conf | incorrect)', 'FontSize', 26);
            if (mm == 1)
                ylabel('p(high conf | correct)', 'FontSize', 26);
            end
        elseif (mm == stim_nb)
            xlabel('p(high conf | incorrect)', 'FontSize', 26);
            ylabel('p(high conf | correct)', 'FontSize', 26);            
        end
        
    end
    
end


% -> ********************************************** <-
% -> nested functions

% -> convert SRC format to ROC format: (stim, conf_level, t2_hit, t2_fa)
function data_ROC = src_to_roc(data_SRC)

    lcl_col_stim = 1;
    lcl_col_resp = 2;
    lcl_col_conf_levl = 3;
    lcl_col_conf_prob = 4;

    lcl_stim_lst = unique(data_SRC(:, lcl_col_stim));
    lcl_stim_nb = length(lcl_stim_lst);

    lcl_resp_lst = unique(data_SRC(:, lcl_col_resp));

    rtg_mat = data_SRC;
    
    nb_conds = size(data_SRC, 1);
    
    my_rating_inds = unique(data_SRC(:, lcl_col_conf_levl));
    my_conf_levels_nb = length(my_rating_inds);
    my_rating_bnds_nb = my_conf_levels_nb - 1;

    data_ROC = NaN(lcl_stim_nb*my_rating_bnds_nb, 4);
    rng_all = 1:nb_conds;

    for lcl_ss = 1:lcl_stim_nb
        lcl_stim_val = lcl_stim_lst(lcl_ss);

        % -> prob of each response
        lcl_left_inds = ((data_SRC(:,lcl_col_stim) == lcl_stim_val) & ...
            (data_SRC(:,lcl_col_resp) == lcl_resp_lst(1)));
        lcl_righ_inds = ((data_SRC(:,lcl_col_stim) == lcl_stim_val) & ...
            (data_SRC(:,lcl_col_resp) == lcl_resp_lst(2)));

        lcl_resp_L = sum(data_SRC(lcl_left_inds, lcl_col_conf_prob));
        lcl_resp_R = sum(data_SRC(lcl_righ_inds, lcl_col_conf_prob));

        rtg_mat(lcl_left_inds, lcl_col_conf_prob) = ...
            data_SRC(lcl_left_inds, lcl_col_conf_prob) / lcl_resp_L;
        rtg_mat(lcl_righ_inds, lcl_col_conf_prob) = ...
            data_SRC(lcl_righ_inds, lcl_col_conf_prob) / lcl_resp_R;


        if (lcl_ss <= (lcl_stim_nb/2))
            resp_corr = lcl_resp_lst(1);
            resp_inco = lcl_resp_lst(2);
        else
            resp_corr = lcl_resp_lst(2);
            resp_inco = lcl_resp_lst(1);
        end

        lcl_inco_inds = ((rtg_mat(:,lcl_col_stim) == lcl_stim_val) & ...
            (rtg_mat(:,lcl_col_resp) == resp_inco));
        lcl_cons_inds = ((rtg_mat(:,lcl_col_stim) == lcl_stim_val) & ...
            (rtg_mat(:,lcl_col_resp) == resp_corr));
        rng_inco = rng_all(lcl_inco_inds);
        rng_cons = rng_all(lcl_cons_inds);
        for lcl_rr = 1:my_rating_bnds_nb
            lcl_type2_cr  = sum(rtg_mat(rng_inco(1:lcl_rr), lcl_col_conf_prob));
            lcl_type2_miss = sum(rtg_mat(rng_cons(1:lcl_rr), lcl_col_conf_prob));
            lcl_type2_fa  = 1 - lcl_type2_cr;
            lcl_type2_hit = 1 - lcl_type2_miss;

            data_ROC((lcl_ss-1)*my_rating_bnds_nb + lcl_rr, :) = ...
                [lcl_stim_val, my_rating_inds(lcl_rr), lcl_type2_hit, lcl_type2_fa];
        end

    end
end
    
end

% -> THE END <-
