% CNCB toolbox(Confidence Noise Confidence Boost) -- v0.1
%
% cncb_group
%   This function groups trials from a confidence ratings experiment.
%   The idea is to replace the raw data with group of trials that have the
%   same properties (same stimuli and same responses).
%
% INPUT:
%   'raw_data': matrix of:
%       1st col: stimulus intensity (e.g. -1.0 for 'B' and 1.0 for 'A')
%       2nd col: perceptual decision (1 = 'A',  0 = 'B')
%       3rd col: confidence rating value (numeric)
%
% OPTIONAL PARAMETERS:
%   'confidence_is_continuous': are confidence judgments on a continuous scale?
%                (true / false). Default is 'false' (i.e. discrete ratings).
%
%   'confidence_cont_range': if continuous, range of confidence ratings
%                [min, max]. Default is [0, 1].
%
%   'confidence_cont_nb_levels': if continuous, the number of points used to
%                estimate the confidence function. Default is 100.
%
%   'confidence_disc_labels': if discrete, a vector of labels for the confidence
%                levels. Only numeric labels are supported so far.
%
% OUTPUT:
%   'grouped_data': grouped data per (stim, resp, conf):
%       1st col: stimulus intensity
%       2nd col: perceptual decision
%       3rd col: confidence level
%       4th col: count of confidence choices per (stim, resp, conf)
%
%
%
% EXAMPLES OF USE:
%   grouped_data = cncb_group(raw_data);
%   grouped_data = cncb_group(raw_data, ...
%        'confidence_is_continuous', true, ...
%        'confidence_half_scale', false, ...
%        'confidence_cont_range', [0, 100], ...
%        'confidence_cont_nb_levels', 100);
%   grouped_data = cncb_group(raw_data, ...
%        'confidence_is_continuous', false, ...
%        'confidence_disc_labels', [1, 2, 3, 4]);
%   
%
% 26-AUG-2021 - pascal mamassian: created
% 16-FEB-2022 - pm: added continuous ratings
% 25-OCT-2023 - pm: fixed continuous ratings
% 08-MAR-2024 - pm: group continuous ratings by quantiles
% 31-MAR-2024 - pm: add warning if only one confidence level


function grouped_data_SRC = cncb_group(raw_data_SRC, varargin)

    % -> default optional arguments
    dflt_is_continuous        = false;  % continuous (T) or discrete (F) rating scale
    dflt_conf_cont_half_scale = false;
    dflt_conf_cont_range      = [0, 1];
    dflt_conf_cont_nb_levels  = 100;    % used to fit continuous data
    dflt_conf_disc_labels     = [];

    % -> parse all arguments
    ip = inputParser;
    ip.StructExpand = false;
    addRequired(ip, 'raw_data_SRC', @isnumeric);
    addParameter(ip, 'confidence_is_continuous', dflt_is_continuous, @islogical);
    addParameter(ip, 'confidence_half_scale', dflt_conf_cont_half_scale, @islogical);
    addParameter(ip, 'confidence_cont_range', dflt_conf_cont_range, @isnumeric);
    addParameter(ip, 'confidence_cont_nb_levels', dflt_conf_cont_nb_levels, @isnumeric);
    addParameter(ip, 'confidence_disc_labels', dflt_conf_disc_labels, @isnumeric);

    parse(ip, raw_data_SRC, varargin{:});
    conf_continuous = ip.Results.confidence_is_continuous;
    conf_half_scale = ip.Results.confidence_half_scale;
    conf_range = ip.Results.confidence_cont_range;
    conf_cont_nb_levels = ip.Results.confidence_cont_nb_levels;
    conf_disc_labels = ip.Results.confidence_disc_labels;

    col_stim = 1;
    col_resp = 2;
    col_conf_levl = 3;
    col_conf_prob = 4;  % prob(conf|stim) or count
    
    % -> if raw data have only 3 columns, add a 4th one
    if (size(raw_data_SRC, 2) < col_conf_prob)
        nb_lines = size(raw_data_SRC, 1);
        raw_data_SRC = [raw_data_SRC, ones(nb_lines, 1)];
    end

    [stim_lst, ~, stim_ic] = unique(raw_data_SRC(:, col_stim));
    stim_nb = length(stim_lst);

    [resp_lst, ~, resp_ic] = unique(raw_data_SRC(:, col_resp));
    resp_nb = length(resp_lst);

    if (conf_continuous)
        % -> divide confidence judgments in quantiles (bins containing an equal
        %    number of trials)
        if (sum((raw_data_SRC(:, col_conf_prob) ~= 1)) ~= 0)

            % -> if the last column is not all '1', assume grouping is already done
            [conf_lst, ~, ~] = unique(raw_data_SRC(:, col_conf_levl));
            conf_bnd_labels = conf_lst';
            conf_bnd_edges = [conf_range(1), conf_bnd_labels];

        else
            conf_vals = raw_data_SRC(:, col_conf_levl);

            if (conf_half_scale)
                % -> if half-scale, there might be a large accumulation of
                %    confidence values on lower-bound, so perform quantile
                %    on the rest of the distribution
                conf_vals = conf_vals(conf_vals ~= min(conf_vals));
            end

            conf_bnd_edges = quantile(conf_vals, conf_cont_nb_levels - 1);
            if (conf_bnd_edges(end) > conf_range(end))
                fprintf('Error: confidence range is not large enough\n');
            end
            conf_bnd_edges = [conf_range(1), conf_bnd_edges, conf_range(end)];
            conf_bnd_labels = conf_bnd_edges(2:end);

            % -> check that all edges are different
            aa = (diff(conf_bnd_labels) == 0);
            my_count = 1;
            while ((sum(aa) ~= 0) && (my_count < 100))
                first_ind = find(aa, 1, 'first');
                if (first_ind ~= 1)
                    conf_bnd_labels(first_ind) = ...
                        (conf_bnd_labels(first_ind) + ...
                         conf_bnd_labels(first_ind-1)) / 2;
                end
                aa = (diff(conf_bnd_labels) == 0);
                my_count = my_count + 1;
            end
        end

    else
        if (isempty(conf_disc_labels))
            [conf_bnd_labels, ~, rating_ic] = unique(raw_data_SRC(:, col_conf_levl));
        else
            conf_bnd_labels = conf_disc_labels;
            [~, rating_ic] = ismember(raw_data_SRC(:, col_conf_levl), conf_bnd_labels);
        end
    end
    conf_levels_nb = length(conf_bnd_labels);

    if (conf_levels_nb == 1)
        fprintf('Warning: only 1 confidence level in dataset\n');
    end

    nb_conds = stim_nb * resp_nb * conf_levels_nb;

    grouped_data_SRC = zeros(nb_conds, 4);
    
    kk = 0;
    for ss = 1:stim_nb
        stim_val = stim_lst(ss);
        
        for rr = 1:resp_nb
            resp_val = resp_lst(rr);

            for cc = 1:conf_levels_nb
                conf_val = conf_bnd_labels(cc);

                % -> search for number of confidence trials in this category
                if (conf_continuous)
                    conf_low = conf_bnd_edges(cc);
                    conf_hig = conf_bnd_edges(cc+1);
                    ratings_inds = ((raw_data_SRC(:, col_conf_levl) > conf_low) & ...
                        (raw_data_SRC(:, col_conf_levl) <= conf_hig));
                    if (cc == 1)
                        ratings_inds = (ratings_inds | ...
                            (raw_data_SRC(:, col_conf_levl) == conf_low));
                    end
                    conf_num = sum(raw_data_SRC((stim_ic==ss) & (resp_ic==rr) & ...
                        ratings_inds, col_conf_prob));
                else
                    conf_num = sum(raw_data_SRC((stim_ic==ss) & (resp_ic==rr) & ...
                        (rating_ic==cc), col_conf_prob));
                end

                kk = kk + 1;

                grouped_data_SRC(kk, col_stim) = stim_val;
                grouped_data_SRC(kk, col_resp) = resp_val;
                grouped_data_SRC(kk, col_conf_levl) = conf_val;
                grouped_data_SRC(kk, col_conf_prob) = conf_num;
            
            end
        end
    end
end

% -> THE END <-
