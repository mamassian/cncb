% convert confidence data from CNCB to ML12 format
%
% ML12 format used in:
%   Maniscalco, B., & Lau, H. C. (2012). A signal detection theoretic 
%      approach for estimating metacognitive sensitivity from confidence 
%      ratings. Consciousness and Cognition, 21(1), 422-430. 
%      http://doi.org/10.1016/j.concog.2011.09.021
%
% 26-AUG-2021 - pascal mamassian
% 10-MAR-2024 - pm: cleaned up

function  [nR_S1, nR_S2] = CNCBtoML12(data_SRC)

    col_stim = 1;
    col_resp = 2;
    col_conf_levl = 3;
    col_conf_prob = 4;

    [stim_lst, ~, stim_ic] = unique(data_SRC(:, col_stim));
    stim_nb = length(stim_lst);
    if (stim_nb > 2)
        fprintf('ERROR in CNCBtoML12: more than 2 stimuli\n');
        nR_S1 = NaN;
        nR_S2 = NaN;
        return;
    end

    [~, ~, resp_ic] = unique(data_SRC(:, col_resp));

    [rating_inds, ~, rating_ic] = unique(data_SRC(:, col_conf_levl));
    conf_levels_nb = length(rating_inds);

    % -> use cell padding technique advised in 'type2_SDT_MLE'
    cellpadding = 1 / (2*conf_levels_nb);
    
    % -> convert to trial count format...    
    nR_S1 = NaN(1, 2*conf_levels_nb);
    nR_S2 = NaN(1, 2*conf_levels_nb);

    % -> get tallies of "S1" rating responses for S1 and S2 stim
    for ii = 1:conf_levels_nb
        nR_S1(ii) = sum(data_SRC((stim_ic==1) & (resp_ic==1) & ...
            (rating_ic==(conf_levels_nb+1-ii)), col_conf_prob));
        nR_S2(ii) = sum(data_SRC((stim_ic==2) & (resp_ic==1) & ...
            (rating_ic==(conf_levels_nb+1-ii)), col_conf_prob));
    end
    
    % -> get tallies of "S2" rating responses for S1 and S2 stim
    for ii = 1:conf_levels_nb
        nR_S1(ii+conf_levels_nb) = sum(data_SRC((stim_ic==1) & ...
            (resp_ic==2) & (rating_ic==ii), col_conf_prob));
        nR_S2(ii+conf_levels_nb) = sum(data_SRC((stim_ic==2) & ...
            (resp_ic==2) & (rating_ic==ii), col_conf_prob));
    end
    
    if any(nR_S1==0) || any(nR_S2==0)
        nR_S1 = nR_S1 + cellpadding;
        nR_S2 = nR_S2 + cellpadding;
    end

end
