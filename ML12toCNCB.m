% convert confidence data from ML12 to CNCB format
%
% ML12 format used in:
%   Maniscalco, B., & Lau, H. C. (2012). A signal detection theoretic 
%      approach for estimating metacognitive sensitivity from confidence 
%      ratings. Consciousness and Cognition, 21(1), 422-430. 
%      http://doi.org/10.1016/j.concog.2011.09.021
%
% 31-MAR-2024 - pascal mamassian

function  data_SRC = ML12toCNCB(nR_S1, nR_S2)

    col_stim = 1;
    col_resp = 2;
    col_conf_levl = 3;
    col_conf_prob = 4;

    stim_lst = [-1, 1];
    stim_nb = length(stim_lst);

    resp_lst = [0, 1];
    resp_nb = length(resp_lst);

    % -> infer the number of confidence levels
    conf_levels_nb = (length(nR_S1) + length(nR_S1)) / (stim_nb * resp_nb);
    if (floor(conf_levels_nb) ~= conf_levels_nb)
        fprintf('ERROR in ML12toCNCB: input not well formatted\n');
        data_SRC = NaN;
        return;
    end
    conf_levels_lst = 1:conf_levels_nb;

    conds_nb = stim_nb * resp_nb * conf_levels_nb;
    data_SRC = NaN(conds_nb, 4);

    ml12_ind = 1;
    for rr = 1:resp_nb
        for cc = 1:conf_levels_nb
            for ss = 1:stim_nb

                if (rr == 1)
                    conf_ind = conf_levels_nb - cc + 1;
                else
                    conf_ind = cc;
                end
                cncb_ind = (rr-1)*conf_levels_nb + conf_ind;

                if (ss == 1)
                    data_set = nR_S1;
                else
                    data_set = nR_S2;
                    cncb_ind = cncb_ind + resp_nb*conf_levels_nb;
                end

                conf_count = data_set(ml12_ind);

                data_SRC(cncb_ind, col_stim) = stim_lst(ss);
                data_SRC(cncb_ind, col_resp) = resp_lst(rr);
                data_SRC(cncb_ind, col_conf_levl) = conf_levels_lst(conf_ind);
                data_SRC(cncb_ind, col_conf_prob) = conf_count;

            end
            ml12_ind = ml12_ind + 1;
        end
    end

end
