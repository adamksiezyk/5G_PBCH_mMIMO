function RB_indices = getPDCCHRBIndices(...
    coreset, search_space, aggregation_level_index)
% Returns the PDCCH RB indices (per aggregation level) found in the
% candidates. 3GPP 38.213 10.1, 38.211.
% Inputs:
%   coreset                     : a PDCCH.CORESET class
%   search_space                : a PDCCH.SearchSpace class
%   aggregation_level_index     : a number representing the aggregation
%   level index
% Outputs:
%   RB_indices  : a matrix representing the Resource Block indices for each
%   candidate

    aggregation_levels = [1 2 4 8 16];  % 3GPP 38.211 Table 7.3.2.1-1: Supported PDCCH aggregation levels
    REG_bundle_mapping = coreset.getREGBundleIndices();
    
    N_CCEs = aggregation_levels(aggregation_level_index);
    RB_indices = zeros(...
        N_CCEs * coreset.N_REGs_per_CCE / coreset.N_sym, ...
        coreset.N_sym, search_space.NumCandidates(aggregation_level_index));
    candidates = getCCEIndices(search_space, coreset, aggregation_level_index);
    N_candidates = size(candidates, 2);

    for j = 1:N_candidates              % loop over candidates
        CCE_indices = candidates(:, j);

        % Map REG bundles to RBs
        f = coreset.getCCEMapping();
        REG_bundles = REG_bundle_mapping(:,f+1);
        % Map CCEs to REG bundles -> RBs
        CCEs = reshape(REG_bundles, coreset.N_REGs_per_CCE, []);
        CCE_RB_indices = CCEs(:, CCE_indices);

        % Map the RBs to symbols in slot in a frequency first order
        RB_indices(:, :, j) = sort(reshape(CCE_RB_indices, ...
            coreset.N_sym, [])');
    end
end

function CCE_indices = getCCEIndices(...
    search_space, coreset, aggregation_level_idx)
% Returns CCE indices for common search space for each candidate
% according to 3GPP 38.213 10.1 and 3GPP 38.211 7.3.2
    n_CI = 0;   % 0 for common search space (3GPP 38.213 10.1)
    N_CCE = coreset.N_sym * coreset.N_RB / coreset.N_REGs_per_CCE;
    Y_p = 0;     % 0 for common search space (3GPP 38.213 10.1)

    % 3GPP 38.211 7.3.2
    aggregation_levels = [1 2 4 8 16];
    M_p = search_space.NumCandidates(aggregation_level_idx);
    L = aggregation_levels(aggregation_level_idx);

    CCE_indices = zeros(L, M_p);
    for ms = 0:M_p-1
        for idx = 1:L
            CCE_indices(idx, ms+1) = L*( mod(Y_p + ...
                floor(ms*N_CCE/(L*M_p)) + n_CI, ...
                floor(N_CCE/L)) ) + idx;
        end
    end
end