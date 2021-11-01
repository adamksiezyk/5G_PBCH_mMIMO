function [all_symbol_indices, all_DMRS_indices] = getPDCCHResources(...
    carrier, pdcch)
%GETPDCCHRESOURCES Summary of this function goes here
%   Detailed explanation goes here   
    if isOccasion(carrier, pdcch)
        [all_symbol_indices, all_DMRS_indices] = getAllPDCCHREIndices(...
            pdcch.CORESET, pdcch.SearchSpace);
    end
end

function is_occasion = isOccasion(carrier, pdcch)
% Returns true if a PDCCH monitoring occasion exists in current slot.
% 3GPP 38.213 10.1
    ss = pdcch.SearchSpace;
    is_occasion = (mod(carrier.NFrame * carrier.SlotsPerFrame + ...
        carrier.NSlot - ss.SlotOffset, ss.SlotPeriod) == 0);
end

function [symbol_RE_indices, DMRS_RE_indices] = getAllPDCCHREIndices(...
    coreset, search_space)
% Returns all PDCCH symbol and DMRS indices
    % Get all candidates RE indices for PDCCH and PDCCH DM-RS
    aggregation_levels = [1 2 4 8 16];              % 3GPP 38.211 Table 7.3.2.1-1: Supported PDCCH aggregation levels
    N_AL = length(aggregation_levels);
    symbol_indices_per_RB = [0 2:4 6:8 10:11] + 1;  % Exclude DM-RS
    N_symbol_per_RB = length(symbol_indices_per_RB);
    DMRS_indices_per_RB = [1 5 9] + 1;              % DM-RS, 3GPP 38.211 7.4.1.3.2
    N_DMRS_per_RB = length(DMRS_indices_per_RB);
    symbol_RE_indices = cell(N_AL, 1);
    DMRS_RE_indices = cell(N_AL, 1);

    % Loop over aggregation levels
    for i = 1:N_AL
        pdcch_RB_indices_per_candidate = PDCCH.getPDCCHRBIndices(...
            coreset, search_space, i);

        % Loop over candidates
        N_candidates = search_space.NumCandidates(i);
        tmp_symbol_indices = zeros(aggregation_levels(i) * ...
            coreset.N_REGs_per_CCE * N_symbol_per_RB, ...
            N_candidates);
        tmp_DMRS_indices = zeros(aggregation_levels(i) * ...
            coreset.N_REGs_per_CCE * N_DMRS_per_RB, ...
            N_candidates);
        for j = 1:N_candidates
            candidate_RB_indices = reshape(...
                pdcch_RB_indices_per_candidate(:,:,j), [], 1);
            % Expand to RE level from RB indexes (frequency-first
            % within an RB)
            for idx = 1:length(candidate_RB_indices)
                % PDCCH RE indices
                tmp_symbol_indices((idx-1)*N_symbol_per_RB+(1:N_symbol_per_RB),j) = ...
                    ((candidate_RB_indices(idx)-1)*12+symbol_indices_per_RB).';
                % PDCCH DM-RS RE indices
                tmp_DMRS_indices((idx-1)*N_DMRS_per_RB+(1:N_DMRS_per_RB),j) = ...
                    ((candidate_RB_indices(idx)-1)*12+DMRS_indices_per_RB).';
            end
        end
        symbol_RE_indices{i} = tmp_symbol_indices;
        DMRS_RE_indices{i} = tmp_DMRS_indices;
    end
end