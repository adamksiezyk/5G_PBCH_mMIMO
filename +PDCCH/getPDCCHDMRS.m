function DMRS = getPDCCHDMRS(carrier, pdcch, aggregation_level_index)
% Returns the PDCCH DM-RS symbols for all candidates according to 3GPP
%38.211 7.4.1.3.
% Inputs:
%   carrier                     : a Carrier class
%   pdcch                       : a PDCCH.PDCCHParameters class
%   aggregation_level_index     : a number representing the aggregation
%   level index

    aggregation_levels = [1 2 4 8 16];  % 3GPP 38.211 Table 7.3.2.1-1: Supported PDCCH aggregation levels
    DMRS_indices_per_RB = [1 5 9] + 1;  % DM-RS, 3GPP 38.211 7.4.1.3.2
    N_DMRS_per_RB = length(DMRS_indices_per_RB);
    
    pdcch_RB_indices_per_candidate = PDCCH.getPDCCHRBIndices(...
        pdcch.CORESET, pdcch.SearchSpace, aggregation_level_index);
    DMRS = zeros(aggregation_levels(aggregation_level_index) * ...
        pdcch.CORESET.N_REGs_per_CCE * N_DMRS_per_RB, ...
        pdcch.SearchSpace.NumCandidates(aggregation_level_index));
    for j = 1:pdcch.SearchSpace.NumCandidates(aggregation_level_index)
        candidate_RB_indices = pdcch_RB_indices_per_candidate(:,:,j);
        DMRS(:,j) = createDMRS(carrier, pdcch, candidate_RB_indices);
    end
end

function DMRS = createDMRS(carrier, pdcch, pdcch_RB_indices)
    modulation_order = 4;   % QPSK
    N_slot = mod(carrier.NSlot, carrier.SlotsPerFrame);
    N_sym = pdcch.CORESET.N_sym;
    N_start_sym = pdcch.SearchSpace.StartSymbolWithinSlot;
    pdcch_PRB_indices = pdcch_RB_indices(:,1) - ...
        N_start_sym * carrier.N_RB;
    N_PRB = max(pdcch_PRB_indices);
    symbols_per_RB = 3*2;  % 3REs per RB, 2 bits per DM-RS symbol
    N_DMRS = symbols_per_RB * N_PRB;
    DMRS = complex(zeros(3*length(pdcch_PRB_indices), N_sym));
    for sym = 1:N_sym
        % Create binary DMRS sequence
        l = sym - 1  + N_start_sym;
        seq = createDMRSSequence(l, N_slot, carrier.NCellID, N_DMRS);
        seq_RB = reshape(seq, symbols_per_RB, []);
        sequence = seq_RB(:, pdcch_PRB_indices);
        % Modulated sequence
        DMRS(:, sym) = utils.modulateQAM(sequence(:), modulation_order);
    end
    DMRS = DMRS(:);
end

function c = createDMRSSequence(l, n_s, N_cell_id, M_PN)
% 3GPP 38.211 7.4.1.3.1
    c_init = mod(2^17*(14*n_s+l+1)*(2*N_cell_id+1)+2*N_cell_id, 2^31);
    c = utils.generatePRBS(M_PN, c_init);
end