classdef PDCCHParameters
    %PDCCHPARAMETERS Represents the PDCCH parameters
    %   3GPP 38.213 section 13 and 38.211 section 7.3.2.3
    
    properties
        CORESET {mustBeA(CORESET, "PDCCH.CORESET")} = PDCCH.CORESET;        % Represents the Control Resource Set
        SearchSpace {mustBeA(SearchSpace, "PDCCH.SearchSpace")} = ...
            PDCCH.SearchSpace                                               % Represents the Search Space
        RNTI {mustBeInteger} = 0                                            % Represents the RNTI TS 38.211 Section 7.3.2.3
        DMRS_scrambling_ID {mustBeInRange(DMRS_scrambling_ID, 0, 65535)}    % Represents the DMRS scrambling ID n_ID TS 38.211 Section 7.3.2.3
    end
end
