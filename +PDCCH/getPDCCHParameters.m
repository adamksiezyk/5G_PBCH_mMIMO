function pdcch = getPDCCHParameters(SCS_SSB, min_channel_BW, MIB, iSSB, ...
    cellid)
%GETPDCCHPARAMETERS Returns a struct of PDCCH parameters based on 3GPP
%38.213 section 13 and 38.211 section 7.3.2.3
% Inputs:
%   SCS_SSB         : a number representing the SSB SCS in Hz
%   min_channel_BW  : a number represetning the minimum CORESET#0 channel
%   bandwidth in Hz
%   MIB             : a MIB class representing the Master Information Block
%   iSSB            : a number representing the SSB index
%   cellid          : a number representing the Cell ID
% Outputs:
%   pdcch   : a PDCCHParameters object representing the PDCCH parameters

    msb_idx = floor(MIB.PDCCHConfigSIB1/16);    % 4 MSB of PDCCHConfigSIB1
    lsb_idx = mod(MIB.PDCCHConfigSIB1, 16);     % 4 LSB of PDCCHConfigSIB1
    SCS_common = MIB.SubcarrierSpacingCommon * 1e3;
    SCS_pair = [SCS_SSB, SCS_common];
    N_subframes_per_frame = 10;
    N_slots_per_frame = SCS_common / (15*1e3) * N_subframes_per_frame;

    % Get CORESET#0 information
    [N_RB_CORESET, N_sym_CORESET, CORESET_RB_offset, pattern] = ...
        PDCCH.getCORESET0Resources(msb_idx, SCS_pair, min_channel_BW, ...
        MIB.kSSB);
    
    % Find Type0-PDCCH monitoring occasions. 3GPP 28.213 13
    [n0, nC, is_occasion, frame_offset] = ...
        PDCCH.getPDCCH0MonitoringOccasions(lsb_idx, iSSB, SCS_pair, ...
        pattern, N_sym_CORESET, MIB.NFrame);
    
    % Configure Search Space
    search_space = PDCCH.SearchSpace;
    search_space.StartSymbolWithinSlot = nC;
    if pattern == 1
        search_space.SlotPeriod = 2*N_slots_per_frame;
        search_space.SlotOffset = n0;
        search_space.Duration = 2;
    else % patterns 2 and 3
        search_space.SlotPeriod = 2*N_slots_per_frame;
        search_space.SlotOffset = n0;
        search_space.Duration = 1;
    end
    search_space.NumCandidates = [0 0 8 4 1]; % TS 38.213 Table 10.1-1
    
    % PDCCH parameters
    pdcch = PDCCH.PDCCHParameters;
    pdcch.N_RB_start = 0;
    pdcch.N_RB = N_RB_CORESET;
    pdcch.SearchSpace = search_space;
    pdcch.RNTI = 0;                     % TS 38.211 Section 7.3.2.3
    pdcch.DMRS_scrambling_ID = cellid;  % TS 38.211 Section 7.3.2.3
end

