function [mon_slots, mon_symbols, mon_subcarriers] = ...
    getMonitoringOccasionsResources(lsb_idx, iSSB, SCS_pair, N_frame, ...
    pattern, N_sym_CORESET, N_RB_CORESET, N_RB_CORESET_min, ...
    CORESET_RB_offset, N_symbols, signal_info, carrier)
%GETMONITORINGOCCASIONSRESOURCES Returns the Type0-PDCCH monitoring
%occasions based on 3GPP 38.213 13 and Tables from 13-11 to 13-15
% Inputs:
%   lsb_idx             : a number representing the 4 LSBs of the
%   pdcch-ConfigSIB1
%   iSSB                : a number representing the SSB index
%   SCS_pair            : two numbers representing the SCS and SCS_common 
%   in Hz
%   N_frame             : a number representing the System Frame Number
%   pattern             : a number representing the multiplexing pattern
%   N_CORESET_symbols   : a number representing the number of symbols in
%   CORESET#0
%   N_RB_CORESET        : a number representing the number of Resource
%   Blocks in CORESET#0
%   N_RB_CORESET_min    : a number representing the minimum number of
%   Resource Blocks to cover CORESET#0
%   CORESET_RB_offset   : a number representing the Resource Block offset
%   from the start of the resource grid
%   N_symbols           : a number representing the number of symbols in
%   the signal
%   signal_info         : a SignlInfo class
%   carrier             : a Carrier class
% Outputs:
%   mon_slots           : a vector representing the monitoring slots
%   mon_symbols         : a vector representing the monitoring symbols
%   mon_subcarriers     : a vector representing the monitoring subcarriers


    [n0, nC, is_occasion, frame_offset] = ...
        PDCCH.getPDCCH0MonitoringOccasions(lsb_idx, iSSB, SCS_pair, ...
        pattern, N_sym_CORESET, N_frame);
    
    N_slots = ceil(N_symbols / carrier.SymbolsPerSlot);
    N_mon_slots = 2;
    mon_slots = n0 + (0:N_mon_slots-1)' + (0:2*carrier.SlotsPerFrame:(N_slots-n0-1));
    mon_slots = mon_slots(:)';
    mon_symbols = mon_slots*carrier.SymbolsPerSlot + (1:carrier.SymbolsPerSlot)';
    mon_symbols = mon_symbols(:)';
    mon_symbols(mon_symbols > N_symbols) = [];
    
    mon_subcarriers = (N_RB_CORESET_min - 20*SCS_pair(1) / ...
        SCS_pair(2))*signal_info.N_subcarriers_per_RB/2 + CORESET_RB_offset * ...
        signal_info.N_subcarriers_per_RB + [1:N_RB_CORESET*12];
end