function [n0, nC, is_occasion, frame_offset] = ...
    getPDCCH0MonitoringOccasions(idx, iSSB, SCS, pattern, ...
    N_sym_CORESET, SFN)
%GETPDCCH0MONITORINGOCCASIONS Returns the Type0-PDCCH monitoring occasions
%based on 3GPP 38.213 13 and Tables from 13-11 to 13-15
% Inputs:
%   idx                 : a number representing the 4 LSBs of the
%   pdcch-ConfigSIB1
%   iSSB                : a number representing the SSB index
%   SCS                 : a number representing the SCS in Hz
%   pattern             : a number representing the multiplexing pattern
%   N_CORESET_symbols   : a number representing the number of symbols in
%   CORESET#0
%   SFN                 : a number representing the system frame number
% Outputs:
%   n0              : a number representing the slot number of the
%   CORESET#0 location
%   nC              : a number representing the OFDM symbol start index of
%   the CORESET#0 location
%   is_occasion     : a boolean that indicates if the frame is a PDCCH
%   monitoring occasion
%   fram_offset     : a number representing the frame offset of the
%   CPRESET#0 location

    SCS_SSB = SCS(1);
    SCS_common = SCS(2);
    
    if pattern == 1      
        if (SCS_SSB == 15 || SCS_SSB == 30) % Frequency Range 1
            if mod(iSSB, 2) == 0
                nC = 0;
            else
                nC = N_sym_CORESET;
            end
            
            tab = [ 0,    1, 2,    3, 4,    5, 6,    7, 8, 9, 10, 11, 12, 13, 14, 15;...
                    0,    0, 2,    2, 5,    5, 7,    7, 0, 5,  0,  0,  2,  2,  5,  5;...
                    1,    2, 1,    2, 1,    2, 1,    2, 1, 1,  1,  1,  1,  1,  1,  1;...
                    1,  1/2, 1,  1/2, 1,  1/2, 1,  1/2, 2, 2,  1,  1,  1,  1,  1,  1;...
                    0, nC, 0, nC, 0, nC, 0, nC, 0, 0,  1,  2,  1,  2,  1,  2];
        else % Frequency Range 2
            if mod(iSSB,2) == 0 % Even SSB index
                firs_sym_idx_1 = 0;
                firs_sym_idx_2 = 0;
            else
                firs_sym_idx_1 = 7;
                firs_sym_idx_2= N_sym_CORESET;
            end
            
            tab = [ 0,     1,   2,     3, 4,     5,     6,     7,     8,    9,    10,    11, 12, 13, 14, 15;...
                    0,     0, 2.5,   2.5, 5,     5,     0,   2.5,     5,  7.5,   7.5,   7.5,  0,  5, NaN(1,2);...
                    1,     2,   1,     2, 1,     2,     2,     2,     2,    1,     2,     2,  1,  1, NaN(1,2);...
                    1,   1/2,   1,   1/2, 1,   1/2,   1/2,   1/2,   1/2,    1,   1/2,   1/2,  2,  2, NaN(1,2);...
                    0, firs_sym_idx_1,   0, firs_sym_idx_1, 0, firs_sym_idx_1, firs_sym_idx_2, firs_sym_idx_2, firs_sym_idx_2, 0, firs_sym_idx_1, firs_sym_idx_2,  0,  0, NaN(1,2)];
        end
        
        % Extract info from table
        t = tab(:,idx==tab(1,:));
        O = t(2);
        M = t(4);
        nC = t(5);
        
        mu = log2(SCS_common/15e3);
        slots_per_frame = 10*2^mu;
        
        % Slot of monitoring occasion
        slot = O*2^mu + floor(iSSB*M);
        frame_offset = floor(slot/slots_per_frame);
        n0 = mod(slot, slots_per_frame);
        
        if xor(mod(frame_offset, 2), mod(SFN, 2))
            is_occasion = false;
        else
            is_occasion = true;
        end
        
    else
        error("Not implemented");
    end
end