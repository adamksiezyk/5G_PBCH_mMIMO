function timing_offset = estimateTimingOffset(sample_rate, N_sym_long, ...
    N_sym, SSB_case, frequency, offset, i_SSB)
%ESTIMATETIMINGOFFSET Summary of this function goes here
% Inputs:
%   sample_rate     : a number representing the sample rate
%   N_sym_long      : a number representing the long OFDM symbol samples
%   amount
%   N_sym           : a number representing the OFDM symbol samples amount
%   SSB_case        : a string representing the SSB case
%   frequency       : a number representing the center frequency
%   offset          : a number representing the PSS timing index
%   i_SSB           : a number representing the SSB index
% Outputs:
%   timing_offset   : a number representing the frame timing offset

    % Adjust timing offset to the start of the SS block. This step removes
    % the extra offset introduced in the reference grid during PSS search,
    % which contained the PSS in the second OFDM symbol.
    offset = offset + N_sym_long;
    
    % Timing offset is adjusted so that the received grid starts at the
    % frame head i.e. adjust the timing offset for the difference between
    % the first symbol of the strongest SSB, and the start of the frame
    SSB_start_symbols = utils.getSSBStartSymbols(SSB_case, frequency);
    SSB_start_symbol = SSB_start_symbols(i_SSB+1);
    
    % Adjust for whole slots
    symbols_per_slot = 14;
    slot_offset = floor(SSB_start_symbol/symbols_per_slot);
    samples_per_slot = N_sym_long + 13*N_sym;
    timing_offset = offset - (slot_offset*samples_per_slot);
    
    % Adjust for remaining OFDM symbols and round offset if not integer
    symbol_offset = mod(SSB_start_symbol, symbols_per_slot);
    symbol_lengths_in_slot = [N_sym_long, N_sym*ones(1, 13)];
    timing_offset = round(timing_offset - sum(...
        symbol_lengths_in_slot(1:symbol_offset)));
    
    % Apply modulo 20 ms periodicity in case offset is negative
    if timing_offset < 0
        timing_offset = mod(timing_offset, sample_rate*20e-3);
    end
end