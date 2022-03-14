function [fractional_CFO, CP_corr, indv_CFO] = estimateFractionalCFO(waveform, ...
    sample_rate, scs)
%ESTIMATEFRACTIONALCFO Estimated the fractional center frequency offset
% Inputs:
%   waveform                : a column vector representing the signal
%   sample_rate             : a number representing the sample rate in Hz
%   scs                     : a number representing the subcarrier spacing
%   in Hz
%   ssb_frequency_offset    : a number representing the SSB frequency
%   offset in relation to the center frequecy in Hz
% Outputs:
%   fractional_CFO  : a number representing the estimated fractional center
%   frequency offset
%   CP_corr         : a vector representing the calculated cyclic prefix
%   based correlation
%   indv_CFO        : a vector representing the individual fractional CFO
%   values

    % Signal autocorelation
    N_CP = utils.getCPLength(scs, sample_rate);
    N_FFT = 2048;
    N_sym = N_FFT + N_CP;
    N_corr = sample_rate * 10e-3;   % 1 5G frame
    CP_corr = zeros(1, N_corr);
    for n = 1:N_corr
        n1 = n:n+N_CP-1; 
        n2 = n1 + N_FFT;
        CP_corr(n) = abs(sum( waveform(n1) .* conj(waveform(n2))));
    end

    %find a few max peaks
    CP_corr_copy = CP_corr;
    npeak = 50;
    CP_pos = zeros(1, npeak);
    for n = 1:npeak
        [~, CP_pos(n)] = max(CP_corr_copy);
        %don't take indexes from beyond array length
        start = max(1, CP_pos(n)-N_sym);
        stop = min(length(CP_corr_copy), CP_pos(n)+N_sym);
        %clean this symbol to ensure that we don't take it again
        CP_corr_copy(start:stop) = zeros(1, stop-start+1);
    end
    
    %calculate CFO for choosen symbols
    CP_pos = sort(CP_pos);
    indv_CFO = zeros(1, npeak);
    for n = 1:npeak 
        n1 = CP_pos(n):CP_pos(n)+N_CP-1;
        n2 = n1 + N_FFT;
        indv_CFO(n) = sample_rate/N_FFT * mean(angle(conj(waveform(n1)) .* waveform(n2)) / (2*pi));
    end
    
    %correct signal with mean of different CFO values
    fractional_CFO = mean(indv_CFO);
end

