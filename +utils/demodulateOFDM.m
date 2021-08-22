function grid = demodulateOFDM(waveform, N_CPs, N_FFT)
%DEMODULATEOFDM Demodulates the OFDM modulated waveform into a resource
%grid. The number of deomodulated symbols is determined by the N_CPs vector
%length
% Inputs:
%   waveform        : a vector representing the waveform
%   N_CPs           : a vector representing the OFDM symbols cyclic prefix
%   lengths 
%   N_FFT           : a number representing the modulation FFT size
% Outputs:
%   grid    : a matrix representing the OFDM demodulated resource grid

    N_OFDM_symbols = length(N_CPs);                 % Number of OFDM symbols
    grid = zeros(N_FFT, N_OFDM_symbols);            % Resource grid
    pos = 1;                                        % Current OFDM symbol position
    for n = 1:N_OFDM_symbols
        N_sym = N_CPs(n) + N_FFT;                               % Symbol length
        samples = waveform(pos+(N_CPs(n):+N_sym-1));            % Remove CP
        grid(:, n) = fftshift(fft(samples)) / sqrt(N_FFT);      % Demodulate OFDM
        pos = pos + N_sym;                                      % Move to next symbol
    end
end