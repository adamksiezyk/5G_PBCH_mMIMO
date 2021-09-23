function waveform = modulateOFDM(fft_grid, signal_info, N_CPs)
%MODULATEOFDM Summary of this function goes here
%   Detailed explanation goes here

    N_syms = size(fft_grid, 2);     % Number of symbols
    N_FFT = signal_info.N_FFT;  % FFT size
    
    % Construct time sequences
    time_seqs = ifft(ifftshift(fft_grid, 1), N_FFT, 1);
    
    % Add CP and serialize
    waveform = zeros(1, N_syms*N_FFT+sum(N_CPs));
    pos = 1;
    for i = 1:N_syms
        N_CP = N_CPs(i);
        waveform(pos+(0:N_CP+N_FFT-1)) = ...
            [time_seqs(end-(N_CP-1):end, i); time_seqs(:, i)].';
        pos = pos + N_CP + N_FFT;
    end
end

