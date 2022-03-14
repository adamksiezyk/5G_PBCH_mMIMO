function [pss_indices, NID2, corrPSS] = detectAndDecodePSS(waveform, N_FFT, ...
    N_CP, threshold)
%DETECTANDDECODEPSS Detect and decode PSS sequences in waveform
% Inputs:
%   waveform            : a vector representing the 5G signal
%   N_FFT               : a number representing the FFT size
%   N_CP                : a number representing the CP lenght
%   threshold           : a number representing the threshold for sequence
%   correlation detection, if threshold == 0 then only the best PSS is
%   detected
% Outputs:
%   pss_indices     : a vector representing the detected PSS indexes
%   NID2            : a vector representing the detected PSS NID2
%   corrPSS         : a matrix (n x m) representing the correlation
%   sequence for the n-th NID2

    if isrow(waveform)
        waveform = waveform.';
    end
    N_PSS = 127;            % PSS length
    N_sym = N_FFT + N_CP;
    
    pss_ref = [PSS.generatePSS(0).', PSS.generatePSS(1).', ...
        PSS.generatePSS(2).'];
    pss_start = N_FFT/2 - floor(N_PSS/2);
    m_seqs = zeros(N_FFT, 3);
    m_seqs(pss_start:pss_start+N_PSS-1, :) = pss_ref;

    time_seqs = ifft(ifftshift(m_seqs, 1), N_FFT, 1) * sqrt(N_FFT);
    time_seqs = [time_seqs(end-N_CP+1:end, :); time_seqs];

    corrPSS = [abs(conv(waveform, conj(time_seqs(end:-1:1, 1)))), ...
        abs(conv(waveform, conj(time_seqs(end:-1:1, 2)))), ...
        abs(conv(waveform, conj(time_seqs(end:-1:1, 3))))];

    if threshold == 0
        [max_val, time_offsets] = max(corrPSS);
        [~, max_NID2] = max(max_val);
        pos = time_offsets(max_NID2);
        pss_indices = pos - N_sym + 1;
        NID2 = max_NID2 - 1;
    else
        %find index of every detected pss sequence 
        idx = 1;
        indices = double.empty;
        NID2 = double.empty;
        for pos = 2:length(corrPSS)
            [max_corr, max_NID2] = max(corrPSS(pos, :));
            if ((max_corr > threshold) && (max_corr > corrPSS(pos-1, max_NID2))...
                    && (max_corr > corrPSS(pos+1, max_NID2)))
                indices(idx) = pos;
                NID2(idx) = max_NID2 - 1;
                idx = idx + 1;
            end
        end
        pss_indices = indices - N_sym + 1;
    end
end