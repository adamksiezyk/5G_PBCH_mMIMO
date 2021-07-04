%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pss_decoding(waveform, nfft, sample_rate)
% Input:
%   waveform - 5g downlink data vector
%   nfft - data fft size
%   nb_rb - number of allocated resource blocks
%   threshold - threshold for sequence correlation detection
% Output:
%   indxs - positions of every detected pss sequence
%   NID2 - signal correct NID2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indxes, NID2] = pss_detect(waveform, nfft, threshold, offset)
%PSS_DETECT Detect PSS sequences in waveform
% Inputs:
%   waveform    : a vector representing the 5G signal
%   nfft        : a number representing the FFT size
%   threshold   : a number representing the threshold for sequence
%   correlation detection
%   offset      : a number representing the subcarrier offset of the PSS
%   sequence in the 5G signal
% Outputs:
%   indexes     : a vector representing the detected PSS indexes
%   NID2        : a number representing the NID2 of the detected PSS

    if (isrow(waveform))
        waveform = waveform.';
    end
    N_PSS = 127;    % PSS length
    
    pss_start = nfft/2 - floor(N_PSS/2) + offset;
    m_seqs = zeros(nfft, 3);
    m_seqs(pss_start:pss_start+N_PSS-1, :) = pss_ref;

    time_seqs = ifft(ifftshift(m_seqs, 1), nfft, 1) * sqrt(nfft);

    corrPSS = [abs(conv(waveform, conj(time_seqs(end:-1:1, 1)))), ...
        abs(conv(waveform, conj(time_seqs(end:-1:1, 2)))), ...
        abs(conv(waveform, conj(time_seqs(end:-1:1, 3))))];
    
    figure; 
    subplot(311); plot(corrPSS(:, 1));
    title('Cross-correlation with sequence NID2 = 0'); xlabel('sample index n');
    subplot(312); plot(corrPSS(:, 2));
    title('Cross-correlation with sequence NID2 = 1'); xlabel('sample index n');
    subplot(313); plot(corrPSS(:, 3));
    title('Cross-correlation with sequence NID2 = 2'); xlabel('sample index n');

    [~, NID2] = max(max(corrPSS));
    corrPSS = corrPSS(:, NID2);
    NID2 = NID2-1;

    %find index of every detected pss sequence 
    idx = 1;
    indxes = double.empty;
    for pos = 2:length(corrPSS)
        if ((corrPSS(pos) > threshold) && (corrPSS(pos) > corrPSS(pos-1)) && (corrPSS(pos) > corrPSS(pos+1)))
            indxes(idx) = pos;
            idx = idx + 1;
        end
    end
end