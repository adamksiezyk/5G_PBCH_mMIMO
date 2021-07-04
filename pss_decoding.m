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
function [indxs, NID2] = pss_decoding(waveform, nfft, threshold, cfo_test)

    %generate sequence for every nid2 and map them 
    pss_start = nfft/2-63 + cfo_test;
    m_seqs = zeros(nfft, 3);
    for NID = 0:2
        m_seqs(pss_start:pss_start+126, NID+1) = pss_generate(NID).';
    end

    %Transform sequences to time domain
    time_seqs = zeros(nfft, 3);
    for NID = 0:2
        time_seqs(:, NID+1) = ifft(ifftshift(m_seqs(:, NID+1))) * sqrt(nfft);
    end

    %signal cross-corelation with m-sequences
    for k = 1:3
         ccPSS(:, k) = abs(conv(waveform.', conj(time_seqs(end:-1:1, k))));
    end
    figure; 
    subplot(311); plot(ccPSS(:, 1));
    title('Cross-correlation with sequence NID2 = 0'); xlabel('sample index n');
    subplot(312); plot(ccPSS(:, 2));
    title('Cross-correlation with sequence NID2 = 1'); xlabel('sample index n');
    subplot(313); plot(ccPSS(:, 3));
    title('Cross-correlation with sequence NID2 = 2'); xlabel('sample index n');

    %find correct NID
    [~, NID2] = max(max(ccPSS));
    %for CFO estimation this could plot itself many times so be aware
    %figure; plot(ccPSS(:, NID2)); title('Detected PSS sequence');
    ccPSS = ccPSS(:, NID2);
    NID2 = NID2-1;

    %find index of every detected pss sequence 
    idx = 1;
    indxs = double.empty;
    for pos = 2:length(ccPSS)
        if ((ccPSS(pos) > threshold) && (ccPSS(pos) > ccPSS(pos-1)) && (ccPSS(pos) > ccPSS(pos+1)))
            indxs(idx) = pos;
            idx = idx + 1;
        end
    end
end