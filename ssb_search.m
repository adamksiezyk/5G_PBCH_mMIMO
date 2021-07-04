function search_results = ssb_search(waveform, nfft, raster)
%SSB_SEARCH Searches through all potential SSB frequency positions
% Inputs:
%   waveform    : a vector representing the signal
%   nfft        : a number representing the FFT size
%   raster      : a vector representing the subcarrier offset where the
%   seatch will be performed
%   
% Outputs:
%   search_results : a Nx2 matrix, where N = length(raster) representing
%   the results of correlating a reference PSS signal shifted by raster(n)
%   with the given waveform and NID2 of the reference PSS signal

    if (isrow(waveform))
        waveform = waveform.';
    end
    N_PSS = 127;
    pss_ref = [pss_generate(0).', pss_generate(1).', pss_generate(2).'];
    search_results = zeros(length(raster), 2);
    status = '';
    fprintf("Searching for SSB frequency position:\t");
    
    for n = 1:length(raster)
        fprintf(repmat('\b', 1, strlength(status)));
        status = sprintf("%d/%d", n, length(raster));
        fprintf(status);
        
        % Generate sequence for every NID2 and map them 
        pss_start = nfft/2 - floor(N_PSS/2) + raster(n);
        m_seqs = zeros(nfft, 3);
        m_seqs(pss_start:pss_start+N_PSS-1, :) = pss_ref;

        % Transform sequences to time domain
        time_seqs = ifft(ifftshift(m_seqs, 1), nfft, 1) * sqrt(nfft);

        % Signal cross-corelation with m-sequences
        corrPSS = [abs(conv(waveform, conj(time_seqs(end:-1:1, 1)))), ...
            abs(conv(waveform, conj(time_seqs(end:-1:1, 2)))), ...
            abs(conv(waveform, conj(time_seqs(end:-1:1, 3))))];

        % find correct NID
        [search_results(n, 1), search_results(n, 2)] = max(max(corrPSS));
    end
    fprintf("\n");
end

