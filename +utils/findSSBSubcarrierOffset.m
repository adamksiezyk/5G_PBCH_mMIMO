function search_results = findSSBSubcarrierOffset(waveform, nfft, ncp, raster)
%FINDSSB Searches through all potential SSB frequency positions
% Inputs:
%   waveform    : a vector representing the signal
%   nfft        : a number representing the FFT size
%   ncp         : a number representing the CP size
%   raster      : a vector representing the subcarrier offset where the
%   seatch will be performed
%   
% Outputs:
%   search_results : a Nx2 matrix, where N = length(raster) representing
%   the results of correlating a reference PSS signal shifted by raster(n)
%   with the given waveform and NID2 of the reference PSS signal

    search_results = zeros(length(raster), 2);
    status = '';
    fprintf("Searching for SSB frequency position:\t");
    
    for n = 1:length(raster)
        fprintf(repmat('\b', 1, strlength(status)));
        status = sprintf("%d/%d", n, length(raster));
        fprintf(status);
        
        waveform_shift = waveform .* exp(-1j*2*pi*raster(n)/nfft*(1:length(waveform)));
        [~, ~, corrPSS] = PSS.detectAndDecodePSS(waveform_shift, nfft, ncp, 0);

        % find correct NID
        [search_results(n, 1), search_results(n, 2)] = max(max(corrPSS));
    end
    fprintf("\n");
end

