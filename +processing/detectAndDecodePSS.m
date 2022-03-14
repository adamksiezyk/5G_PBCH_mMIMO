function [pss_indices, NID2] = detectAndDecodePSS(waveform, N_FFT, ...
    N_CP, threshold, subcarrier_offset, show_plots_)
%DETECTANDDECODEPSS Detect and decode PSS sequences in waveform
% Inputs:
%   waveform            : a vector representing the 5G signal
%   N_FFT               : a number representing the FFT size
%   N_CP                : a number representing the CP lenght
%   threshold           : a number representing the threshold for sequence
%   correlation detection, if threshold == 0 then only the best PSS index
%   and it's NID2 are returned
%   subcarrier_offset   : a number representing the subcarrier offset of
%   the PSS sequence in the 5G signal
%   show_plots_         : a boolean representing the show plots flag
%   default is true
% Outputs:
%   pss_indices     : a vector (or a single number if threshold == 0)
%   representing the detected PSS indexes
%   NID2            : a vector (or a single number if threshold == 0)
%   representing the detected PSS NID2

    if nargin<6
        show_plots = true;
    else
        show_plots = show_plots_;
    end

    waveform_shift = waveform .* exp(-1j*2*pi*subcarrier_offset/N_FFT*(1:length(waveform)));
    [pss_indices, NID2, corrPSS] = PSS.detectAndDecodePSS(waveform_shift, ...
        N_FFT, N_CP, threshold);

    if show_plots
        figure;
        subplot(311);
        hold on;
        plot(corrPSS(:, 1));
        title('Cross-correlation with sequence N_{ID}^{(2)} = 0');
        xlabel('sample index n');
        subplot(312);
        hold on;
        plot(corrPSS(:, 2));
        title('Cross-correlation with sequence N_{ID}^{(2)} = 1');
        xlabel('sample index n');
        subplot(313);
        hold on;
        plot(corrPSS(:, 3));
        title('Cross-correlation with sequence N_{ID}^{(2)} = 2');
        xlabel('sample index n');
    end

    if show_plots
        for i = 1:length(NID2)
            subplot(3, 1, NID2(i));
            plot(pss_indices(i), corrPSS(i), 'kx', 'LineWidth', 2, 'MarkerSize', 8);
            legend("correlation result", "N_{ID}^{(2)}");
        end
    end
end