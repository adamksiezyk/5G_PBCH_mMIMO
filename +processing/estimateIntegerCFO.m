function int_CFO = estimateIntegerCFO(waveform, signal_info, show_plots_)
%ESTIMATEINTEGERCFO Returns the estimated Center Frequency Offset
% Inputs:
% Inputs:
%   waveform        : a column vector that represents the waveform
%   signal_info     : a SignalInfo object
%   show_plots_     : a boolean, if true plots are shown
% Outputs:
%   int_CFO     : a number that represents the estimated integer CFO

    if nargin < 3
        show_plots = false;
    else
        show_plots = show_plots_;
    end

    offsets = -10:10;
    N = min([length(waveform), 25e-3*signal_info.fs]);
    search_results = utils.findSSBSubcarrierOffset(waveform(1:N), ...
        signal_info.N_FFT, signal_info.SSB.subcarrier_offset+offsets);

    [max_val, max_ind] = max(search_results(:, 1));
    int_CFO = offsets(max_ind);
    
    if show_plots
        figure;
        hold on;
        plot(offsets, search_results(:, 1));
        plot(int_CFO, max_val, 'kx', 'LineWidth', 2, 'MarkerSize', 8);
        title("SSB frequency offset search results");
        ylabel("Maximum of xcorr");
        xlabel("Subcarrier offset");
        legend("correlation result", "intefer CFO");
        fprintf("Press ENTER to continue ...\n");
        pause;
    end
end

