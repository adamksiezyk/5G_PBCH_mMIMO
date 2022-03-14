function SSB_subcarrier_offset = findSSBSubcarrierOffset(waveform, ...
    signal_info, show_plots_)
%FINDSSBPOSITION SSBs can be placed in a fixed frequency raster called the
%GSCN. 3GPP 38.104 Table 5.4.3.1-1
% Inputs:
%   waveform        : a column vector that represents the waveform
%   signal_info     : a SignalInfo object
%   show_plots_     : a boolean, if true plots are shown
% Outputs:
%   SSB_subcarrier offset   : a number representing the SSB subcarrier
%   offset

    if nargin < 3
        show_plots = false;
    else
        show_plots = show_plots_;
    end
    
    raster = utils.getSSBRaster(signal_info.fc, signal_info.BW, ...
        signal_info.SCS);
    N = min([length(waveform), 25e-3*signal_info.fs]);
    search_results = utils.findSSBSubcarrierOffset(waveform(1:N), ...
        signal_info.N_FFT, signal_info.N_CP, raster);
    [max_val, max_ind] = max(search_results(:, 1));
    SSB_subcarrier_offset = raster(max_ind);

    if show_plots
        figure;
        hold on;
        stem(raster, search_results(:, 1));
        plot(SSB_subcarrier_offset, max_val, 'kx', 'LineWidth', 2, ...
            'MarkerSize', 8);
        title("SSB frequency position search results");
        ylabel("Maximum of xcorr");
        xlabel("Subcarrier offset");
        legend("correlation result", "SSB position");
        fprintf("Press ENTER to continue ...\n");
        pause;
    end
end

