function fractional_CFO = estimateFractionalCFO(waveform, signal_info, ...
    show_plots_)
%ESTIMATEFRACTIONALCFO Returns the estimated Fractional Center Frequency
%Offset
% Inputs:
%   waveform        : a column vector that represents the waveform
%   signal_info     : a SignalInfo object
%   show_plots_     : a boolean, if true plots are shown
% Outputs:
%   fractional_CFO  : a number that represents the estimated fractional CFO

    if nargin < 3
        show_plots = false;
    else
        show_plots = show_plots_;
    end

    [fractional_CFO, CP_corr, indv_CFO] = utils.estimateFractionalCFO(waveform, ...
        signal_info.fs, signal_info.SCS);

    if show_plots
        figure;
        plot(abs(CP_corr));
        title('Signal autocorrelation using CP of SSB');
        xlabel('Sample index n');
        
        figure;
        plot(indv_CFO, 'o-');
        title('Estimated fractional CFO');
        xlabel('sample index n');
        ylabel('[Hz]');
        
        fprintf("Press ENTER to continue ...\n");
        pause;
    end
end

