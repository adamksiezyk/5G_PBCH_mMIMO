function [fractional_CFO, CP_corr] = estimateFractionalCFO(waveform, ...
    sample_rate, scs, ssb_frequency_offset)
%ESTIMATEFRACTIONALCFO Estimated the fractional center frequency offset
% Inputs:
%   waveform                : a column vector representing the signal
%   sample_rate             : a number representing the sample rate in Hz
%   scs                     : a number representing the subcarrier spacing
%   in Hz
%   ssb_frequency_offset    : a number representing the SSB frequency
%   offset in relation to the center frequecy in Hz
% Outputs:
%   fractional_CFO  : a number representing the estimated fractional center
%   frequency offset
%   CP_corr         : a vector representing the calculated cyclic prefix
%   based correlation

    % Shift the SSB to center frequency
    t = (0:length(waveform)-1) / sample_rate;
    waveform_center = waveform .* exp(-1i*2*pi*ssb_frequency_offset*t);
    
    % Downsample
    N_FFT_ds = 256;
    sample_rate_ds = N_FFT_ds * scs;
    [~, N_CP_ds] = utils.getCPLength(scs, sample_rate_ds);
    N_sym_ds = N_CP_ds + N_FFT_ds;
    waveform_ds = resample(waveform_center, sample_rate_ds, sample_rate);
    
    % Find PSS time indices
    pss_indices = PSS.detectAndDecodePSS(waveform_ds, N_FFT_ds, ...
        N_CP_ds, 0, 0, false);
    time_offset = pss_indices(1);
    waveform_shift = waveform_ds(time_offset:end);

    % Multiply the waveform by itself delayed by N_sym samples and conjugated
    delayed = [zeros(1, N_FFT_ds), waveform_shift(1:end-N_FFT_ds)];
    cpProduct = waveform_shift .* conj(delayed);

    % Apply a moving sum filter with a window size equal to the CP length
    CP_corr = filter(ones([N_CP_ds 1]), 1, cpProduct);

    % Moving sum over 4 OFDM symbols (SSB size)
    CP_corr_delayed = CP_corr;
    for k = 1:3
        CP_corr_delayed = [zeros(1, N_sym_ds), ...
            CP_corr_delayed(1:end-N_sym_ds)];
        CP_corr = CP_corr + CP_corr_delayed;
    end
    
    % Estimate frequency offset
    CP_corr_idx = N_FFT_ds + N_CP_ds + 3*N_sym_ds; % 4*N_sym = N_FFT (delayed correlation) + N_CP (sum filter CP) + 3*N_sym (moving sum filter of 4 symbols)
    fractional_CFO =  scs * angle(mean(CP_corr(CP_corr_idx))) / (2*pi);
end

