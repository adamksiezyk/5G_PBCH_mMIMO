clear variables; close all;
%% Load waveform and initial parameters
fprintf(" -- Loading waveform and initial parameters --\n");

show_plots = false;

% Signal parameters
load('signals/2021-05-21_11-05-13_GPS_Fc3440.0_G50.0_Bw50.0_Fn000_ch2_spline_61.44MHz.mat');
fs = 61.44e6;
waveform = waveform(1:25.0e-3*fs).';    % First 2.5 frames
fc = 3440e6;
BW = 50e6;
N_FFT = 2048;
SCS_SSB = 30e3;
threshold = 2e5;                                % Threshold for SSB detection
SSB_case = "Case C";
available_channel_BWs = [10, 15, 20, 25, 30, 40, 50, 60, 70, 80, ...
    90, 100] * 1e6;                             % 3GPP 38.101-1 5.3.5-1
SSB_subcarrier_offset = -507;                   % SSB subcarrier offset
int_CFO = 0;                                    % Integer Center Frequency Offset

[SSB_start, synthetic_SSBs] = processing.synthesis(waveform, SCS_SSB, ...
    BW, fc, fs, N_FFT, threshold, SSB_subcarrier_offset, show_plots);