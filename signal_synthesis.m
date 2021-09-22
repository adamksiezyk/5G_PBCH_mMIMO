clear variables; close all;
%% Load waveform and initial parameters
fprintf(" -- Loading waveform and initial parameters --\n");

show_plots = true;

% Signal parameters
load('signals/2021-05-21_11-05-13_GPS_Fc3440.0_G50.0_Bw50.0_Fn000_ch2_spline_61.44MHz.mat');
waveform = waveform(1:25.0e-3*61.44e6).';   % First 2.5 frames
N = length(waveform);
signal_info = SignalInfo;
signal_info.fs = 61.44e6;
signal_info.fc = 3440e6;
signal_info.BW = 50e6;
signal_info.N_FFT = 2048;
signal_info.SCS = 30e3;

signal_info.SSB.SSB_case = "Case C";
signal_info.SSB.subcarrier_offset = -507;
signal_info.SSB.L_SSB = utils.getLSSB(signal_info.fc);

PSS_detection_threshold = 2e5;      % Threshold for SSB detection
int_CFO = 0;                        % Integer Center Frequency Offset

%% Plot the signal PSD and Spectrogram
    if show_plots
        fprintf(" -- Signal PSD and Spectrigram --\n");
        figure;
        pwelch(waveform, 2048, 2048-1024, 2048, signal_info.fs, 'centered');

        figure;
        f = signal_info.fs/2048 * (-1024:1023);
        spectrogram(waveform, 2048, 1024, f, signal_info.fs);
        fprintf("Press ENTER to continue ...\n");
        pause;
    end

%% Find SSB in frequency domain
if isempty(signal_info.SSB.subcarrier_offset)
    fprintf(" -- Find SSB in frequency domian --\n");
    signal_info.SSB.subcarrier_offset = ...
        processing.findSSBSubcarrierOffset(waveform, signal_info, ...
        show_plots);
end

%% Estimation and correct the integer CFO
if (exist('int_CFO', 'var') ~= 1)
    fprintf(" -- Integer CFO estimation and correction --\n");
    int_CFO = processing.estimateIntegerCFO(waveform, signal_info, ...
        show_plots);
end
fprintf("Integer CFO: %d\n", int_CFO);
waveform = waveform .* exp(-1j*2*pi*int_CFO/signal_info.N_FFT*(0:N-1));

%% Estimation and correct the fractional CFO
fprintf(" -- Fractional CFO estimation and correction --\n");
fractional_CFO = processing.estimateFractionalCFO(waveform, signal_info, ...
    show_plots);
fprintf("Fractional CFO = %.2f\n", fractional_CFO);
waveform = waveform .* exp(-1j*2*pi*fractional_CFO/signal_info.fs*(0:N-1));

%% Find SSB position in time domain and detect PSS
fprintf(" -- Find SSB position in time domain --\n");
[SSB_indices, SSBs, NID2_list, NID1_list] = processing.detectSSBs(...
    waveform, signal_info, PSS_detection_threshold, show_plots);
% PCI calculation
SSB_amount = length(SSB_indices);
cell_ids = zeros(1, SSB_amount);
for i = 1:SSB_amount
    cell_ids(i) = 3 * NID1_list(i) + NID2_list(i);
    fprintf("CellID = %d\n", cell_ids(i));
end

%% Decode BCH
[MIBs, PBCH_bits, iSSBs] = processing.decodePBCH(signal_info, SSBs, ...
    cell_ids, show_plots);

%% Synthesise SSBs
synthetic_SSBs = processing.synthesiseSSBs(signal_info, PBCH_bits, NID2_list, NID1_list, iSSBs, show_plots);