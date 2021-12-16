clear variables; close all;
%% Load waveform and initial parameters
fprintf(" -- Loading waveform and initial parameters --\n");

show_plots = false;

% Signal parameters
load('signals/ssb_5ms_periodicity/2021-12-09_12-31-43_GPS_Fc3440.0_G45.0_Bw50.0_Fn000_ch1.mat');
waveform = waveform(1:5.0e-3*61.44e6);   % First 2.5 frames
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

PSS_detection_threshold = 0;    % Threshold for SSB detection
int_CFO = 0;                    % Integer Center Frequency Offset

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

%% Process each SSB
SSB_idx = SSB_indices(1);
SSB_grid = SSBs(1, :, :);
NID2 = NID2_list(1);
NID1 = NID1_list(1);

% Decode Cell ID
cell_id = 3 * NID1 + NID2;
fprintf("CellID = %d\n", cell_id);

% Decode BCH
[MIB, tr_block, iSSB, HFR] = processing.decodePBCH(signal_info, SSB_grid, ...
    cell_id, show_plots);

% Reconstruct the encoded BCH transport block
cd_block = nrBCH(tr_block, MIB.NFrame, HFR, signal_info.SSB.L_SSB, ...
    MIB.kSSB, cell_id);

% Synthesise SSB
synthetic_SSB = processing.synthesiseSSB(signal_info, cd_block, ...
    NID2, NID1, iSSB, show_plots);

% Delete PBCH DM-RS
% SSB_start = signal_info.N_FFT/2 - signal_info.SSB.N_subcarriers_SSB/2 + ...
%             signal_info.SSB.subcarrier_offset;
% PBCH_DMRS_pos = signal_info.SSB.PBCH_DMRS_position(cell_id);
% PBCH_DMRS_pos = PBCH_DMRS_pos + ...
%     floor(PBCH_DMRS_pos / signal_info.SSB.N_subcarriers_SSB) * ...
%     (signal_info.N_FFT-signal_info.SSB.N_subcarriers_SSB);
% PBCH_DMRS_pos = PBCH_DMRS_pos + SSB_start - 1;
% synthetic_SSB(PBCH_DMRS_pos) = zeros(1, length(PBCH_DMRS_pos));

% Synthesize waveform
N_CPs = ones(1, signal_info.SSB.N_symbols_SSB) * signal_info.N_CP;
synthetic_waveform = utils.modulateOFDM(synthetic_SSB, signal_info, N_CPs);

if show_plots
    fprintf(" -- Synthesised signal Spectrigram --\n");
    figure;
    f = signal_info.fs/2048 * (-1024:1023);
    spectrogram(synthetic_waveform, 2048, 1024, f, signal_info.fs);
    fprintf("Press ENTER to continue ...\n");
    pause;
end