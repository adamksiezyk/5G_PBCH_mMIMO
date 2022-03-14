clear variables; close all;
%% Load waveform and initial parameters
fprintf(" -- Loading waveform and initial parameters --\n");

show_plots = true;
signal_choice = 1;
if (signal_choice == 1)
    % Signal deflections from radar measurements
    load('signals/ssb_5ms_periodicity/2021-12-09_12-31-43_GPS_Fc3440.0_G45.0_Bw50.0_Fn000_ch1.mat');
%     waveform = waveform(1:25.0e-3*61.44e6).';    % First 3 frames
%     waveform = waveform.';
    
    signal_info = SignalInfo;
    signal_info.fs = 61.44e6;
    signal_info.fc = 3440e6;
    signal_info.BW = 50e6;
    signal_info.N_FFT = 2048;
    signal_info.SCS = 30e3;
    
    signal_info.SSB.SSB_case = "Case C";
    signal_info.SSB.subcarrier_offset = -507;
    
    PSS_detection_threshold = 0;    % Threshold for SSB detection
    int_CFO = 0;                    % Integer Center Frequency Offset
    
    available_channel_BWs = [10, 15, 20, 25, 30, 40, 50, 60, 70, 80, ...
        90, 100] * 1e6;             % 3GPP 38.101-1 5.3.5-1
elseif (signal_choice == 2)
    % Orange Krków 5G signal
    load('signals/NR_DL_2159.91MHz_10MHz.mat'); %ncellid = 440
    signal_info = SignalInfo;
    signal_info.fs = 15360000;
    signal_info.fc = 2159.91e6;
    signal_info.BW = 10e6;
    signal_info.N_FFT = 1024;
    signal_info.SCS = 15e3;
    
    signal_info.SSB.SSB_case = "Case C";
    signal_info.SSB.subcarrier_offset = -48;
    
    PSS_detection_threshold = 60;   % Threshold for SSB detection
    int_CFO = 0;                    % Integer Center Frequency Offset
elseif (signal_choice == 3)
    % Synthetic signal with beamforming and SIB1
    load('signals/matlab_5g_downlink_5MHzBW_30kHzSCS.mat');
    waveform = tmp_waveform.';
    
    signal_info = SignalInfo;
    signal_info.fs = 15360000;
    signal_info.fc = 3e9;
    signal_info.BW = 10e6;
    signal_info.N_FFT = 512;
    signal_info.SCS = 30e3;
    
    signal_info.SSB.SSB_case = "Case B";
    signal_info.SSB.subcarrier_offset = 0;
    
    PSS_detection_threshold = 0;    % Threshold for SSB detection
    int_CFO = 0;                    % Integer Center Frequency Offset
    
    available_channel_BWs = [5, 10, 40] * 1e6; % 3GPP 38.101-1 5.3.5-1
end

N = length(waveform);
signal_info.SSB.L_SSB = utils.getLSSB(signal_info.fc);

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

%% Estimation and correct the fractional CFO
fprintf(" -- Fractional CFO estimation and correction --\n");
fractional_CFO = processing.estimateFractionalCFO(waveform, signal_info, ...
    show_plots);
fprintf("Fractional CFO = %.2f\n", fractional_CFO);
waveform = waveform .* exp(-1j*2*pi*fractional_CFO/signal_info.fs*(0:N-1));

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

%% Find SSB position in time domain and detect PSS
fprintf(" -- Find SSB position in time domain --\n");
[SSB_indices, SSBs, NID2_list, NID1_list] = processing.detectSSBs(...
    waveform, signal_info, PSS_detection_threshold, show_plots);

%% Process each SSB
for i = 1:length(SSB_indices)
    %%%%%%%%%% Decode PBCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SSB_idx = SSB_indices(i);
    SSB_grid = SSBs(i, :, :);
    NID2 = NID2_list(i);
    NID1 = NID1_list(i);

    % Decode Cell ID
    cell_id = 3 * NID1 + NID2;
    fprintf("CellID = %d\n", cell_id);

    % Decode BCH
    [MIB, ~, iSSB] = processing.decodePBCH(signal_info, SSB_grid, ...
        cell_id, show_plots);
    
    if isempty(MIB)
        % Skip when CRC error detected
        continue;
    end
    
    %%%%%%%%%% Check if CORESET#0 is present %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if all(MIB.kSSB > 23)
        % Another SSB with CORESET#0 information can be found according to 3GPP 38.213 13
        fprintf("No CORESET#0 present\n");
        continue;
    end
    
    
    %%%%%%%%%% Demodulate OFDM resource grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(" -- Demodulate OFDM resource grid --\n");
    
    % Estimate and correct timing offset to the frame start
    timing_offset = utils.estimateTimingOffset(signal_info.fs, ...
        signal_info.N_sym_long, signal_info.N_sym, ...
        signal_info.SSB.SSB_case, signal_info.fc, SSB_idx-1, iSSB);
    waveform_frame = waveform(1+timing_offset:end);
    
    % Estimate and correct frequency offset to the common resource block
    SCS_common = MIB.SubcarrierSpacingCommon * 1e3;
    f_shift = MIB.kSSB * SCS_common;
    waveform_frame = waveform_frame .* ...
        exp(1j*2*pi*f_shift/signal_info.fs*(0:length(waveform_frame)-1));
    
    msb_idx = floor(MIB.PDCCHConfigSIB1/16);   % 4 MSB of PDCCHConfigSIB1
    lsb_idx = mod(MIB.PDCCHConfigSIB1, 16);    % 4 LSB of PDCCHConfigSIB1
    SCS_pair = [signal_info.SCS, SCS_common];
    min_channel_BW = min(available_channel_BWs);
    
    % Demodulate the resource grid
    [N_RB_CORESET, N_sym_CORESET, CORESET_RB_offset, pattern] = ...
        PDCCH.getCORESET0Resources(msb_idx, SCS_pair, min_channel_BW, ...
        MIB.kSSB);
    N_slots_in_signal = floor(length(waveform_frame) / ...
        (signal_info.N_sym_long+13*signal_info.N_sym));
    N_slots_per_frame = signal_info.N_subframes_per_frame * ...
        utils.getSlotsPerSubframe(SCS_common);
    repeat = min(N_slots_in_signal, N_slots_per_frame);
    N_CPs = repmat([signal_info.N_CP_long, ...
        signal_info.N_CP*ones(1, 13)], 1, repeat);          % Minimum number of symbols to cover CORESET#0
    c0 = CORESET_RB_offset + 10*signal_info.SCS/SCS_common; % CORESET RB offset from carrier center (+padding)
    N_RB_CORESET_min = 2*max(c0, N_RB_CORESET-c0);          % Minimum number of RBs to cover CORESET#0
    N_subcarriers = N_RB_CORESET_min * signal_info.N_subcarriers_per_RB;    
    subcarriers = (-N_subcarriers/2+1:N_subcarriers/2)+signal_info.N_FFT/2;
%     resource_grid = utils.demodulateOFDM(waveform_frame, N_CPs, ...
%         signal_info.N_FFT);
%     resource_grid = resource_grid(subcarriers, :);
    resource_grid = nrOFDMDemodulate(waveform_frame.', N_RB_CORESET_min, ...
        30, 0, 'SampleRate', signal_info.fs, 'CarrierFrequency', 0);
    
    if show_plots
        figure;
        imagesc(abs(resource_grid));
        axis xy;
        title("Frame resource grid");
        xlabel("symbols");
        ylabel("subcarriers");
        fprintf("Press ENTER to continue ...\n");
        pause;
    end
    

    %%%%%%%%%% Decode PDCCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(" -- Find PDCCH --\n");
    pdcch = PDCCH.getPDCCHParameters(signal_info.SCS, min_channel_BW, ...
        MIB, iSSB, cell_id);
    
    % Create a carrier that spans the wole CORESET#0 BWP
    c0Carrier = Carrier;
    c0Carrier.SubcarrierSpacing = MIB.SubcarrierSpacingCommon;
    c0Carrier.N_RB_start = pdcch.CORESET.N_RB_start;
    c0Carrier.N_RB = pdcch.CORESET.N_RB;
    c0Carrier.NSlot = pdcch.SearchSpace.SlotOffset;
    c0Carrier.NFrame = MIB.NFrame;
    c0Carrier.NCellID = cell_id;
    
    % Find Type0-PDCCH monitoring occasions. 3GPP 28.213
    N_symbols = size(resource_grid, 2);
    [mon_slots, mon_symbols, mon_subcarriers] = ...
        PDCCH.getMonitoringOccasionsResources(lsb_idx, iSSB, SCS_pair, ...
        MIB.NFrame, pattern, N_sym_CORESET, N_RB_CORESET, N_RB_CORESET_min, ...
        CORESET_RB_offset, N_symbols, signal_info, c0Carrier);
    mon_grid = resource_grid(mon_subcarriers, ...
        mon_symbols);
    
    % Get PDCCH resources
    [pdcch_indices, pdcch_dmrs_indices]  = PDCCH.getPDCCHResources(...
        c0Carrier, pdcch);
    
    % Decode PDCCH
    [dci, dci_crc] = processing.decodePDCCH(mon_grid, c0Carrier, pdcch, ...
    pdcch_dmrs_indices, pdcch_indices, mon_slots, show_plots);
    if dci_crc ~= 0
        fprintf("Could not decode DCI\n");
        continue;
    end
    
    % Get PDSCH parameters
    [pdsch, K_0] = PDSCH.getPDSCHResources(dci, pdcch.CORESET.N_RB, ...
        MIB.DMRSTypeAPosition, pattern);
end