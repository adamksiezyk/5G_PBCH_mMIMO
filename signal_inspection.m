clear variables; close all;
%% Load waveform and initial parameters
fprintf(" -- Loading waveform and initial parameters --\n");

show_plots = true;
signal_choice = 5;
if (signal_choice == 1)
    % Signal deflections from radar measurements
    load('signals/2021-05-21_11-05-13_GPS_Fc3440.0_G50.0_Bw50.0_Fn000_ch2_spline_61.44MHz.mat');
    sample_rate = 61.44e6;
    waveform = waveform(1:25.0e-3*sample_rate).';    % First 2.5 frames
%     waveform = waveform.';
    fc = 3440e6;
    BW = 50e6;
    N_FFT = 2048;
    SCS_SSB = 30e3;
    threshold = 2e5;
    SSB_case = "Case C";
    available_channel_BWs = [10, 15, 20, 25, 30, 40, 50, 60, 70, 80, ...
        90, 100] * 1e6; % 3GPP 38.101-1 5.3.5-1
    SSB_subcarrier_offset = -507; % SSB frequency position
    int_CFO = 0;    % Integer Center Frequency Offset
elseif (signal_choice == 2)
    % Line of sight signal from radar measurements with deleted uplink
    load('signals/Fn001_chan1.mat');
    sample_rate = 61.44e6;
    waveform = chan1;
    fc = 3440e6;
    BW = 50e6;
    N_FFT = 2048;
    SCS_SSB = 30e3;
    threshold = 2e5;
    SSB_case = "Case C";
    available_channel_BWs = [10, 15, 20, 25, 30, 40, 50, 60, 70, 80, ...
        90, 100] * 1e6; % 3GPP 38.101-1 5.3.5-1
    SSB_subcarrier_offset = -507; % SSB frequency position
    int_CFO = 0;    % Integer Center Frequency Offset
elseif (signal_choice == 3)
    % Orange Krków 5G signal
    load('signals/NR_DL_2159.91MHz_10MHz.mat'); %ncellid = 440
    fc = 2159.91e6;
    BW = 10e6;
    N_FFT = 1024;
    sample_rate = 15360000;
    SCS_SSB = 15e3;
    threshold = 60; %
    SSB_subcarrier_offset = -48;  % SSB frequency position
    int_CFO = 0;    % Integer Center Frequency offset
elseif (signal_choice == 4)
    % Radar measurements with deleted uplink
    load('signals/dron_monostatycznie_beamforming/signals_Fn000-0025_chan2/2021-05-21_15-45-36_GPS_Fc3440.0_G50.0_Bw50.0_Fn025_ch2.mat');
    sample_rate = 61.44e6;
    % waveform = waveform(1:40e-3*sample_rate);
    fc = 3440e6;
    BW = 50e6;
    N_FFT = 2048;
    SCS_SSB = 30e3;
    threshold = 0;
    SSB_case = "Case C";
    available_channel_BWs = [10, 15, 20, 25, 30, 40, 50, 60, 70, 80, ...
        90, 100] * 1e6;             % 3GPP 38.101-1 5.3.5-1
    SSB_subcarrier_offset = -507;   % SSB frequency position
    int_CFO = 0;                    % Integer Center Frequency Offset
elseif (signal_choice == 5)
    % Synthetic signal with beamforming and SIB1
    load('signals/matlab_5g_downlink_5MHzBW_30kHzSCS.mat');
    sample_rate = 15360000;
    waveform = tmp_waveform.';
    fc = 3e9;
    BW = 10e6;
    N_FFT = 512;
    SCS_SSB = 30e3;
    threshold = 8;
    SSB_case = "Case B";
    available_channel_BWs = [5, 10, 40] * 1e6; % 3GPP 38.101-1 5.3.5-1
    SSB_subcarrier_offset = 0;
    int_CFO = 0;
end

N_subframes_per_frame = 10;
N_symbols_per_slot = 14;
N_subcarriers_per_RB = 12;
N_subcarriers_SSB = 240;
N_subcarriers_PSS = 127;
N = length(waveform);
[N_CP_long, N_CP] = utils.getCPLength(SCS_SSB, sample_rate);
N_sym_long = N_FFT + N_CP_long;
N_sym = N_FFT + N_CP;
N_RBs = utils.getRBAmount(SCS_SSB, BW);
L_SSB = utils.getLSSB(fc);

%% Plot the signal PSD and Spectrogram
if show_plots
    fprintf(" -- Signal PSD and Spectrigram --\n");
    figure;
    pwelch(waveform, 2048, 2048-1024, 2048, sample_rate, 'centered');

    figure;
    f = sample_rate/2048 * (-1024:1023);
    spectrogram(waveform, 2048, 1024, f, sample_rate);
    fprintf("Press ENTER to continue ...\n");
    pause;
end

%% Find SSB in frequency domain
if (exist('SSB_subcarrier_offset', 'var') ~= 1)
    fprintf(" -- Find SSB in frequency domian --\n");
    raster = utils.getSSBRaster(fc, BW, SCS_SSB);
    search_results = utils.findSSB(waveform, N_FFT, raster);
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

%% Estimation and correct the integer CFO
if (exist('int_CFO', 'var') ~= 1)
    fprintf(" -- Integer CFO estimation and correction --\n");
    offsets = -10:10;
    search_results = utils.findSSB(waveform, N_FFT, ...
        SSB_subcarrier_offset+offsets);

    [max_val, max_ind] = max(search_results(:, 1));
    int_CFO = offsets(max_ind);
    fprintf("Integer CFO: %d\n", int_CFO);
    waveform = waveform .* exp(-1j*2*pi*int_CFO/N_FFT*(0:N-1));

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

%% Estimation and correct the fractional CFO
fprintf(" -- Fractional CFO estimation and correction --\n");
ssb_frequency_offset = SSB_subcarrier_offset * SCS_SSB;
[factional_CFO, CP_corr] = utils.estimateFractionalCFO(waveform, ...
    sample_rate, SCS_SSB, ssb_frequency_offset);
fprintf("Fractional CFO = %.2f\n", factional_CFO);
waveform = waveform .* exp(-1j*2*pi * factional_CFO/sample_rate *(0:N-1));

if show_plots
    figure;
    plot(abs(CP_corr));
    title('Signal autocorrelation using CP of SSB');
    xlabel('Sample index n');
    fprintf("Press ENTER to continue ...\n");
    pause;
end

%% Find SSB position in time domain and detect PSS
fprintf(" -- Find SSB position in time domain --\n");
[pss_indices, NID2] = PSS.detectAndDecodePSS(waveform, N_FFT, N_CP, ...
    threshold, SSB_subcarrier_offset, show_plots);
if show_plots
    fprintf("Press ENTER to continue ...\n");
    pause;
end

%% Process each SSB
for i = 1:length(pss_indices)
    %%%%%%%%%%% Deocode Cell ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(" -- Decoding Cell ID --\n");
    pss_idx = pss_indices(i);
    NID2_ = NID2(i);
    fprintf("NID2 = %d\n", NID2_);
    
    % Second CFO estimation using PSS
    n1 = pss_idx:pss_idx+N_CP-1;
    n2 = n1 + N_FFT;
    pss_CFO = sample_rate/N_FFT * mean(angle(conj(waveform(n1)) .* ...
        waveform(n2)) / (2*pi));

    % Current SSB correction
    block_size = (pss_idx:pss_idx+4*N_sym-1);
    waveform(block_size) = waveform(block_size) .* exp(-1j*2*pi * ...
        pss_CFO/sample_rate *(0:length(block_size)-1));

    % SSB OFDM demodulation
    SSB_CPs = ones(1, 4)*N_CP;
    grid = utils.demodulateOFDM(waveform(pss_idx:end), SSB_CPs, N_FFT);
    SSB = grid((-N_subcarriers_SSB/2+1:N_subcarriers_SSB/2)+N_FFT/2+SSB_subcarrier_offset, :);

    % Decoding SSS, searching for NID1
    NID1 = SSS.decodeSSS(SSB(57:183, 3).', NID2_, show_plots);
    fprintf("NID1 = %d\n", NID1);

    % PCI calculation
    cellid = 3 * NID1 + NID2_;
    fprintf("CellID = %d\n", cellid);
    
    if show_plots
        fprintf("Press ENTER to continue ...\n");
        pause;
    end
    
    
    %%%%%%%%%% Decode PBCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get PBCH indices
    [pbch_pos, pbch_dmrs_pos] = PBCH.getPBCHPosition(cellid);
    
    % Find correct issb for PBCH DM-RS
    fprintf(" -- Find i_SSB --\n");
    pbch_dmrs = SSB(pbch_dmrs_pos).';
    [iSSB, snr] = PBCH.decodePBCHDMRS(cellid, pbch_dmrs, ...
        show_plots);
    
    if show_plots
        fprintf("Press ENTER to continue ...\n");
        pause;
    end

    % Channel estimation
    fprintf(" -- Channel estimation and correction and PBCH decoding --\n");
    ref_pbch_dmrs = PBCH.generatePBCHDMRS(cellid, iSSB);
    h_est = pbch_dmrs .* conj(ref_pbch_dmrs);
    % Extend channel estimation to PBCH samples
    h_est = h_est + [0;0;0];
    h_est = h_est(:);
    
    % PBCH correction
    pbch = SSB(pbch_pos);
    pbch_eq = pbch ./ h_est;
    
    % Show PBCH constellation diagram
    if show_plots
        constellation_ref = [1+1i, -1+1i, -1-1i, 1-1i] ./ sqrt(2);
        figure;
        hold on;
        plot(pbch_eq, 'o');
        plot(constellation_ref, 'kx', 'LineWidth', 2, 'MarkerSize', 8);
        title('PBCH constellation diagram');
        xlabel('In-Phase');
        ylabel('Quadrature');
        fprintf("Press ENTER to continue ...\n");
        pause;
    end

    % Decode PBCH
    modulation_order = 2^2;
    modulation_vector = [3, 1, 0, 2];
    n_var = 1e-10;
    
    v = mod(iSSB, 8); % ??? 8, 4 or L_SSB?
    pbchBits_ = comm.internal.qam.demodulate(pbch_eq, 2^2, 'custom', ...
            [3 1 0 2], 1, 'bit', 1e-10);
    pbch_bits = nrPBCHDecode(pbch_eq, cellid, v, 1e-2);
    
    % Decode BCH
    fprintf(" -- BCH and MIB decoding --\n");
    polarListLength = 8;
    [~, BCH_CRC, tr_block, SFN_4_LSB, n_half_frame, kSSB_MSB] = ...
        BCH.decodeBCH(pbch_bits, L_SSB, cellid);
    
    % Display the BCH CRC and ssb index
    fprintf("BCH CRC: %d\n", BCH_CRC);
    fprintf("SSB index: %d\n", v);
    
    if BCH_CRC ~= 0
        fprintf("Error detected\n");
        continue;
    end
    
    % Decode MIB
    MIB = BCH.decodeMIB(tr_block, SFN_4_LSB, kSSB_MSB);

    fprintf("MIB:\n");
    disp(MIB);
    fprintf("Press ENTER to continue ...\n");
    pause;
    
    
    %%%%%%%%%% Check if CORESET#0 is present %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if all(MIB.kSSB > 23)
        % Another SSB with CORESET#0 information can be found according to 3GPP 38.213 13
        fprintf("No CORESET#0 present\n");
        continue;
    end
    
    
    %%%%%%%%%% Demodulate OFDM resource grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(" -- Demodulate OFDM resource grid --\n");
    
    % Estimate and correct timing offset to the frame start
    timing_offset = utils.estimateTimingOffset(sample_rate, N_sym_long, ...
        N_sym, SSB_case, fc, pss_idx-1, iSSB);
    waveform_frame = waveform(1+timing_offset:end);
    
    % Estimate and correct frequency offset to the common resource block
    SCS_common = MIB.SubcarrierSpacingCommon * 1e3;
    f_shift = MIB.kSSB * SCS_common;
    waveform_frame = waveform_frame .* ...
        exp(1j*2*pi*f_shift/sample_rate*(0:length(waveform_frame)-1));
    
    msb_idx = floor(MIB.PDCCHConfigSIB1/16);   % 4 MSB of PDCCHConfigSIB1
    lsb_idx = mod(MIB.PDCCHConfigSIB1, 16);    % 4 LSB of PDCCHConfigSIB1
    SCS_pair = [SCS_SSB, SCS_common];
    min_channel_BW = min(available_channel_BWs);
    
    % Get CORESET#0 information
    [N_RB_CORESET, N_sym_CORESET, CORESET_RB_offset, pattern] = ...
        PDCCH.getCORESET0Resources(msb_idx, SCS_pair, min_channel_BW, ...
        MIB.kSSB);
    
    % Demodulate the resource grid
    N_slots_in_signal = floor(length(waveform_frame) / ...
        (N_sym_long+13*N_sym));
    N_slots_per_frame = N_subframes_per_frame * ...
        utils.getSlotsPerSubframe(SCS_common);
    repeat = min(N_slots_in_signal, N_slots_per_frame);
    N_CPs = repmat([N_CP_long, N_CP*ones(1, 13)], 1, repeat);   % Minimum number of symbols to cover CORESET#0
    c0 = CORESET_RB_offset + 10*SCS_SSB/SCS_common;             % CORESET RB offset from carrier center (+padding)
    N_RB_CORESET_min = 2*max(c0, N_RB_CORESET-c0);              % Minimum number of RBs to cover CORESET#0
    N_subcarriers = N_RB_CORESET_min * N_subcarriers_per_RB;    
    subcarriers = (-N_subcarriers/2+1:N_subcarriers/2)+N_FFT/2;
    resource_grid = utils.demodulateOFDM(waveform_frame, N_CPs, N_FFT);
    resource_grid = resource_grid(subcarriers, :);
    
    figure;
    imagesc(abs(resource_grid));
    axis xy;
    title("Frame resource grid");
    xlabel("symbols");
    ylabel("subcarriers");
    fprintf("Press ENTER to continue ...\n");
    pause;
    
    
    %%%%%%%%%% Find PDCCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(" -- Find PDCCH --\n");
    
    % Find Type0-PDCCH monitoring occasions. 3GPP 28.213 13
    [n0, nC, is_occasion, frame_offset] = ...
        PDCCH.getPDCCH0MonitoringOccasions(lsb_idx, iSSB, SCS_pair, ...
        pattern, N_sym_CORESET, MIB.NFrame);
    
    monitoring_symbols = [n0, n0+1]*N_symbols_per_slot + nC + ...
        [1:N_symbols_per_slot].' + ...
        frame_offset*N_slots_per_frame*N_symbols_per_slot;
    monitoring_symbols = monitoring_symbols(:).';
    monitoring_subcarriers = (N_RB_CORESET_min - 20*SCS_SSB/SCS_common)*...
        N_subcarriers_per_RB/2 + CORESET_RB_offset*N_subcarriers_per_RB+...
        [1:N_RB_CORESET*12];
    monitoring_grid = resource_grid(monitoring_subcarriers, ...
        monitoring_symbols);
    
    % Get PDCCH parameters
    pdcch = PDCCH.getPDCCHParameters(SCS_SSB, min_channel_BW, MIB, ...
        iSSB, cellid);
    
    % Create a carrier that spans the wole CORESET#0 BWP
    c0Carrier = Carrier;
    c0Carrier.SubcarrierSpacing = MIB.SubcarrierSpacingCommon;
    c0Carrier.N_RB_start = pdcch.N_RB_start;
    c0Carrier.N_RB = pdcch.N_RB;
    c0Carrier.NSlot = pdcch.SearchSpace.SlotOffset;
    c0Carrier.NFrame = MIB.NFrame;
    c0Carrier.NCellID = cellid;
    
    
end