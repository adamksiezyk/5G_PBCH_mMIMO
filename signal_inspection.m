clear variables; close all;
%% Load waveform and initial parameters
fprintf("Loading waveform and initial parameters\n");

show_plots = true;
signal_choice = 4;
if (signal_choice == 1)
    load('signals/2021-05-21_11-05-13_GPS_Fc3440.0_G50.0_Bw50.0_Fn000_ch2_spline_61.44MHz.mat');
    sample_rate = 61.44e6;
    %waveform = waveform(1:25.0e-3*sample_rate).';    % First 2.5 frames
    waveform = waveform.';
    fc = 3440e6;
    bw = 50e6;
    nfft = 2048;
    scs = 30e3;
    threshold = 2e5;
    ssb_pos = -507; % SSB frequency position
    int_CFO = 0;    % Integer Center Frequency Offset
elseif (signal_choice == 2)
    load('signals/Fn001_chan1.mat');
    sample_rate = 61.44e6;
    waveform = chan1;
    fc = 3440e6;
    bw = 50e6;
    nfft = 2048;
    scs = 30e3;
    threshold = 2e5;
    ssb_pos = -507; % SSB frequency position
    int_CFO = 0;    % Integer Center Frequency Offset
elseif (signal_choice == 3)
    load('Signals/NR_DL_2159.91MHz_10MHz.mat'); %ncellid = 440
    fc = 2159.91e6;
    bw = 10e6;
    nfft = 1024;
    sample_rate = 15360000;
    scs = 15e3;
    threshold = 60; %
    ssb_pos = -48;  % SSB frequency position
    int_CFO = 0;    % Integer Center Frequency offset
elseif (signal_choice == 4)
    load('signals/matlab_5g_downlink_5MHzBW_30kHzSCS.mat');
    sample_rate = 15360000;
    waveform = tmp_waveform.';
    fc = 3e9;
    bw = 10e6;
    nfft = 512;
    scs = 30e3;
    threshold = 8;
    ssb_pos = 0;
    int_CFO = 0;
end

N = length(waveform);
[N_CP_long, N_CP] = utils.getCPLength(scs, sample_rate);
N_sym_long = nfft + N_CP_long;
N_sym = nfft + N_CP;
N_RB = utils.getRBAmount(scs, bw);
N_SSB = 240;
N_PSS = 127;
L_SSB = utils.getLSSB(fc);

%% Plot the signal PSD and Spectrogram
if show_plots
    fprintf("Signal PSD and Spectrigram\n");
    figure;
    pwelch(waveform, 2048, 2048-1024, 2048, sample_rate, 'centered');

    figure;
    f = sample_rate/2048 * (-1024:1023);
    spectrogram(waveform, 2048, 1024, f, sample_rate);
    fprintf("Press ENTER to continue ...\n");
    pause;
end

%% Signal autocorelation
fprintf("Signal autocorelation\n");
N_corr = min(N-N_sym);%, 40e-3*sample_rate);
corrCP = zeros(1, N_corr);
for n = 1:N_corr
    n1 = n:n+N_CP-1; 
    n2 = n1 + nfft;
    a1 = abs(sum( waveform(n1) .* conj(waveform(n2))));
    b1 = abs(sum( waveform(n1) .* conj(waveform(n1)) + ...
        sum(waveform(n2) .* conj(waveform(n2)))));
    corrCP(n) = (a1*a1) / (b1*b1);
end
if show_plots
    figure;
    plot(corrCP);
    title('Signal autocorrelation using CP');
    xlabel('Sample index n');
    fprintf("Press ENTER to continue ...\n");
    pause;
end

%% Fractional CFO estimation and correction
fprintf("Fractional CFO estimation and correction\n");
% find npeak max peaks
npeak = 50;
range = [4e-5, 5e-5];
cp_pos = int64.empty;
corr_copy = corrCP;
idx = 1;
while (true)
    [val, pos] = max(corr_copy);
    if (val == 0 || idx == npeak)
        break;
    end
    prev_symbol = max(1, pos-N_sym);
    next_symbol = min(length(corr_copy), pos+N_sym);
    % Take only if symbol has a "neighbour"
    
    % VERSION WITH RANGE
    % if (val < range(1) || val > range(2) || ...
    %         corrCP(prev_symbol) < range(1) || corrCP(prev_symbol) > range(2) || ...
    %         corrCP(next_symbol) < range(1) || corrCP(next_symbol) > range(2))
    
    % VERSION WITH ABOUT SAME VALUE
    if (corrCP(prev_symbol) < val-0.1*val || ...
            corrCP(prev_symbol) > val+0.1*val || ...
            corrCP(next_symbol) < val-0.1*val || ...
            corrCP(next_symbol) > val+0.1*val)
        % Clean this symbol to ensure that we don't take it again
        corr_copy(prev_symbol:next_symbol) = zeros(1, ...
            next_symbol-prev_symbol+1);
        continue;
    end
    
    cp_pos(idx) = pos;
    % Clean this symbol to ensure that we don't take it again
    corr_copy(prev_symbol:next_symbol) = zeros(1, ...
        next_symbol-prev_symbol+1);
    idx = idx + 1;
end
cp_pos = sort(cp_pos); 

% calculate CFO for choosen symbols
indv_CFO = zeros(1, length(cp_pos));
for n = 1:length(cp_pos) 
    n1 = cp_pos(n):cp_pos(n)+N_CP-1;
    n2 = n1 + nfft;
    indv_CFO(n) = sample_rate/nfft * mean(angle(conj(waveform(n1)) .* ...
        waveform(n2)) / (2*pi));
end

% correct signal with mean of different CFO values
frc_CFO = mean(indv_CFO);
fprintf("Fractional CFO = %.2f\n", frc_CFO);
waveform = waveform .* exp(-1j*2*pi * frc_CFO/sample_rate *(0:N-1));

if show_plots
    figure;
    plot(cp_pos, indv_CFO, 'o-');
    title('Estimated fractional CFO');
    xlabel('sample index n');
    ylabel('[Hz]');

    hold on;
    yline(frc_CFO, 'k-');
    legend("individual fractional CFO", "mean fractional CFO");
    fprintf("Press ENTER to continue ...\n");
    pause;
end

%% Finding SSB in frequency domain
if (exist('ssb_pos', 'var') ~= 1)
    fprintf("Find SSB in frequency domian\n");
    raster = utils.getSSBRaster(fc, bw, scs);
    search_results = findSSB(waveform, nfft, raster);
    [max_val, max_ind] = max(search_results(:, 1));
    ssb_pos = raster(max_ind);

    if show_plots
        figure;
        hold on;
        stem(raster, search_results(:, 1));
        plot(ssb_pos, max_val, 'kx', 'LineWidth', 2, 'MarkerSize', 8);
        title("SSB frequency position search results");
        ylabel("Maximum of xcorr");
        xlabel("Subcarrier offset");
        legend("correlation result", "SSB position");
        fprintf("Press ENTER to continue ...\n");
        pause;
    end
end

%% Integer CFO estimation and correction
if (exist('int_CFO', 'var') ~= 1)
    fprintf("Integer CFO estimation and correction\n");
    offsets = -10:10;
    search_results = findSSB(waveform, nfft, ssb_pos+offsets);

    [max_val, max_ind] = max(search_results(:, 1));
    int_CFO = offsets(max_ind);
    fprintf("Integer CFO: %d\n", int_CFO);
    waveform = waveform .* exp(-1j*2*pi * int_CFO/nfft * (0:N-1));

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

%% Finding SSB position in time domain and detecting PSS
fprintf("Find SSB position in time domain\n");
[pss_indexes, NID2] = PSS.detectAndDecodePSS(waveform, nfft, threshold, ...
    ssb_pos, show_plots);
pss_indexes = pss_indexes - N_sym + 1;
if show_plots
    fprintf("Press ENTER to continue ...\n");
    pause;
end

%% Decoding each PBCH, MIB
NID1 = zeros(1, length(pss_indexes));
cellid = zeros(1, length(pss_indexes));
issb = zeros(1, length(pss_indexes));
snr = zeros(1, length(pss_indexes));
for i = 1:length(pss_indexes)
    fprintf("Decoding Cell ID\n");
    pss_pos = pss_indexes(i);
    fprintf("NID2 = %d\n", NID2(i));
    
    % Second CFO estimation using PSS
    n1 = pss_pos:pss_pos+N_CP-1;
    n2 = n1 + nfft;
    pss_CFO = sample_rate/nfft * mean(angle(conj(waveform(n1)) .* ...
        waveform(n2)) / (2*pi));

    % Current SSB correction
    block_size = (pss_pos:pss_pos+4*N_sym-1);
    waveform(block_size) = waveform(block_size) .* exp(-1j*2*pi * ...
        pss_CFO/sample_rate *(0:length(block_size)-1));

    SSB = zeros(240, 4);
    for n = 1:4
        pos = pss_pos + N_sym*(n-1);
        samples = waveform(pos+N_CP:pos+N_sym-1);
        spectrum = fftshift(fft(samples)) / sqrt(nfft);
        SSB(:, n) = spectrum(...
            nfft/2-N_SSB/2+1+ssb_pos:nfft/2+N_SSB/2+ssb_pos);
    end

    % Decoding SSS, searching for NID1
    NID1(i) = SSS.decodeSSS(SSB(57:183, 3).', NID2(i), show_plots);
    fprintf("NID1 = %d\n", NID1(i));

    % PCI calculation
    cellid(i) = 3 * NID1(i) + NID2(i);
    fprintf("CellID = %d\n", cellid(i));
    
    if show_plots
        fprintf("Press ENTER to continue ...\n");
        pause;
    end
    
    % Get PBCH indices
    [pbch_pos, pbch_dmrs_pos] = PBCH.getPBCHPosition(cellid(i));
    
    % Find correct issb for PBCH DM-RS
    fprintf("Find i_SSB\n");
    pbch_dmrs = SSB(pbch_dmrs_pos).';
    [issb(i), snr(i)] = PBCH.decodePBCHDMRS(cellid(i), pbch_dmrs, ...
        show_plots);
    
    if show_plots
        fprintf("Press ENTER to continue ...\n");
        pause;
    end

    % Channel estimation
    fprintf("Channel estimation and correction and PBCH decoding\n");
    ref_pbch_dmrs = PBCH.generatePBCHDMRS(cellid(i), issb(i));
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
    
    v = mod(issb(i), 8); % ??? 8, 4 or L_SSB?
    pbchBits_ = comm.internal.qam.demodulate(pbch_eq, 2^2, 'custom', ...
            [3 1 0 2], 1, 'bit', 1e-10);
    pbchBits = nrPBCHDecode(pbch_eq, cellid(i), v, 1e-2);
    
    % Decode BCH
    fprintf("BCH and MIB decoding\n");
    polarListLength = 8;
    [~, BCH_CRC, tr_block, SFN_4_LSB, n_half_frame, kSSB_MSB] = ...
        BCH.decodeBCH(pbchBits, L_SSB, cellid(i));
    
    % Display the BCH CRC and ssb index
    fprintf("BCH CRC: %d\n", BCH_CRC);
    fprintf("SSB index: %d\n", v);
    
    % Decode MIB
    MIB(i) = BCH.decodeMIB(tr_block, SFN_4_LSB, kSSB_MSB);

    % Display the MIB structure
    fprintf("MIB:\n");
    disp(MIB(i));
    fprintf("Press ENTER to continue ...\n");
    pause;
end

