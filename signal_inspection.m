clear variables; close all;
%% Load waveform and initial parameters
fprintf("Loading waveform and initial parameters\n");

signal_choice = 1;
if (signal_choice == 1)
    load('signals/2021-05-21_11-05-13_GPS_Fc3440.0_G50.0_Bw50.0_Fn000_ch2_spline_61.44MHz.mat');
    sample_rate = 61.44e6;
    waveform = waveform(1:25.0e-3*sample_rate).';    % First 2.5 frames
    fc = 3440e6;
    bw = 50e6;
    nfft = 2048;
    scs = 30e3;
    threshold = 2e5;
    % ssb_offset = ;
elseif (signal_choice == 2)
    load('Signals/NR_DL_2159.91MHz_10MHz.mat'); %ncellid = 440
    fc = 2159.91;
    bw = 10e6;
    nfft = 1024;
    sample_rate = 15360000;
    scs = 15e3;
    threshold = 60; %
    % ssb_offset = -48;
end

N = length(waveform);
[N_CP_long, N_CP] = get_n_cp(scs, sample_rate);
N_sym_long = nfft + N_CP_long;
N_sym = nfft + N_CP;
N_RB = get_n_rb(scs, bw);
N_SSB = 240;
N_PSS = 127;

%% Plot the signal PSD and Spectrogram
fprintf("Signal PSD and Spectrigram\n");
figure;
pwelch(waveform, 2048, 2048-1024, 2048, sample_rate, 'centered');

figure;
f = sample_rate/2048 * (-1024:1023);
spectrogram(waveform, 2048, 1024, f, sample_rate);
fprintf("Press ENTER to continue ...\n");
pause;

%% Signal autocorelation
fprintf("Signal autocorelation\n");
N_corr = min(N-N_sym, 40e-3*sample_rate);
corrCP = zeros(1, N_corr);
for n = 1:N_corr
    n1 = n:n+N_CP-1; 
    n2 = n1 + nfft;
    a1 = abs(sum( waveform(n1) .* conj(waveform(n2))));
    b1 = abs(sum( waveform(n1) .* conj(waveform(n1)) + sum(waveform(n2) .* conj(waveform(n2)))));
    corrCP(n) = (a1*a1) / (b1*b1);
end
figure;
plot(corrCP);
title('Signal autocorrelation using CP');
xlabel('Sample index n');
fprintf("Press ENTER to continue ...\n");
pause;

%% Fractional CFO estimation and correction
fprintf("Fractional CFO estimation and correction\n");
% find npeak max peaks
npeak = 50;
range = [4e-5, 5e-5];
cp_pos = int64.empty;
corr_copy = corrCP;
% for n = 1:npeak
%     [aa, cp_pos(n)] = max(corr_copy);
%     %don't take indexes from beyond array length
%     start = max(1, cp_pos(n)-N_sym);
%     stop = min(length(corr_copy), cp_pos(n)+N_sym);
%     %clean this symbol to ensure that we don't take it again
%     corr_copy(start:stop) = zeros(1, stop-start+1);
% end
idx = 1;
while (true)
    [val, pos] = max(corr_copy);
    if (val == 0 || idx == npeak)
        break;
    end
    prev_symbol = max(1, pos-N_sym);
    next_symbol = min(length(corr_copy), pos+N_sym);
    % Take only if symbol has a "neighbour"
    if (val < range(1) || val > range(2) || ...
            corrCP(prev_symbol) < range(1) || corrCP(prev_symbol) > range(2) || ...
            corrCP(next_symbol) < range(1) || corrCP(next_symbol) > range(2))
        % Clean this symbol to ensure that we don't take it again
        corr_copy(prev_symbol:next_symbol) = zeros(1, next_symbol-prev_symbol+1);
        continue;
    end
    cp_pos(idx) = pos;
    % Clean this symbol to ensure that we don't take it again
    corr_copy(prev_symbol:next_symbol) = zeros(1, next_symbol-prev_symbol+1);
    idx = idx + 1;
end
cp_pos = sort(cp_pos); 

% calculate CFO for choosen symbols
indv_CFO = zeros(1, length(cp_pos));
for n = 1:length(cp_pos) 
    n1 = cp_pos(n):cp_pos(n)+N_CP-1;
    n2 = n1 + nfft;
    indv_CFO(n) = sample_rate/nfft * mean(angle(conj(waveform(n1)) .* waveform(n2)) / (2*pi));
end

figure;
plot(cp_pos, indv_CFO, 'o-');
title('Estimated fractional CFO');
xlabel('sample index n');
ylabel('[Hz]');

% correct signal with mean of different CFO values
frc_CFO = mean(indv_CFO);
fprintf("Fractional CFO = %.2f\n", frc_CFO);
waveform = waveform .* exp(-1j*2*pi * frc_CFO/sample_rate *(0:N-1));
fprintf("Press ENTER to continue ...\n");
pause;

%% Find SSB in frequency domain
fprintf("Find SSB in frequency domian\n");
raster = ssb_raster(fc, bw, scs);
search_results = ssb_search(waveform, nfft, raster);
[~, max_ind] = max(search_results(:, 1));
ssb_pos = raster(max_ind);

figure;
plot(raster, search_results(:, 1));
title("SSB frequency position search results");
ylabel("Maximum of xcorr");
xlabel("Subcarrier offset");
fprintf("Press ENTER to continue ...\n");
pause;

%% Integer CFO estimation and correction
fprintf("Integer CFO estimation and correction\n");
offsets = -10:10;
search_results = ssb_search(waveform, nfft, ssb_pos+offsets);

[~, max_ind] = max(search_results(:, 1));
int_CFO = offsets(max_ind);
fprintf("Integer CFO: %d\n", int_CFO);
waveform = waveform .* exp(-1j*2*pi * int_CFO/nfft * (0:N-1));

figure;
plot(offsets, search_results(:, 1));
title("SSB frequency offset search results");
ylabel("Maximum of xcorr");
xlabel("Subcarrier offset");
fprintf("Press ENTER to continue ...\n");
pause;

%% Find SSB position in time domain and detect PSS
fprintf("Find SSB position in time domain\n");
[pss_indexes, NID2] = pss_detect_and_decode(waveform, nfft, threshold, ssb_pos);
pss_indexes = pss_indexes - N_sym + 1;
fprintf("Press ENTER to continue ...\n");
pause;

%% Detect SSS
NID1 = zeros(1, length(NID2));
cellid = zeros(1, length(NID2));
for i = 1:length(pss_indexes)
    pss_pos = pss_indexes(i);
    fprintf("NID2 = %d\n", NID2(i));
    
    % Second CFO estimation using PSS
    n1 = pss_pos:pss_pos+N_CP-1;
    n2 = n1 + nfft;
    pss_CFO = sample_rate/nfft * mean(angle(conj(waveform(n1)) .* waveform(n2)) / (2*pi));

    % Current SSB correction
    block_size = (pss_pos:pss_pos+4*N_sym-1);
    waveform(block_size) = waveform(block_size) .* exp(-1j*2*pi * pss_CFO/sample_rate *(0:length(block_size)-1));

    SSB = zeros(240, 4);
    for n = 1:4
        pos = pss_pos + N_sym*(n-1);
        samples = waveform(pos+N_CP:pos+N_sym-1);
        spectrum = fftshift(fft(samples)) / sqrt(nfft);
        SSB(:, n) = spectrum(nfft/2-N_SSB/2+1+ssb_pos:nfft/2+N_SSB/2+ssb_pos);
    end

    %searching for NID1
    NID1(i) = sss_decode(SSB(57:183, 3).', NID2(i));
    fprintf("NID1 = %d\n", NID1(i));

    %PCI calculation
    cellid(i) = 3 * NID1(i) + NID2(i);
    fprintf("CellID = %d\n", cellid(i));
    fprintf("Press ENTER to continue ...\n");
    pause;
    
    % ##################
    % ##################
    % ##################
    
    if(1)
    
    SS_block = SSB;
    do_FEQ = 1;
    
    %find correct issb for extracted SS block
    issb = pbchdmrs_decode( cellid(i), SS_block );
    pause

    %find indicates for PBCH
    pbch_pos = PBCHposition( cellid(i) );

    %channel estimation 
    ref_dmrs = generate_dmrs( cellid(i), issb );
    rx_dmrs = SS_block(PBCHDMRSposition( cellid(i) )).';
    hest = rx_dmrs .* conj(ref_dmrs);
    %extend channel estimation to PBCH samples
    hest = hest + [0;0;0];
    pbch_hest = hest(:);

    %PBCH correction
    pbch_sym = SS_block(pbch_pos);
    if(do_FEQ)
        pbch_eq = pbch_sym ./ pbch_hest;
    else
        pbch_eq = pbch_sym;
    end
    %show PBCH constellation
    figure; plot(pbch_eq, 'o');
    title('PBCH constellation'); xlabel('In-Phase'); ylabel('Quadrature');
    pause

    %MIB decoding using 5G toolbox functions
    v = mod(issb, 4);
    pbchBits = nrPBCHDecode(pbch_eq, cellid(i), v, 1e-2);
    
    polarListLength = 8;
    [~, crcBCH, trblk, sfn4lsb, nHalfFrame, msbidxoffset] = ...
        nrBCHDecode(pbchBits, polarListLength, 4, cellid(i));
    
    % Display the BCH CRC and ssb index
    disp([' BCH CRC: ' num2str(crcBCH)]);
    disp([' SSB index: ' num2str(v)]);

    k_SSB = msbidxoffset * 16;
    commonSCSs = [15 30];
    
    % Create a structure of MIB fields from the decoded MIB bits. The BCH
    % transport block 'trblk' is the RRC message BCCH-BCH-Message, consisting
    % of a leading 0 bit then 23 bits corresponding to the MIB
    mib.NFrame = bi2de([trblk(2:7); sfn4lsb] .','left-msb');
    mib.SubcarrierSpacingCommon = commonSCSs(trblk(8) + 1);
    mib.k_SSB = k_SSB + bi2de(trblk(9:12).','left-msb');
    mib.DMRSTypeAPosition = 2 + trblk(13);
    mib.PDCCHConfigSIB1 = bi2de(trblk(14:21).','left-msb');
    mib.CellBarred = trblk(22);
    mib.IntraFreqReselection = trblk(23);

    % Display the MIB structure
    disp(' BCH/MIB Content:')
    disp(mib);
    pause

    end    
end