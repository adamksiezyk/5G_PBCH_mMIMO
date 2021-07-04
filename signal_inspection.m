clear variables; close all;
%% Load waveform and initial parameters
fprintf("Loading waveform and initial parameters\n");

signal_choice = 1;
if (signal_choice == 1)
    load('signals/2021-05-21_11-05-13_GPS_Fc3440.0_G50.0_Bw50.0_Fn000_ch2_spline_61.44MHz.mat');
    sample_rate = 61.44e6;
    waveform = waveform(1:20.5e-3*sample_rate).';    % First 2.5 frames
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
sym_len_long = nfft + N_CP_long;
sym_len = nfft + N_CP;
N_RB = get_n_rb(scs, bw);
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
N_corr = min(N-sym_len, 40e-3*sample_rate);
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
cp_pos = zeros(1, npeak);
corr_copy = corrCP;
for n = 1:npeak
    [aa, cp_pos(n)] = max(corr_copy);
    %don't take indexes from beyond array length
    start = max(1, cp_pos(n)-sym_len);
    stop = min(length(corr_copy), cp_pos(n)+sym_len);
    %clean this symbol to ensure that we don't take it again
    corr_copy(start:stop) = zeros(1, stop-start+1);
end
cp_pos = sort(cp_pos); 
% calculate CFO for choosen symbols
indv_CFO = zeros(1, npeak);
for n = 1:npeak 
    n1 = cp_pos(n):cp_pos(n)+N_CP-1;
    n2 = n1 + nfft;
    indv_CFO(n) = sample_rate/nfft * mean(angle(conj(waveform(n1)) .* waveform(n2)) / (2*pi));
end

figure;
plot(cp_pos(1:npeak), indv_CFO(1:npeak), 'o-');
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
ssb_offset = offsets(max_ind);
int_CFO = ssb_offset;
fprintf("Integer CFO: %d\n", int_CFO);
waveform = waveform .* exp(-1j*2*pi * int_CFO/nfft * (0:N-1));

figure;
plot(offsets, search_results(:, 1));
title("SSB frequency offset search results");
ylabel("Maximum of xcorr");
xlabel("Subcarrier offset");
fprintf("Press ENTER to continue ...\n");
pause;

%% Find SSB position in time domain
fprintf("Find SSB position in time domain\n");
[pss_pos, NID2] = pss_decoding(waveform, nfft, threshold, ssb_pos);
