function [pss_indices, synthetic_SSBs] = synthesis(waveform, SCS_SSB, BW, ...
    fc, fs, N_FFT, threshold, SSB_subcarrier_offset, show_plots)
%SYNTHESIS Detects the SSBs in the waveform and returns the syntesised SSB
%spectrum
% Inputs:
%   waveform                : a column vector that represents the 5G
%   waveform
%   SCS_SSB                 : a number that represents the SSB SCS in Hz
%   BW                      : a number that represents the signal bandwidth
%   in Hz
%   fc                      : a number that represents the carrier
%   frequency in Hz
%   fs                      : a number that represents the sampling
%   frequency in Hz
%   N_FFT                   : a number that represents the FFT size
%   threshold               : a number that represents the PSS detection
%   threshold, if equal to 0 the best PSS is chosen
%   SSB_subcarrier_offset   : a number that represents the SSB subcarrier
%   offset from the center subcarrier
%   show_plots              : a boolean, if true plots are shown
% Outputs:
%   pss_indices     : a vector representing the starting index of an SSB
%   synthetic_SSBs  : a n x N_FFT x 4 matrix representing the SSBs
%   spectrum. 1-dim - SSB id, 2-dim - subcarriers, 3-dim - symbols

    % 5G parameters
    N_subframes_per_frame = 10;
    N_symbols_per_slot = 14;
    N_subcarriers_per_RB = 12;
    N_subcarriers_SSB = 240;
    N_symbols_SSB = 4;
    N_subcarriers_PSS = 127;
    N = length(waveform);
    [N_CP_long, N_CP] = utils.getCPLength(SCS_SSB, fs);
    N_sym_long = N_FFT + N_CP_long;
    N_sym = N_FFT + N_CP;
    N_RBs = utils.getRBAmount(SCS_SSB, BW);
    L_SSB = utils.getLSSB(fc);

    % Plot the signal PSD and Spectrogram
    if show_plots
        fprintf(" -- Signal PSD and Spectrigram --\n");
        figure;
        pwelch(waveform, 2048, 2048-1024, 2048, fs, 'centered');

        figure;
        f = fs/2048 * (-1024:1023);
        spectrogram(waveform, 2048, 1024, f, fs);
        fprintf("Press ENTER to continue ...\n");
        pause;
    end

    % Find SSB in frequency domain
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

    % Estimation and correct the integer CFO
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

    % Estimation and correct the fractional CFO
    fprintf(" -- Fractional CFO estimation and correction --\n");
    ssb_frequency_offset = SSB_subcarrier_offset * SCS_SSB;
    [factional_CFO, CP_corr] = utils.estimateFractionalCFO(waveform, ...
        fs, SCS_SSB, ssb_frequency_offset);
    fprintf("Fractional CFO = %.2f\n", factional_CFO);
    waveform = waveform .* exp(-1j*2*pi * factional_CFO/fs *(0:N-1));

    if show_plots
        figure;
        plot(abs(CP_corr));
        title('Signal autocorrelation using CP of SSB');
        xlabel('Sample index n');
        fprintf("Press ENTER to continue ...\n");
        pause;
    end

    % Find SSB position in time domain and detect PSS
    fprintf(" -- Find SSB position in time domain --\n");
    [pss_indices, NID2_list] = PSS.detectAndDecodePSS(waveform, N_FFT, N_CP, ...
        threshold, SSB_subcarrier_offset, show_plots);
    if show_plots
        fprintf("Press ENTER to continue ...\n");
        pause;
    end

    % Process each SSB
    synthetic_SSBs = zeros(length(pss_indices), N_FFT, 4);
    for i = 1:length(pss_indices)
        %%%%%%%%%%% Deocode Cell ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(" -- Decoding Cell ID --\n");
        pss_idx = pss_indices(i);
        NID2 = NID2_list(i);
        fprintf("NID2 = %d\n", NID2);

        % Second CFO estimation using PSS
        n1 = pss_idx:pss_idx+N_CP-1;
        n2 = n1 + N_FFT;
        pss_CFO = fs/N_FFT * mean(angle(conj(waveform(n1)) .* ...
            waveform(n2)) / (2*pi));

        % Current SSB correction
        block_size = (pss_idx:pss_idx+4*N_sym-1);
        waveform(block_size) = waveform(block_size) .* exp(-1j*2*pi * ...
            pss_CFO/fs *(0:length(block_size)-1));

        % SSB OFDM demodulation
        SSB_CPs = ones(1, 4)*N_CP;
        grid = utils.demodulateOFDM(waveform(pss_idx:end), SSB_CPs, N_FFT);
        SSB = grid((-N_subcarriers_SSB/2+1:N_subcarriers_SSB/2)+N_FFT/2+SSB_subcarrier_offset, :);

        % Decoding SSS, searching for NID1
        NID1 = SSS.decodeSSS(SSB(57:183, 3).', NID2, show_plots);
        fprintf("NID1 = %d\n", NID1);

        % PCI calculation
        cellid = 3 * NID1 + NID2;
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
        [iSSB, ~] = PBCH.decodePBCHDMRS(cellid, pbch_dmrs, ...
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

        v = mod(iSSB, L_SSB);
        pbch_bits = nrPBCHDecode(pbch_eq, cellid, v, 1e-2);

        % Check BCH CRC
        fprintf(" -- BCH decoding --\n");
        polarListLength = 8;
        [~, BCH_CRC, ~, ~, ~, ~] = BCH.decodeBCH(pbch_bits, L_SSB, cellid);
        fprintf("CRC = %d\n", BCH_CRC);
        if BCH_CRC ~= 0
            fprintf("Error detected\n");
            continue;
        end

        %%%%%%%%%%% Synthesise PBCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3GPP 38.211 7.4.3, (resource allocation and scaling factor)
        synthetic_SSB = zeros(N_subcarriers_SSB, 4);
        synthetic_SSB(57:183, 1) = 1 * PSS.generatePSS(NID1);
        synthetic_SSB(57:183, 3) = 2 * SSS.generateSSS(NID1, NID2).';
        synthetic_SSB(pbch_pos) = 3 * nrPBCH(pbch_bits, cellid, v);
        synthetic_SSB(pbch_dmrs_pos) = 4 * PBCH.generatePBCHDMRS(cellid, iSSB);

        if show_plots
            figure;
            subplot(121);
            imagesc(abs(SSB));
            axis xy;
            xticklabels([0, 1, 2, 3]);
            xlabel("Symbols");
            ylabel("Subcarriers");
            title("Detected SSB");
            subplot(122);
            imagesc(abs(synthetic_SSB));
            axis xy;
            xticklabels([0, 1, 2, 3]);
            xlabel("Symbols");
            ylabel("Subcarriers");
            title("Synthesised SSB");
            fprintf("Press ENTER to continue ...\n");
            pause;
        end
        SSB_start = N_FFT/2 - N_subcarriers_SSB/2 + SSB_subcarrier_offset;
        synthetic_SSBs(i, (0:N_subcarriers_SSB-1)+SSB_start, :) = ...
            synthetic_SSB;
    end
end

