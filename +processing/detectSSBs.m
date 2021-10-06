function [SSB_indices, SSBs, NID2_list, NID1_list] = detectSSBs(waveform, ...
    signal_info, threshold, show_plots_)
%DETECTSSBS Returns the SSBs detected in the given waveform
% Inputs:
%   waveform        : a column vector that represents the waveform
%   signal_info     : a SignalInfo object
%   threshold       : a number representing the PSS detection threshold, if
%   0, the best PSS is returned
%   show_plots_     : a boolean, if true plots are shown
% Outputs:
%   SSB_indices     : a vector representing the detected SSB indices
%   SSBs            : a 3D matrix representing the found SSBs (1-dim - SSB
%   id, 2-dim - subcarriers, 3-dim - symbols)
%   NID2_list       : a vector representing the detected NID2 list
%   NID1_list       : a vector representing the detected NID1 list

    if nargin < 4
        show_plots = false;
    else
        show_plots = show_plots_;
    end

    [SSB_indices, NID2_list] = PSS.detectAndDecodePSS(waveform, ...
        signal_info.N_FFT, signal_info.N_CP, threshold, ...
        signal_info.SSB.subcarrier_offset, show_plots);
    if show_plots
        fprintf("Press ENTER to continue ...\n");
        pause;
    end
    
    % Deocode Cell ID
    SSB_amount = length(SSB_indices);
    SSBs = zeros(SSB_amount, SSB.N_subcarriers_SSB, SSB.N_symbols_SSB);
    NID1_list = zeros(1, SSB_amount);
    for i = 1:SSB_amount
        SSB_idx = max(SSB_indices(i), 1);
        NID2 = NID2_list(i);

        % Second CFO estimation using PSS
        n1 = SSB_idx:SSB_idx+signal_info.N_CP-1;
        n2 = n1 + signal_info.N_FFT;
        PSS_CFO = signal_info.fs/signal_info.N_FFT * ...
            mean(angle(conj(waveform(n1)) .* waveform(n2)) / (2*pi));

        % Current SSB correction
        block = SSB_idx:SSB_idx+SSB.N_symbols_SSB*signal_info.N_sym-1;
        waveform(block) = waveform(block) .* ...
            exp(-1j*2*pi * PSS_CFO/signal_info.fs * (0:length(block)-1));

        % SSB OFDM demodulation
        SSB_CPs = ones(1, SSB.N_symbols_SSB) * signal_info.N_CP;
        grid = utils.demodulateOFDM(waveform(SSB_idx:end), SSB_CPs, ...
            signal_info.N_FFT);
        grid = grid((-SSB.N_subcarriers_SSB/2+1:SSB.N_subcarriers_SSB/2) + ...
            signal_info.N_FFT/2 + signal_info.SSB.subcarrier_offset, :);
        SSBs(i, :, :) = grid;
        

        % Decoding SSS, searching for NID1
        NID1_list(i) = SSS.decodeSSS(grid(...
            SSB.SSS_subcarriers, SSB.SSS_symbols).', NID2, show_plots);
        fprintf("NID1 = %d\n", NID1_list(i));

        if show_plots
            fprintf("Press ENTER to continue ...\n");
            pause;
        end
    end
end

