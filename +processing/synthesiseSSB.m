function synthetic_grid = synthesiseSSB(signal_info, PBCH_bits, ...
    NID2, NID1, iSSB, show_plots_)
%SYNTHESIS Detects the SSBs in the waveform and returns the syntesised SSB
%spectrum. 3GPP 38.211 7.4.3, (resource allocation and scaling factor)
% Inputs:
%   signal_info     : a SignalInfo object
%   PBCH_bits       : a 2D matrix that represents the PBCH bits (1-dim -
%   PBCH id, 2-dim - bits)
%   NID2            : a vector that represents the detected NID2s
%   NID1            : a vector that represents the detected NID1s
%   iSSB            : a vector that represents the detected iSSBs
%   show_plots_     : a boolean, if true plots are shown
% Outputs:
%   synthetic_SSB   : a N_FFT x N_symbols_SSB matrix representing the SSBs
%   spectrum. 1-dim - SSB id, 2-dim - subcarriers, 3-dim - symbols

    if nargin < 6
        show_plots = false;
    else
        show_plots = show_plots_;
    end

    cell_id = 3*NID1 + NID2;
    v = mod(iSSB, signal_info.SSB.L_SSB);
    PBCH_pos = signal_info.SSB.PBCH_position(cell_id);
    PBCH_DMRS_pos = signal_info.SSB.PBCH_DMRS_position(cell_id);

    % Synthesize the SSB
    fprintf(" -- SSB synthesis --\n");
    synthetic_SSB = zeros(SSB.N_subcarriers_SSB, SSB.N_symbols_SSB);
    synthetic_SSB(57:183, 1) = 1 * PSS.generatePSS(NID2);
    synthetic_SSB(57:183, 3) = 2 * SSS.generateSSS(NID1, NID2).';
    synthetic_SSB(PBCH_pos) = 3 * nrPBCH(PBCH_bits, cell_id, v);
    synthetic_SSB(PBCH_DMRS_pos) = 4 * PBCH.generatePBCHDMRS(cell_id, iSSB);

    if show_plots
        figure;
        imagesc(abs(synthetic_SSB));
        axis xy;
        xticklabels([0, 1, 2, 3]);
        xlabel("Symbols");
        ylabel("Subcarriers");
        title("Synthesised SSB");
        fprintf("Press ENTER to continue ...\n");
        pause;
    end

    % Place the SSB in the spectrum
    synthetic_grid = zeros(signal_info.N_FFT, signal_info.SSB.N_symbols_SSB);
    synthetic_grid(signal_info.get_SSB_subcarriers, :) = synthetic_SSB;
end

