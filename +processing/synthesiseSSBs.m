function synthetic_SSBs = synthesiseSSBs(signal_info, PBCH_bits, ...
    NID2_list, NID1_list, iSSBs, show_plots_)
%SYNTHESIS Detects the SSBs in the waveform and returns the syntesised SSB
%spectrum. 3GPP 38.211 7.4.3, (resource allocation and scaling factor)
% Inputs:
%   signal_info     : a SignalInfo object
%   PBCH_bits       : a 2D matrix that represents the PBCH bits (1-dim -
%   PBCH id, 2-dim - bits)
%   NID2_list       : a vector that represents the detected NID2s
%   NID1_list       : a vector that represents the detected NID1s
%   iSSBs           : a vector that represents the detected iSSBs
%   show_plots_     : a boolean, if true plots are shown
% Outputs:
%   pss_indices     : a vector representing the starting index of an SSB
%   synthetic_SSBs  : a n x N_FFT x 4 matrix representing the SSBs
%   spectrum. 1-dim - SSB id, 2-dim - subcarriers, 3-dim - symbols

    if nargin < 6
        show_plots = false;
    else
        show_plots = show_plots_;
    end

    % Process each SSB
    SSB_amount = size(PBCH_bits, 1);
    synthetic_SSBs = zeros(SSB_amount, signal_info.N_FFT, ...
        SSB.N_symbols_SSB);
    for i = 1:SSB_amount
        bits = PBCH_bits(i, :);
        NID2 = NID2_list(i);
        NID1 = NID1_list(i);
        cell_id = 3*NID1 + NID2;
        iSSB = iSSBs(i);
        v = mod(iSSB, signal_info.SSB.L_SSB);
        PBCH_pos = signal_info.SSB.PBCH_position(cell_id);
        
        fprintf(" -- SSB synthesis --\n");
        synthetic_SSB = zeros(SSB.N_subcarriers_SSB, SSB.N_symbols_SSB);
        synthetic_SSB(57:183, 1) = 1 * PSS.generatePSS(NID2);
        synthetic_SSB(57:183, 3) = 2 * SSS.generateSSS(NID1, NID2).';
        synthetic_SSB(PBCH_pos) = 3 * nrPBCH(bits.', cell_id, v);
%         synthetic_SSB(pbch_dmrs_pos) = 4 * PBCH.generatePBCHDMRS(cellid, iSSB);

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
        synthetic_SSBs(i, signal_info.get_SSB_subcarriers, :) = ...
            synthetic_SSB;
    end
end

