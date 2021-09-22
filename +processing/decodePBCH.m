function [MIBs, PBCH_bits, iSSBs] = decodePBCH(signal_info, SSBs, ...
    cell_ids, show_plots_)
%DECODEPBCH Summary of this function goes here
% Inputs:
%   signal_info     : a SignalInfo object
%   SSBs            : a 3D matrix representing the SSB grids
%   cell_ids        : a vector representing the decoded cell IDs
%   show_plots_     : a boolean, if true plots are shown
% Outputs:
%   MIBs        : a vector representing the decoded MIBs
%   PBCH_bits   : a matrix representing the demodulated PBCH bits
%   iSSBs       : a vector that represents the detected iSSBs

    if nargin < 3
        show_plots = false;
    else
        show_plots = show_plots_;
    end
    
    SSB_amount = size(SSBs, 1);
    MIBs = BCH.MIB.empty(SSB_amount, 0);
    iSSBs = zeros(1, SSB_amount);
    PBCH_bits = zeros(SSB_amount, 864);
    for i = 1:SSB_amount
        cell_id = cell_ids(i);
        SSB_grid = SSBs(i, :, :);
        
        % Get PBCH indices
        [PBCH_pos, PBCH_DMRS_pos] = PBCH.getPBCHPosition(cell_id);

        % Find correct issb for PBCH DM-RS
        fprintf(" -- Find i_SSB --\n");
        PBCH_DMRS = SSB_grid(PBCH_DMRS_pos).';
        [iSSBs(i), ~] = PBCH.decodePBCHDMRS(cell_id, PBCH_DMRS, ...
            show_plots);

        if show_plots
            fprintf("Press ENTER to continue ...\n");
            pause;
        end

        % Channel estimation
        fprintf(" -- Channel estimation and correction and PBCH decoding --\n");
        ref_PBCH_DMRS = PBCH.generatePBCHDMRS(cell_id, iSSBs(i));
        h_est = PBCH_DMRS .* conj(ref_PBCH_DMRS);
        % Extend channel estimation to PBCH samples
        h_est = h_est + [0;0;0];
        h_est = h_est(:);

        % PBCH correction
        PBCH_grid = SSB_grid(PBCH_pos);
        PBCH_eq = PBCH_grid ./ h_est;

        % Show PBCH constellation diagram
        if show_plots
            constellation_ref = [1+1i, -1+1i, -1-1i, 1-1i] ./ sqrt(2);
            figure;
            hold on;
            plot(PBCH_eq, 'o');
            plot(constellation_ref, 'kx', 'LineWidth', 2, 'MarkerSize', 8);
            title('PBCH constellation diagram');
            xlabel('In-Phase');
            ylabel('Quadrature');
            fprintf("Press ENTER to continue ...\n");
            pause;
        end

        % Decode PBCH
        v = mod(iSSBs(i), signal_info.SSB.L_SSB);
        PBCH_bits(i, :) = nrPBCHDecode(PBCH_eq, cell_id, v, 1e-2).';

        % Check BCH CRC
        fprintf(" -- BCH decoding --\n");
        [~, BCH_CRC, tr_block, SFN_4_LSB, n_half_frame, kSSB_MSB] = ...
            BCH.decodeBCH(PBCH_bits(i, :), signal_info.SSB.L_SSB, cell_id);
        fprintf("CRC = %d\n", BCH_CRC);
        if BCH_CRC ~= 0
            fprintf("Error detected\n");
            continue;
        end
        
        % Decode MIB
        MIBs(i) = BCH.decodeMIB(tr_block, SFN_4_LSB, kSSB_MSB);

        fprintf("MIB:\n");
        disp(MIBs(i));
        if show_plots
            fprintf("Press ENTER to continue ...\n");
            pause;
        end
    end
end

