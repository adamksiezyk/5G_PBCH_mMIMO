function [MIB, tr_block, iSSB, HFR] = decodePBCH(signal_info, SSB_grid, ...
    cell_id, show_plots_)
%DECODEPBCH Decodes the Physical Broadcast Channel.
% Inputs:
%   signal_info     : a SignalInfo object
%   SSBs            : a 3D matrix representing the SSB grids
%   cell_id         : a vector representing the decoded cell IDs
%   show_plots_     : a boolean, if true plots are shown
% Outputs:
%   MIB         : a vector representing the decoded MIBs
%   tr_block    : a matrix representing the BCH transport block
%   iSSB        : a vector that represents the detected iSSBs
%   HFR         : a logical representing the half frame bit

    if nargin < 3
        show_plots = false;
    else
        show_plots = show_plots_;
    end
    
    % Get PBCH indices
    PBCH_pos = PBCH.getPBCHPosition(cell_id);
    PBCH_DMRS_pos = PBCH.getPBCHDMRSPosition(cell_id);

    % Find correct issb for PBCH DM-RS
    fprintf(" -- Find i_SSB --\n");
    PBCH_DMRS = SSB_grid(PBCH_DMRS_pos).';
    [iSSB, ~] = PBCH.decodePBCHDMRS(cell_id, PBCH_DMRS, show_plots);

    if show_plots
        fprintf("Press ENTER to continue ...\n");
        pause;
    end

    % Channel estimation
    fprintf(" -- Channel estimation and correction and PBCH decoding --\n");
    ref_PBCH_DMRS = PBCH.generatePBCHDMRS(cell_id, iSSB);
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
    v = mod(iSSB, signal_info.SSB.L_SSB);
    PBCH_bits = PBCH.decodePBCH(PBCH_eq, cell_id, v, 1e-2).';

    % Check BCH CRC
    fprintf(" -- BCH decoding --\n");
    [~, BCH_CRC, tr_block, SFN_4_LSB, HFR, kSSB_MSB] = ...
        BCH.decodeBCH(PBCH_bits, signal_info.SSB.L_SSB, cell_id);
    fprintf("CRC = %d\n", BCH_CRC);
    if BCH_CRC ~= 0
        fprintf("Error detected, CRC != 0\n");
        MIB = cell(0);
    else
        % Decode MIB
        MIB = BCH.decodeMIB(tr_block, SFN_4_LSB, kSSB_MSB);

        fprintf("MIB:\n");
        disp(MIB);
        if show_plots
            fprintf("Press ENTER to continue ...\n");
            pause;
        end
    end
end

