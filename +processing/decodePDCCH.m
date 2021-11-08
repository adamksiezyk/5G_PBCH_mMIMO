function [dci, dci_crc] = decodePDCCH(grid, carrier, pdcch, ...
    pdcch_dmrs_indices, pdcch_indices, mon_slots, show_plots_)
%DECODEPDCCH Decodes the PDCCH and returns the DCI and CRC value.
% Inputs:
%   grid                : a matrix representing the resources grid
%   grid(subcariers, symbols)
%   carrier             : a Carrier class
%   pdcch               : a PDCCHParameters class
%   pdcch_dmrs_indices  : a cell representing the PDCCH DM-RS indices per
%   aggregation level and candidate
%   pdcch_indices       : a cell representing the PDCCH indices per
%   aggregation level and candidate
%   mon_slots           : a vector representing the monitoring slots
%   show_plots_         : a boolean, if true plots are shown
% Outputs:
%   dci         : a DCI class
%   dci_crs     : a number represeting the calculated CRC value

    if nargin < 7
        show_plots = false;
    else
        show_plots = show_plots_;
    end
    
    dci_crc = 1;
    for monitoring_slot = 0:length(mon_slots)-1
        % Get PDCCH candidates according to TS 38.213 Section 10.1
        slot_grid = grid(:, (1:carrier.SymbolsPerSlot) + ...
            carrier.SymbolsPerSlot*monitoring_slot, :);
        slot_grid = slot_grid/max(abs(slot_grid(:)));

        for al = 1:5
            for c = 1:pdcch.SearchSpace.NumCandidates(al)
                % Estimate channel using PDCCH DM-RS        
                pdcch_dmrs = slot_grid(pdcch_dmrs_indices{al}(:, c));
                pdcch_dmrs_ref = PDCCH.getPDCCHDMRS(carrier, pdcch, al);
                pdcch_dmrs_ref = pdcch_dmrs_ref(:, c);
                h_est = pdcch_dmrs .* conj(pdcch_dmrs_ref);
                h_est = h_est + [0, 0, 0];
                h_est = h_est(:);

                % Equalize PDCCH
                pdcch_symbols = slot_grid(pdcch_indices{al}(:, c));
                pdcch_symbols_eq = pdcch_symbols ./ h_est;

                % Decode PDCCH
                empty_dci = PDCCH.DCI(pdcch.CORESET.N_RB);
                N_dci = empty_dci.getDCISize();
                [dci_bits, dci_crc] = PDCCH.decodePDCCH(pdcch_symbols_eq, N_dci,...
                    pdcch.DMRS_scrambling_ID, pdcch.RNTI, 1e-2);

                % Decode DCI
                dci = PDCCH.decodeDCI(dci_bits, empty_dci);
                if dci_crc == 0
                    break;
                end
            end
            if dci_crc == 0
                break;
            end
        end
        if dci_crc == 0
            break;
        end
    end
    
    if dci_crc == 0 && show_plots
        constellation_ref = [1+1i, -1+1i, -1-1i, 1-1i] ./ sqrt(2);
        figure;
        hold on;
        plot(pdcch_symbols_eq, 'o');
        plot(constellation_ref, 'kx', 'LineWidth', 2, 'MarkerSize', 8);
        title('PDCCH constellation diagram');
        xlabel('In-Phase');
        ylabel('Quadrature');
        fprintf("Press ENTER to continue ...\n");
        pause;
        
        fprintf("PDCCH found in aggregation level: %d and candidate: %d\n", ...
            pdcch.SearchSpace.AggregationLevels(al), c);
        disp(dci);
        fprintf("Press ENTER to continue ...\n");
    end
end

