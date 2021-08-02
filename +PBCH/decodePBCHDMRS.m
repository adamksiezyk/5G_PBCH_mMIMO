function issb = decodePBCHDMRS(cellid, dmrs, show_plots_)
%DECODEPBCHDMRS Returns the SSB index based on the given Cell ID and PBCH
%DM-RS sequence
% Inputs:
%   cellid          : a number representing the Cell ID
%   dmrs            : a vector representing the PBCH DM-RS
%   show_plots_     : a boolean representing the show plots flag default is
%   true
% Outputs:
%   issb    : a number representing the decoded SSB index

    if nargin<3
        show_plots = true;
    else
        show_plots = show_plots_;
    end

    dmrs_estimated = zeros(1, 8);
    for issb = 0:7 %4 slots in SS burst and 2 half-frames
        gen_dmrs = PBCH.generatePBCHDMRS(cellid, issb);
        dmrs_estimated(issb+1) = sqrt(abs(sum(dmrs .* conj(gen_dmrs(1:end))) .^2));
    end
    
    [max_val, issb] = max(dmrs_estimated); 
    issb = issb - 1;
    
    if show_plots
        figure;
        hold on;
        stem(0:7, dmrs_estimated, 'o');
        plot(issb, max_val, 'kx', 'LineWidth', 2, 'MarkerSize', 8);
        title('DM-RS detection');
        xlabel('issb');
        legend("correlation result", "i_{SSB}");
    end
end