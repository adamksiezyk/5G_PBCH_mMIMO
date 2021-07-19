function issb = decodePBCHDMRS(cellid, dmrs)
%DECODEPBCHDMRS Returns the SSB index based on the given Cell ID and PBCH
%DM-RS sequence
% Inputs:
%   cellid  : a number representing the Cell ID
%   dmrs    : a vector representing the PBCH DM-RS
% Outputs:
%   issb    : a number representing the decoded SSB index

    dmrs_estimated = zeros(1, 8);
    for issb = 0:7 %4 slots in SS burst and 2 half-frames
        gen_dmrs = generatePBCHDMRS(cellid, issb);
        dmrs_estimated(issb+1) = sqrt(abs(sum(dmrs .* conj(gen_dmrs(1:end))) .^2));
    end
    
    figure;
    hold on;
    stem(0:7, dmrs_estimated, 'o');
    title('DM-RS detection');
    xlabel('issb');
    
    [max_val, issb] = max(dmrs_estimated); 
    issb = issb - 1;
    
    hold on;
    plot(issb, max_val, 'kx', 'LineWidth', 2, 'MarkerSize', 8);
    legend("correlation result", "i_{SSB}");
end