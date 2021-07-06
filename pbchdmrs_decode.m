%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pbchdmrs_decode(ncellid, sig_grid)
% Input:
%   ncellid - cell ID
%   sig_grid - extracted and demodulated SS block
% Output:
%   issb - SS block correct issb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function issb = pbchdmrs_decode(ncellid, sig_grid)
    
    %[38.211 7.4.3.1-1]
    dmrs_pos = PBCHDMRSposition(ncellid);
    
    for issb = 0:7 %4 slots in SS burst and 2 half-frames
        dmrs = sig_grid(dmrs_pos).';
        gen_dmrs = generate_dmrs(ncellid, issb);
        dmrs_estimated(issb+1) = sqrt(abs(sum(dmrs .* conj(gen_dmrs(1:end))) .^2));
    end
    
    figure; stem(0:7, dmrs_estimated, 'o');
    title('DM-RS detection'); xlabel('issb');
    [aa, issb] = max(dmrs_estimated); 
    issb = issb - 1;
end