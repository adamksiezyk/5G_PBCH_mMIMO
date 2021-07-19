function pbch_dmrs = generatePBCHDMRS(cellid, issb)
%GENERATEPBCHDMRS Generates the PBCH DM-RS sequence for the given Cell ID
%and SSB index based on 3GPP TS 38.211 7.4.1.4.1
% Inputs:
%   cellid      : a number representing the Cell ID
%   issb        : a number representing the SSB index
% Outputs:
%   pbch_dmrs   : a vector representing the PBCH DM-RS sequence

    c_init = 2^11 * (issb+1) * (fix(cellid/4) + 1) + 2^6 * (issb+1) + ...
        (mod(cellid, 4));
    c_seq = generate_pr_sequence(c_init, 2*144);
    pbch_dmrs = zeros(1, 144);
    for n = 1:144
        pbch_dmrs(n) = (1 - 2*c_seq(2*(n-1)+1))/sqrt(2) + ...
            1j*(1 - 2*c_seq(2*(n-1)+2))/sqrt(2);
    end
end