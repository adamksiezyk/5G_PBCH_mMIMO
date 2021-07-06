%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_dmrs(ncellid, issb)
%   generate PBCH DM-RS for given ncellid and issb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = generate_dmrs(ncellid, issb)
    %[38.211 7.4.1.4.1]
    c_init = 2^11 * (issb+1) * (fix(ncellid/4) + 1) + 2^6 * (issb+1) + (mod(ncellid, 4));
    c_seq = generate_pr_sequence(c_init, 2*144);
    %dmrs sequence
    r = zeros(1, 144);
    for n = 1:144
        r(n) = (1 - 2*c_seq(2*(n-1)+1))/sqrt(2) + 1j*(1 - 2*c_seq(2*(n-1)+2))/sqrt(2);
    end
end