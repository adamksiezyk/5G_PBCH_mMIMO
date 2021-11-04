function N = getNPolarEncoded(K, E, n_max)
%GETN Returns the mother code block length for the specified number of
%input bits K, number of rate-matched output bits E and maximum value of
%n N_MAX. This process is the inverse of the one described in 3GPP TS 
%38.212 5.3.1.
% Inputs:
%   K       : a number representing the number of input bits
%   E       : a number representing the number of rate-matched output bits
%   n_max   : a number representing the maximum number of N
% Outputs:
%   N   : a number representing the mother code block length

    cl2e = ceil(log2(E));
    if (E <= (9/8) * 2^(cl2e-1)) && (K/E < 9/16)
        n1 = cl2e-1;
    else
        n1 = cl2e;
    end

    R_min = 1/8;
    n2 = ceil(log2(K/R_min));

    nMin = 5;
    n = max(min([n1 n2 n_max]), nMin);
    N = 2^n;
end

