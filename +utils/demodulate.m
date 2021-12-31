function soft_bits = demodulate(symbols, modulation, n_var)
%DEMODULATE Summary of this function goes here
%   Detailed explanation goes here
    soft_bits = nrSymbolDemodulate(symbols, modulation, n_var);
end

