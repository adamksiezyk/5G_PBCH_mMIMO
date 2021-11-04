function soft_bits = demodulateQAM(symbols, modulation_order, n_var)
%DEMODULATEQAM Summary of this function goes here
%   Detailed explanation goes here
    switch modulation_order
        case 2
            modulation = 'BPSK';
        case 4
            modulation = 'QPSK';
    end
    soft_bits = nrSymbolDemodulate(symbols, modulation, n_var);
end

