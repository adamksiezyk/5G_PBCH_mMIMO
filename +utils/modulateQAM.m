function symbols = modulateQAM(bits, order)
%MODULATEQAM Quadrature Amplitude Modukation according to 3GPP 38.211 5.1
% Inputs:
%   bits    : a vector representing the sequence of bits
%   order   : a number (power of 2) representing the modulation order
% Outputs:
%   symbols     : a vector representing the modulated symbols
    
    if order == 2
        symbols = 1/sqrt(2) * [(1-2*bits) + 1j*(1-2*bits)];
    elseif order == 4
        b1 = bits(1:2:end);
        b2 = bits(2:2:end);
        symbols = 1/sqrt(2) * [(1-2*b1) + 1j*(1-2*b2)];
    end
end

