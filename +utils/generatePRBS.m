function c = generatePRBS(M_PN, c_init)
%GENERATEPRBS Generates the Pseudo-Random Binary Sequence as defined in
%3GPP 38.211 5.2.1.
% Inputs:
%   M_PC    : a number representing the sequence length
%   c_init  : a number representing the initial value for sequence
%generation
% Outputs:
%   c   : a vector representing the generated pseudo-random binary sequence

    % Initialization
    Nc = 1600;
    x1 = zeros(1, M_PN+Nc);
    x1(1) = 1;
    x2 = de2bi(c_init, M_PN+Nc);
    
    % Generate 2 m-sequences
    for n = 32:M_PN+Nc
        x1(n) = mod(x1(n-28) + x1(n-31), 2);
        x2(n) = mod(x2(n-28) + x2(n-29) + x2(n-30) + x2(n-31), 2);
    end
    
    % Generate the pseudo-random sequence
    c = zeros(M_PN, 1);
    for n = 1:M_PN
        c(n) = mod(x1(n+Nc) + x2(n+Nc), 2);
    end
end

