function N_slots = getSlotsPerSubframe(SCS)
%GETSLOTSPERSUBFRAME Returns the number of slots per subframe
% Inputs:
%   SCS     : a number representing the subcarrier spacing in Hz
% Outputs:
%   N_slots : a number representing the number of slots per subframe

    mu = log2(SCS / 15e3);
    N_slots = 10*2^mu;
end

