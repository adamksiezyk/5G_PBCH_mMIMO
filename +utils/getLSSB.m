function L_SSB = getLSSB(fc)
%GETLSSB Returns the maximum number of candidate SSBs within a SS burst (in
%one half frame) based on 3GPP TS 38.213 4.1
% Inputs:
%   fc  : a number representing the carrier frequency in Hz
% Outputs:
%   L_SSB   : a number representing the maximum number of candidate SSBs
%   within a SS burst
    if fc <= 3e9
        L_SSB = 4;
    elseif fc <= 6e9
        L_SSB = 8;
    else
        L_SSB = 64;
    end
end

