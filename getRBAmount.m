function N_RB = getRBAmount(scs, bw)
%GETRBAMOUNT Calculates the number of resource blocks according to 3GPP TS
%38.101 5.3.2-1
% Inputs:
%   scs     : a number representing the subcarrier spacing in Hz
%   bw      : a number representing the bandwidth in Hz
% Outputs:
%   N_RB    : a number representing the number of resource blocks

    if (scs == 15e3)
        if (bw == 5e6)
            N_RB = 25;
        elseif (bw == 10e6)
            N_RB = 52;
        elseif (bw == 15e6)
            N_RB = 79;
        elseif (bw == 20e6)
            N_RB = 106;
        elseif (bw == 25e6)
            N_RB = 133;
        elseif (bw == 30e6)
            N_RB = 160;
        elseif (bw == 40e6)
            N_RB = 216;
        elseif (bw == 50e6)
            N_RB = 270;
        end
    elseif (scs == 30e3)
        if (bw == 5e6)
            N_RB = 11;
        elseif (bw == 10e6)
            N_RB = 24;
        elseif (bw == 15e6)
            N_RB = 38;
        elseif (bw == 20e6)
            N_RB = 51;
        elseif (bw == 25e6)
            N_RB = 65;
        elseif (bw == 30e6)
            N_RB = 78;
        elseif (bw == 40e6)
            N_RB = 106;
        elseif (bw == 50e6)
            N_RB = 133;
        elseif (bw == 60e6)
            N_RB = 162;
        elseif (bw == 80e6)
            N_RB = 217;
        elseif (bw == 90e6)
            N_RB = 245;
        elseif (bw == 100e6)
            N_RB = 273;
        end
    elseif (scs == 60e3)
        if (bw == 10e6)
            N_RB = 11;
        elseif (bw == 15e6)
            N_RB = 18;
        elseif (bw == 20e6)
            N_RB = 24;
        elseif (bw == 25e6)
            N_RB = 31;
        elseif (bw == 30e6)
            N_RB = 38;
        elseif (bw == 40e6)
            N_RB = 51;
        elseif (bw == 50e6)
            N_RB = 65;
        elseif (bw == 60e6)
            N_RB = 79;
        elseif (bw == 80e6)
            N_RB = 107;
        elseif (bw == 90e6)
            N_RB = 121;
        elseif (bw == 100e6)
            N_RB = 135;
        end
    end
end

