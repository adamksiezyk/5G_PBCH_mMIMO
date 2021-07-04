function raster = ssb_raster(fc, bw, scs)
%SSB_RASTER Constructs SSB raster in subcarrier offsets based on
%   3GPP TS 38.104 5.4.3.1-1
% Inputs:
%   fc     : a number representing the carrier frequency in Hz
%   bw     : a number representing the bandwidth in Hz
%   scs    : a number representing the subcarrier spacing in Hz
% Outputs:
%   raster : a vector representing the subcarrier offset where an SSB is
%            expected

    raster = zeros(1, 26638);
    i = 1;
    % For frequency range from 0 to 3000 MHz
    for n = 1:2499
        for m = [1 3 5]
            raster(i) = 1.2e6*n + 50e3*m;
            i = i+1;
        end
    end
    % For frequency range from 3000 to 24250 MHz
    for n = 0:14756
        raster(i) = 3000e6 + n*1.44e6;
        i = i+1;
    end
    % For frequency range from 24250 to 100000 MHz
    for n = 0:4383
        raster(i) = 24250.08e6 + n*17.28e6;
        i = i+1;
    end
    % Extract only values that fit in bandwidth
    raster = raster(raster > (fc - bw/2));
    raster = raster(raster < (fc + bw/2));
    raster = round((raster - fc) / scs);
end

