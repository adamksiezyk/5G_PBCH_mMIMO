function pbch_dmrs_pos = getPBCHDMRSPosition(cellid)
%GETPBCHDMRSPOSITION Returns the PBCH DM-RS position in SSB based on the
%provided Cell ID
% Inputs:
%   cellid  : a number representing the Cell ID
% Outputs:
%   pbch_dmrs_pos   : a vector representing the PBCH DM-RS position in SSB

    ssb_len = 12*20;
    v = mod(cellid, 4);
    %[38.211 7.4.3.1-1]
    pbch_dmrs_pos = [ssb_len:4:ssb_len+236, ...
                ssb_len*2:4:ssb_len*2+44, ...
                ssb_len*2+192:4:ssb_len*2+236, ...
                ssb_len*3:4:ssb_len*3+236].';
    pbch_dmrs_pos = pbch_dmrs_pos + 1 + v;
end