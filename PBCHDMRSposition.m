%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PBCHDMRSposition(ncellid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dmrs_pos = PBCHDMRSposition(ncellid)
    ssb_len = 12*20;
    v = mod(ncellid, 4);
    %[38.211 7.4.3.1-1]
    dmrs_pos = [ssb_len:4:ssb_len+236, ...
                ssb_len*2:4:ssb_len*2+44, ...
                ssb_len*2+192:4:ssb_len*2+236, ...
                ssb_len*3:4:ssb_len*3+236].';
    dmrs_pos = dmrs_pos + 1 + v;
end