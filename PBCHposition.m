%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PBCHposition(ncellid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pbch_pos = PBCHposition(ncellid)
    ssb_len = 12*20;
    v = mod(ncellid, 4);
    %[38.211 7.4.3.1-1]
    dmrs_pos = [ ssb_len:4:ssb_len+236, ...
                 ssb_len*2:4:ssb_len*2+44, ...
                 ssb_len*2+192:4:ssb_len*2+236, ...
                 ssb_len*3:4:ssb_len*3+236];

    pbch = [1 2 3; 0 2 3; 0 1 3; 0 1 2].';
    %and add missing elements 
    pbch_pos = dmrs_pos + pbch(:, v+1);
    pbch_pos = pbch_pos(:) + 1;
end