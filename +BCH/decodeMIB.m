function MIB = decodeMIB(BCCH_BCH_msg, SFN_4_LSB, kSSB_MSB)
%DECODEBCCH Decode the Master Information Block according to 3GPP 38.331
%6.2.2
% Inputs:
%   BCCH_BCH_msg    : a vector representing the BCCH-BCH message
%   SFN_4_LSB       : a vector representing the 4 LSBs of the system frame
%   number
%   kSSB_MSB        : a bit representing the kSSB MSB
% Outputs:
%   MIB     : a structure representing the Master Information Block

    k_SSB = kSSB_MSB * 2^4;     % MSB of a 5 bit long kSSB sequence
    commonSCSs = [15 30];
    
    MIB.NFrame = bi2de([BCCH_BCH_msg(2:7); SFN_4_LSB] .','left-msb');
    MIB.SubcarrierSpacingCommon = commonSCSs(BCCH_BCH_msg(8) + 1);
    MIB.k_SSB = k_SSB + bi2de(BCCH_BCH_msg(9:12).','left-msb');
    MIB.DMRSTypeAPosition = 2 + BCCH_BCH_msg(13);
    MIB.PDCCHConfigSIB1 = bi2de(BCCH_BCH_msg(14:21).','left-msb');
    MIB.CellBarred = BCCH_BCH_msg(22);
    MIB.IntraFreqReselection = BCCH_BCH_msg(23);
end

