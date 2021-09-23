classdef MIB
    %MIB Represents the Master Information Block
    %   3GPP 38.331 6.2.2
    
    properties
        NFrame                      % This field defines the System Frame Number
        SubcarrierSpacingCommon     % This field defines the subcarrier spacing used for SIB1
        kSSB                        % This field defines the frequency domain offset between SSB and the overall resource block grid in number of subcarriers
        DMRSTypeAPosition           % This field defines the position of first DM-RS symbol for downlink (PDSCH) and uplink (PUSCH)
        PDCCHConfigSIB1             % This field is used to configure CORESET#0 and search space#0
        CellBarred                  % This field indicates whether or not UEs in the cell are allowed to access the cell; ‘barred’ indicates, the UEs are not allowed to access the cell
        IntraFreqReselection        % This field controls cell selection/reselection to intra-frequency cells when the highest ranked cell is barred
    end
end

