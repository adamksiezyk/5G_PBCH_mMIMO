classdef PDSCHParameters
    %PDSCHPARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N_RB;                   % Bandwidthpart RB amount
        N_RB_start;             % Index of the starting RB
        RB_allocated;           % Vector of indices of the allocated RBs
        VRBToPRBInterleaving;
        RNTI = 65535;
        mapping_type;
        sym_allocated;          % A vector representing the allocated symbols
        DMRSTypeAPosition;
        modulation = 'QPSK';    % TS 38.214 Section 5.1.3.1
        NumLayers = 1;          % TS 38.214 Section 5.1.6.2
        DMRS = PDSCH.DMRS();
        EnablePTRS;             % TS 38.214 Section 5.1.6.3
    end
end

