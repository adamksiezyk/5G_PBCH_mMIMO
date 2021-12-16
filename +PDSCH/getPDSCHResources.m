function [pdsch, K_0] = getPDSCHResources(dci, N_RB, DMRSTypeAPosition, pattern)
%GETPDSCHRESOURCES Summary of this function goes here
%   Detailed explanation goes here

    pdsch = PDSCH.PDSCHParameters();
    
    % Configure PDSCH from DCI message
    [L_RBs, RB_start] = getRBAllocation(N_RB, dci.FrequencyDomainResourceAssignment);
    pdsch.RB_allocated = RB_start + [0:(L_RBs-1)];
    pdsch.VRBToPRBInterleaving = dci.VRBToPRBMapping;
    pdsch.RNTI = 65535; % SI-RNTI 
    
    % Select applicable PDSCH time domain resource allocation table based
    % on RNTI and CORESET pattern (TS 38.214 Table 5.1.2.1.1-1)
    % SI-RNTI Type0 Common
    time_tables = getPDSCHTimeAllocationTables(); 
    time_table = time_tables{pattern};
    time_table = time_table(time_table.TDDIndex==dci.TimeDomainResourceAssignment,:);
    time_allocations = time_table(time_table.DMRSTypeAPosition==DMRSTypeAPosition,:);
    
    K_0 = time_allocations.K_0;
    
    pdsch.mapping_type = time_allocations.PDSCHMappingType;
    pdsch.sym_allocated = [time_allocations.S time_allocations.L];
   
    % Configure PDSCH from MIB / SSB
%     pdsch.NID = []; % Cell identity (TS 38.211 Section 7.3.1.1)
    pdsch.DMRS.DMRSTypeAPosition = DMRSTypeAPosition;
%     pdsch.DMRS.NIDNSCID = []; % Cell identity (TS 38.211 Section 7.4.1.1)
    
    % Configure PDSCH from other relevant rules
    % TS 38.214 Section 5.1.3.1
    pdsch.modulation = 'QPSK';
    % TS 38.211 Section 7.4.1.1
    pdsch.DMRS.NSCID = 0;
    % TS 38.214 Section 5.1.6.2
    pdsch.NumLayers = 1; 
    pdsch.DMRS.DMRSPortSet = 0;
    pdsch.DMRS.DMRSConfigurationType = 1;
    pdsch.DMRS.DMRSLength = 1;
    L = pdsch.sym_allocated(2);
    if (L==2)
        pdsch.DMRS.NumCDMGroupsWithoutData = 1;
    else
        pdsch.DMRS.NumCDMGroupsWithoutData = 2;
    end
    if (strcmpi(pdsch.mapping_type,'A'))
        pdsch.DMRS.DMRSAdditionalPosition = 2;
    else % 'B'
        switch L
            case {2,4}
                pdsch.DMRS.DMRSAdditionalPosition = 0;
            case 7
                pdsch.DMRS.DMRSAdditionalPosition = 1;
        end
    end
    
    pdsch.DMRS.DMRSReferencePoint = 'PRB0';
    
    pdsch.EnablePTRS = false; % TS 38.214 Section 5.1.6.3
end

function [L_RBs, RB_start] = getRBAllocation(N_RB, RIV)
    L_RBs = floor(RIV / N_RB) + 1;
    RB_start = RIV - ((L_RBs - 1) * N_RB);
    if (L_RBs > N_RB - RB_start)
        L_RBs = N_RB - L_RBs + 2;
        RB_start = N_RB - 1 - RB_start;
    end
end

function tables = getPDSCHTimeAllocationTables()
% 3GPP 38.214 5.1.2.1
    % K_0: slot offset from DCI slot scheduling PDSCH
    % S: Start symbol of allocated PDSCH
    % L: Length in symbols of allocated PDSCH
            
    duplicate = @(x)reshape(repmat(x,2,1),[],1);
    repeat = @(x,n)repmat(x.',n,1);

    % TS 38.214 Table 5.1.2.1.1-2
    TDDIndex = duplicate(0:15);
    DMRSTypeAPosition = repeat([2 3],16);
    %                  row   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
    PDSCHMappingType = duplicate(["A" "A" "A" "A" "A" "B" "B" "B" "B" "B" "B" "A" "A" "A" "B" "B"]);
    K_0 = repeat([0 0],16);
    S = [[ 2  3  2  3  2  3  2  3  2  3  9 10  4  6].'; duplicate([5  5  9 12  1  1  2  4  8])];
    L = [[12 11 10  9  9  8  7  6  5  4  4  4  4  4].'; duplicate([7  2  2  2 13  6  4  7  4])];
    defaultA = table(TDDIndex,DMRSTypeAPosition,PDSCHMappingType,K_0,S,L);

    % TS 38.214 Table 5.1.2.1.1-4
    TDDIndex = [duplicate(0:14); 15];
    DMRSTypeAPosition = [repeat([2 3],15); NaN];
    PDSCHMappingType = [duplicate(["B" "B" "B" "B" "B" "B" "B" "B" "B" "B" "B" "B" "B" "A" "B"]); ""];
    K_0 = [duplicate([0  0  0  0  0  1  1  0  0  0  0  0  0  0  1]); NaN];
    S = [duplicate([2  4  6  8 10  2  4  2  4  6  8 10  2]);  2;  3;  2;  2; NaN];
    L = [duplicate([2  2  2  2  2  2  2  4  4  4  4  4  7]); 12; 11;  4;  4; NaN];
    defaultB = table(TDDIndex,DMRSTypeAPosition,PDSCHMappingType,K_0,S,L);

    % TS 38.214 Table 5.1.2.1.1-5
    TDDIndex = [duplicate(0:4); 5; 6; duplicate(7:15)];
    DMRSTypeAPosition = [repeat([2 3],5); NaN; NaN; repeat([2 3],9)];
    PDSCHMappingType = [duplicate(["B" "B" "B" "B" "B"]); ""; ""; duplicate(["B" "B" "B" "B" "B" "B" "A" "A" "A"])];
    K_0 = [repeat([0 0],5); NaN; NaN; repeat([0 0],9)];
    S = [duplicate([2  4  6  8 10]); NaN; NaN; duplicate([2  4  6  8 10  2]);  2;  3; duplicate([0  2])];
    L = [duplicate([2  2  2  2  2]); NaN; NaN; duplicate([4  4  4  4  4  7]); 12; 11; duplicate([6  6])];
    defaultC = table(TDDIndex,DMRSTypeAPosition,PDSCHMappingType,K_0,S,L);

    tables = {defaultA defaultB defaultC};
end