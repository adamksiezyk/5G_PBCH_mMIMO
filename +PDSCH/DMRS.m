classdef DMRS
    %DMRS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NSCID;             % TS 38.211 Section 7.4.1.1
        DMRSPortSet = 0;
        DMRSConfigurationType = 1;
        DMRSLength = 1;
        NumCDMGroupsWithoutData;
        DMRSTypeAPosition;
        DMRSAdditionalPosition;
        DMRSReferencePoint = 'PRB0';
    end
end

