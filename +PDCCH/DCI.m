classdef DCI
    %DCI Downlink Control Information
    
    properties
        FrequencyDomainResourceAssignment {mustBeInteger};
        TimeDomainResourceAssignment {mustBeInteger};
        VRBToPRBMapping {mustBeInteger};
        ModulationAndCodingSchemes {mustBeInteger};
        RedundancyVersion {mustBeInteger};
        SystemInformationIndicator {mustBeInteger};
        Reserved {mustBeInteger};
        
        field_sizes = struct();
    end
    
    methods
        function obj = DCI(N_CORESET_RB)
            obj.field_sizes.FrequencyDomainResourceAssignment = ceil(...
                log2(N_CORESET_RB*(N_CORESET_RB+1)/2));
            obj.field_sizes.TimeDomainResourceAssignment = 4;
            obj.field_sizes.VRBToPRBMapping = 1;
            obj.field_sizes.ModulationAndCodingSchemes = 5;
            obj.field_sizes.RedundancyVersion = 2;
            obj.field_sizes.SystemInformationIndicator = 1;
            obj.field_sizes.Reserved = 15;
        end
        
        function val = getDCISize(obj)
            values = structfun(@(x)x, obj.field_sizes);
            val = sum(values);
        end
    end
end

