function dci = decodeDCI(dci_bits, empty_dci)
%DECODEDCI Summary of this function goes here
%   Detailed explanation goes here
    
    dci = empty_dci;
    
    field_sizes = structfun(@(x)x, dci.field_sizes);
    field_ranges = [[0; cumsum(field_sizes(1:end-1))]+1, cumsum(field_sizes)];
    field_ranges = num2cell(field_ranges, 2);
    
    get_value = @(x,y)bin2dec(char(x(y(1):y(2)) + '0'));
    values = cellfun(@(x)get_value(dci_bits.', x), field_ranges);
    values = num2cell(values);
    values = cell2struct(values, fieldnames(dci.field_sizes));
    
    dci.FrequencyDomainResourceAssignment = values.FrequencyDomainResourceAssignment;
    dci.TimeDomainResourceAssignment = values.TimeDomainResourceAssignment;
    dci.VRBToPRBMapping = values.VRBToPRBMapping;
    dci.ModulationAndCodingSchemes = values.ModulationAndCodingSchemes;
    dci.RedundancyVersion = values.RedundancyVersion;
    dci.SystemInformationIndicator = values.SystemInformationIndicator;
    dci.Reserved = values.Reserved;
end

