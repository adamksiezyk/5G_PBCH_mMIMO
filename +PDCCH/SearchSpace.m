classdef SearchSpace
    %SEARCHSPACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        StartSymbolWithinSlot {mustBeInRange(StartSymbolWithinSlot, 0, 13)} % Represents the index of the start symbol in a slot
        SlotPeriod {mustBeInteger}                                          % Represents the slot period
        SlotOffset {mustBeInteger}                                          % Represents the slot offset
        Duration {mustBeInRange(Duration, 2, 3)}                            % Represents the amount of symbols
        AggregationLevels = [1 2 4 8 16]                                    % Represents the aggregation levels 3GPP 38.211 Table 7.3.2.1-1
        NumCandidates = [0 0 8 4 1]                                         % Represents the candidates 3GPP 38.213 Table 10.1-1
    end
end

