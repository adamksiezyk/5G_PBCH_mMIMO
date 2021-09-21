classdef SearchSpace
    %SEARCHSPACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        StartSymbolWithinSlot {mustBeInRange(StartSymbolWithinSlot, 0, 13)} % Represents the index of the start symbol in a slot
        SlotPeriod {mustBeInteger}                                          % Represents the slot period
        SlotOffset {mustBeInteger}                                          % Represents the slot offset
        Duration {mustBeInRange(Duration, 2, 3)}                            % Represents the amount of symbols
        NumCandidates {mustBeNumeric} = [0 0 8 4 1]                         % Represents the candidates TS 38.213 Table 10.1-1
    end
end

