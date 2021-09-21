classdef Carrier
    %CARRIER Represents a 5G carrier
    
    properties
        SubcarrierSpacing {mustBeInteger}           % Represents the Subcarrier Spacing
        N_RB_start {mustBeInteger}                  % Represents the start RB index
        N_RB {mustBeInteger}                        % Represents the RB number
        NSlot {mustBeInteger}                       % Represents the current slot
        NFrame {mustBeInRange(NFrame, 0, 1023)}     % Represents the current System Frame Number
        NCellID {mustBeInRange(NCellID, 0, 1023)}   % Represents the cell ID
    end
    
    properties (SetAccess = private)
        %SymbolsPerSlot Number of OFDM symbols in a slot = 14
        SymbolsPerSlot;

        %SlotsPerSubframe Number of slots in a 1 ms subframe
        %   The value is one of {1, 2, 4, 8, 16} depending on
        %   SubcarrierSpacing values {15, 30, 60, 120, 240}, respectively.
        SlotsPerSubframe;

        %SlotsPerFrame Number of slots in a 10 ms frame
        %   The value is one of {10, 20, 40, 80, 160} depending on
        %   SubcarrierSpacing values {15, 30, 60, 120, 240}, respectively.
        SlotsPerFrame;
    end
    
    methods
        function val = get.SymbolsPerSlot(obj)
            % The number of OFDM symbols in a slot
            val = 14;
        end

        function val = get.SlotsPerSubframe(obj)
            % The number of slots in a subframe depends on the subcarrier spacing
            val = double(obj.SubcarrierSpacing) / 15;
        end

        function val = get.SlotsPerFrame(obj)
            % The number of slots in a frame depends on the subcarrier spacing
            val = 10*(double(obj.SubcarrierSpacing) / 15);
        end
    end
end

