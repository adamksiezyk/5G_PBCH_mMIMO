classdef SSB
    %SSB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SSB_case;           % Represents the SSB case: "Case A", "Case B", ...
        subcarrier_offset;  % Represents the SSB offset from the center subcarrier
        L_SSB;              % Represents the maximum number of SSBs in a SS Burst
    end
    
    methods(Static)
        function val = N_symbols_SSB()
            % Represents the number of symbols in a SSB
            val = 4;
        end
        
        function val = N_subcarriers_SSB()
            % Represents the number of subcarriers in a SSB
            val = 240;
        end
        
        function val = N_symbols_PSS()
            % Represents the number of symbols in a PSS
            val = 1;
        end
        
        function val = N_subcarriers_PSS()
            % Represents the number of subcarriers in a PSS
            val = 127;
        end
        
        function val = N_symbols_SSS()
            % Represents the number of symbols in a SSS
            val = 1;
        end
        
        function val = N_subcarriers_SSS()
            % Represents the number of subcarriers in a SSS
            val = 127;
        end
        
        function val = PSS_symbols()
            % Represents the symbol indices in a SSB where a PSS is presents
            val = [1];
        end
        
        function val = PSS_subcarriers()
            % Represents the subcarrier indices in a SSB where a PSS is presents
            val = [57:183];
        end
        
        function val = SSS_symbols()
            % Represents the symbol indices in a SSB where a SSS is presents
            val = [3];
        end
        
        function val = SSS_subcarriers()
            % Represents the subcarrier indices in a SSB where a SSS is presents
            val = [57:183];
        end
        
        function val = PBCH_position(cell_id)
            % Returns the RE indices in a SSB where a PBCH is present
            val = PBCH.getPBCHPosition(cell_id);
        end
        
        function val = PBCH_DMRS_position(cell_id)
            % Returns the RE indices in a SSB where a PBCH is present
            val = PBCH.getPBCHDMRSPosition(cell_id);
        end
    end
end

