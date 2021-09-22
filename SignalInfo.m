classdef SignalInfo
    %SIGNALINFO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fs;         % Represents the sample rate in Hz
        fc;         % Represents the carrier frequency in Hz
        BW;         % Represents the signal bandwidth in Hz, 5e6, 10e6, ...
        N_FFT;      % Represents the FFT size: 512, 1024, ...
        SCS;        % Represents the Subcarrier Spacing in Hz: 15e3, 30e3, 60e3 or 120e3
        SSB = SSB;  % Represents the Synchronization Signal Block
    end
    
    properties (SetAccess = private)
        N_subframes_per_frame;  % Represents the number of subframes in a frame
        N_symbols_per_slot;     % Represents the number of symbols in a slot
        N_subcarriers_per_RB;   % Represents the number of subcarriers in a Resource Block
        N_CP_long;              % Represents the number of samples per long CP
        N_CP;                   % Represents the number of samples per short CP
        N_sym_long;             % Represents the number of samples per long symbol
        N_sym;                  % Represents the number of samples per short symbol
        N_RBs;                  % Represents the number of RBs in the BW
        SSB_frequency_offset;   % Represents the SSB offset from the center frequency in Hz
    end
    
    methods
        function val = get.N_subframes_per_frame(obj)
            val = 10;
        end
        
        function val = get.N_symbols_per_slot(obj)
            val = 14;
        end
        
        function val = get.N_subcarriers_per_RB(obj)
            val = 12;
        end
        
        function val = get.N_CP_long(obj)
            [val, ~] = utils.getCPLength(obj.SCS, obj.fs);
        end
        
        function val = get.N_CP(obj)
            [~, val] = utils.getCPLength(obj.SCS, obj.fs);
        end
        
        function val = get.N_sym_long(obj)
            val = obj.N_FFT + obj.N_CP_long;
        end
        
        function val = get.N_sym(obj)
            val = obj.N_FFT + obj.N_CP;
        end
        
        function val = get.N_RBs(obj)
            val = utils.getRBAmount(obj.SCS, obj.BW);
        end
        
        function val = get.SSB_frequency_offset(obj)
            val = obj.SSB.subcarrier_offset * obj.SCS;
        end
        
        function val = get_SSB_subcarriers(obj)
            SSB_start =  obj.N_FFT/2 - obj.SSB.N_subcarriers_SSB/2 + ...
                obj.SSB.subcarrier_offset;
            val = (0:obj.SSB.N_subcarriers_SSB-1) + SSB_start;
        end
    end
end

