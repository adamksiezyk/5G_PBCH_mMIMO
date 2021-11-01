classdef CORESET
    %CORESET A Control Resource Set representiation according to 3GPP
    %38.211 7.3.2
    
    properties
        id {mustBeInteger};
        N_RB_start {mustBeInteger}  % Represents the index of the starting Resource Block
        N_RB {mustBeInteger}        % Represents the total number of Resource Blocks
        N_sym = 2;                  % a number that represents the time duration in OFDM symbols
        REG_bundle_size = 6;        % a number that represents the Resource Element Group Bundle size
        N_REGs_per_CCE = 6;         % a number that represents the number of REGs per CCE
        interleaved {mustBeNumericOrLogical} = true;
        interleaver_size = 2;
        N_shift = 0;
    end
    
    methods
        function f = getCCEMapping(obj)
            % Returns the REG bundle mapping. 3GPP 38.211 7.3.2
            
            % Construct the control-channel element which consists of 6
            % REGs (1 REG = 1 RB over 1 OFDM symbol
            N_REG = obj.N_RB * obj.N_sym;

            if obj.interleaved
                L = obj.REG_bundle_size;
                R = obj.interleaver_size;
                C = N_REG / (L*R);
                f = zeros(R*C, 1);
                for c = 0:C-1
                    for r = 0:R-1
                        j = c*R + r;
                        f(j+1) = mod(r*C + c + obj.N_shift, N_REG/L);
                    end
                end
            else
                L = 6;
                N_CCE = N_REG / obj.N_REGs_per_CCE;
                f = (0:(N_CCE-1)).';
            end            
        end
        
        function RB_indices = getREGBundleIndices(obj)
            % Returns the RB indices of the REG bundles. 3GPP 38.211 7.3.2
            RB_indices = reshape(obj.getREGIndices(), ...
                obj.REG_bundle_size, []);
        end
        
        function RB_indices = getREGIndices(obj)
            % Returns the RB indices of the REGs. REGs are numbered in a
            % time-first manner! 3GPP 38.211 7.3.2
            REGs = find(ones(obj.N_RB, obj.N_sym));
            RB_indices = reshape(reshape(REGs, [], obj.N_sym).', [], 1);
        end
    end
end

