function SSB_start_symbols = getSSBStartSymbols(SSB_case, frequency)
%GETSSBSTARTSYMBOLS Returns the SSB starting symbols indices. 3GPP 38.213
%4.1
% Inputs:
%   SSB_case    : a string representing the SSB case
%   frequency   : a number represetnting the frequency in Hz
% Outputs:
%   SSB_start_symbols   : a vector representing all possible SSB starting
%   symbols
    
    switch SSB_case
        case {"Case A", "Case C"}
            if frequency <= 3e9
                n = [0, 1];
            elseif frequency <= 6e9
                n = [0, 1, 2, 3];
            else
                error("SSB Case A or C can not have frequency > 6 GHz")
            end
            m = 14;
            i = [2, 8];
        case "Case B"
            if frequency <= 3e9
                n = [0];
            elseif frequency <= 6e9
                n = [0, 1];
            else
                error("SSB Case B can not have frequency > 6 GHz")
            end
            m = 28;
            i = [4, 8, 16, 20];
        case "Case D"
            if frequency <= 6e9
                error("SSB Case D can not have frequency > 6 GHz")
            end
            n = [0:18];
            m = 28;
            i = [4, 8, 16, 20];
        case "Case E"
            if frequency <= 6e9
                error("SSB Case A can not have frequency > 6 GHz")
            end
            n = [0:8];
            m = 56;
            i = [8, 12, 16, 20, 32, 36, 40, 44];
        otherwise
            error("Wrong SSB case provided");
    end
    
    SSB_start_symbols = i.' + m*n;
    SSB_start_symbols = SSB_start_symbols(:);
end

