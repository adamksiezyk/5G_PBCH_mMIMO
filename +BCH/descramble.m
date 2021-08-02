function tr_block = descramble(scr_block, cell_id, L_SSB)
%DESCRAMBLE descramble the scrambled block according to 3GPP 38.212 7.1.2
% Inputs:
%   scr_block   : a vector representing the scrambled block
%   cell_id     : a number representing the cell ID
%   L_SSB       : a number repersenting the maximum number of SSBs in a
%   half frame
% Outputs:
%   tr_block    : a vector represetning the descrambled transport block

    A = length(scr_block);
    c_init = cell_id;
    G = BCH.getPayloadInterleaverPattern();
    isScrambled = true(A, 1);
    SFN_2_3_LSB = [];
    jSFN = 0;
    jHRF = 10;
    jSSB = 11;
    for i = 0:A-1     % 0-based
        if ((i>0 && i<=6) || (i>=24 && i<=27)) % SFN
            % Bits with idx 1...6 are the 6 MSBs of the SFN in the
            % MIB (3GPP 38.331 6.2.2). Bits with idx 24...27 are the 4 LSBs
            % of the SFN conveyed in additional PBCH payload bits (7.1.1)
            if (i==25) || (i==26)
                % 2nd, 3rd LSB of SFN
                isScrambled(G(jSFN+1)+1, 1) = false;
                SFN_2_3_LSB = [SFN_2_3_LSB; G(jSFN+1)+1];
            end
            jSFN = jSFN+1;
        elseif (i==28)                  % HRF
            isScrambled(G(jHRF+1)+1, 1) = false;
        elseif (i>=29 && i<=31)         % Candidate SS/PBCH block index
            if (L_SSB==64)
                isScrambled(G(jSSB+1)+1, 1) = false;
            end
            jSSB = jSSB+1;
        end
    end
    
    if L_SSB == 64
        M = A - 6;
    elseif L_SSB == 20
        M = A - 5;
    elseif L_SSB == 10
        M = A - 4;
    else
        M = A - 3;
    end
    % Get v from scrBlk, which is not scrambled, just interleaved
    v = bi2de(scr_block(SFN_2_3_LSB).', 2, 'left-msb');
    c = utils.generatePRBS((v+1)*M, c_init);
    s = zeros(A, 1);
    s(isScrambled) = c(v*M+1:v*M+M);
    tr_block = mod(scr_block+s, 2);
end
