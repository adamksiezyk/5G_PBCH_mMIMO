function [scr_block, err_flag, BCCH_BCH_msg, SFN_4_LSB, HRF, kSSB_MSB] ...
    = decodeBCH(code_block, L_SSB, cell_id)
%DECODEBCH Decodes the BCH tramsport channel. This process is the inverse
%of BCH genereation described in 3GPP 38.212 7.1.
% Inputs:
%   code_block  : a vector representing the demodulated log-likelihood
%ratios softbits
%   L_SSB       : a number representing the maximal number of PBCH blocks
%in one half frame
%   cell_id     : a number representing the cell ID
% Outputs:
%   scr_block   : a vector representing the channel decoded but scrambled
%block
%   err_flag        : a boolean representing the CRC error flag
%   BCCH_BCH_msg    : a vector representing the decoded BCCH-BCH-message
%   SFN_4_LSB       : a vector representing the 4 least significant bits
%od the system frame number
%   HRF             : a bit representing the half frame bit
%   kSSB_MSB        : a bit representing the kSSB most significant bit

    % Rate recovery 7.1.5
    E = length(code_block); % Input sequence length: 7.1.4
    if E ~= 864
        return
    end
    K = 56;                 % N_payload + N_CRC: 7.1.3
    N = 512;
    rec_block = BCH.recoverRatePolar(code_block, K, N);
    
    % Polar decoding 7.1.4
    n_max = 9;
    iIL = true;
    n_PC = 0;               % Pairity check bit number: 7.1.4
    n_PC_wm = 0;
    
    L = 8;                  % list size to use during Successive
                            % Cancellation List (SCL) decoding
    N_CRC = 24;             % CRC length: 7.1.3
    poly_CRC = '24C';       % CRC polynomial 7.1.3
    pad_CRC = false;         % ??? default, for BCH and UCI
    rnti = 0;               % ??? default, unused
    dec_block = BCH.decodePolar(rec_block, K, E, n_max, n_PC, L, iIL, N_CRC);
    
    % CRC decoding 7.1.3
    [scr_block, err_flag] = utils.decodeCRC(dec_block, poly_CRC);
    
    % Descrambling 7.1.2
    tr_block = BCH.descramble(scr_block, cell_id, L_SSB);
    
    % Decode payload 7.1.1
    [BCCH_BCH_msg, SFN_4_LSB, HRF, kSSB_MSB] = decodePayload(tr_block, ...
        L_SSB);
end

function [BCCH_BCH_msg, SFN_4_LSB, HRF, kSSB_MSB] = decodePayload(...
    tr_block, L_SSB)
%DECODEPAYLOAD Decodes the PBCH payload according to 3GPP 38.212 7.1.1
% Inputs:
%   tr_block    : a vector representing the transport block
%   L_SSB       : a number representing the maximum number of SSBs in a
%   half frame
% Outputs:
%   BCCH_BCH_msg    : a vector represetning the decoded BCCH-BCH message
%   SFN_4_LSB       : a ventor representing the 4 least significant bits of
%   the system frame number
%   HRF             : a bit representing the half radio frame bit
%   kSSB_MSB        : a bit representing the most significant bit of the
%   kSSB

    A = length(tr_block);
    G = BCH.getPayloadInterleaverPattern();
    jSFN = 0;
    jHRF = 10;
    jSSB = 11;
    jOther = 14;
    aBar = zeros(A, 1);
    for i = 0:A-1
        if ((i>0 && i<=6) || (i>=24 && i<=27)) % SFN
            % Bits with idx 1...6 are the 6 MSBs of the SFN in the
            % MIB (3GPP 38.331 6.2.2). Bits with idx 24...27 are the 4 LSBs
            % of the SFN conveyed in additional PBCH payload bits (7.1.1)
            aBar(i+1, 1) = tr_block(G(jSFN+1)+1, 1);
            jSFN = jSFN+1;
        elseif i==28                    % HRF
            aBar(i+1, 1) = tr_block(G(jHRF+1)+1, 1);
        elseif i>=29                    % Candidate SS/PBCH block index
            aBar(i+1, 1) = tr_block(G(jSSB+1)+1, 1);
            jSSB = jSSB+1;
        else                            % Other 
            aBar(i+1, 1) = tr_block(G(jOther+1)+1, 1);
            jOther = jOther+1;
        end
    end

    BCCH_BCH_msg = aBar(1:24, 1);  % BCCH-BCH-Message (containing MIB)
    SFN_4_LSB = aBar(25:28, 1); % SFN 4 LSBs
    HRF = aBar(29, 1);    % HRF
    if L_SSB==64
        % MSB (6,5,4th) of SS Block index
        SSB_idx_3_MSB = aBar(30:32, 1);
    else
        % MSB of Kssb 
        kSSB_MSB = aBar(30, 1); 
    end
end