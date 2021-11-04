function [dci_bits, err_flag] = decodePDCCH(symbols, N_ID, N_RNTI, n_var)
%DECODEPDCCH Decodes the Physical Downlink Channel according to etsi 38.211
% and 38.212
% Inputs:
%   symbols     : a vector representing the modulatied PDCCH symbols
%   N_ID        : a number representing the pdcch-DMTS-ScramblingID
%   N_RNTI      : a number representing the PDCCH C-RNTI (default 0)
% Outputs:
%   dci_bits    : a vector representing the soft demodulated and
%   descrambled Downlink Control Information codeword
%   err_flag    : a number representing the CRC checksum

    % Decode PDCCH - etsi 38.211
    % 7.3.2.4 Demodulation
    modulation_order = 4;   % QPSK
    dci_cw = utils.demodulateQAM(symbols, modulation_order, n_var);
    
    % 7.3.2.3 Descrambling
    M = length(dci_cw);
    c_init = mod(N_RNTI*2^16+N_ID, 2^31);
    c = (-2)*utils.generatePRBS(M, c_init)+1;
    dci_cw = dci_cw .* c;

    
    % Decode DCI codeword
    E = length(dci_cw);
    A = 37;         % Payload size: 7.3.1.2.1 sum of all fields
    N_crc = 24;     % CRC size: 7.3.2
    poly_crc = '24C';
    K = A + N_crc;  % DCI size: 7.3.2
    n_max = 9;      % 7.3.3
    n_PC = 0;       % 7.3.3
    n_PC_wm = 0;    % 7.3.3
    iIL = true;     % 7.3.3
    N = channel.getNPolarEncoded(K,E,n_max);
    
    % 7.3.4 Rate recovering
    rec_block = channel.recoverRatePolar(dci_cw, K, N);

    % 7.3.3 Polar decoding
    L = 8;  % List size to use during Successive Cancellation List decoding
    dec_block = channel.decodePolar(rec_block, K, E, n_max, n_PC, L, iIL, N_crc);

    % 7.3.2 CRC decoding
    crc_block = [ones(24, 1); dec_block];
    si_RNTI = 65535; % 38.321 Table 7.1-1
    [dci_bits, err_flag] = utils.decodeCRC(crc_block, poly_crc, si_RNTI);
    dci_bits = dci_bits(25:end);
end