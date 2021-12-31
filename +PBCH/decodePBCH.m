function BCH_code_block = decodePBCH(PBCH_symbols, cell_id, v, n_var)
%DECODEPBCH Returns the SSB index based on the given Cell ID and PBCH
%DM-RS sequence
% Inputs:
%   PBCH_symbols    : a vector representing the PBCH symbols
%   cell_id         : a number representing the Cell ID
%   v               : a number representing the 3 LSBs of the SSB/beam ID
%   n_var           : a number representing the estimated noise variance
% Outputs:
%   BCH_code_block  : a vector representing the BCH code block
    
    % Demodulate 3GPP 38.211 7.3.3.2
    PBCH_demod = utils.demodulate(PBCH_symbols, 'QPSK', n_var);
    N_PBCH = length(PBCH_demod);
    
    % Generate the PRBS for the given Cell ID and all iSSBs
    c = utils.generatePRBS((v+1)*N_PBCH, cell_id);
    
    % Get the beam specific bipolar sequence
    c_descr = (-2)*c([1:N_PBCH]+v*N_PBCH)+1;
    
    % Descramble
    BCH_code_block = PBCH_demod .* c_descr;
end

