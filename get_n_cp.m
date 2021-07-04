function [N_CP_long, N_CP, N_ECP_long, N_ECP] = get_n_cp(scs, sample_rate)
%GET_N_CP Calculates the length of cyclic prefix, long cyclic prefix
%and extended cyclic prefix in samples based on 3GPP TS 38.211 5.3.1
% Inputs:
%   scs             : a number representing the subcarrier spacing in Hz
%   sample_rate     : a number representing the sampling rate in Hz
% Outputs:
%   N_CP_long   : a number representing the long cyclic prefix length in
%   samples
%   N_CP        : a number representing the cyclic prefix length in samples
%   N_ECP_long  : (present if scs == 60e-3) a number representing the
%   long extending cyclic prefix length in samples
%   N_ECP       : (present if scs == 60e-3) a number representing the
%   extending cyclic prefix length in samples

    if (scs == 15e3)
        CP_len_long = 5.2e-6;
        CP_len = 4.69e-6;
    elseif (scs == 30e3)
        CP_len_long = 2.86e-6;
        CP_len = 2.34e-6;
    elseif (scs == 60e3)
        CP_len_long = 1.69e-6;
        CP_len = 1.17e-6;
        ECP_len_long = 4.17e-6;
        ECP_len = 4.17e-6;
        N_ECP_long = round(ECP_len_long * sample_rate);
        N_ECP = round(ECP_len * sample_rate);
    elseif (scs == 120e3)
        CP_len_long = 1.11e-6;
        CP_len = 0.59e-6;
    elseif (scs == 240e3)
        CP_len_long = 1.11e-6;
        CP_len = 0.29e-6;
    end
    N_CP_long = round(CP_len_long * sample_rate);
    N_CP = round(CP_len * sample_rate);
end

