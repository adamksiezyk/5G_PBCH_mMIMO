function output = recoverRatePolar(input, K, N, varargin)
%MATCHRATEPOLAR Returns the rate-recovered output of the input sequence
%INPUT of length E. The output OUTPUT is of size N x K and represents
%the information block. The process is the inverse of rate matching is
%described in 3GPP TS 38.212 5.4.1. Interleaving of coded bits is applied
%to uplink signals and turned on by the optional bil flag.
% Inputs:
%   code_block_     : a vector representing the code block
%   K               : a number representing the number of information bits
%   N               : a number representing the polar encoded block length
%   IBIL            : a boolean flag that turns on coded bits
%deinterleaving. The default value is flase.

    % Default coded bits deinterleaving flag
    if nargin == 3
        IBIL = false;
    else
        IBIL = varargin{1};
    end
    
    % Coded bits deinterleaving
    if IBIL
        code_block = codedBitsDeinterleaving(input);
    else
        code_block = input;
    end
    
    % Bit selection
    out_bs = bitSelection(code_block, K, N);
    
    % Sub-block deinterleaving
    output = subBlockDeinterleaving(out_bs);
end

function output = codedBitsDeinterleaving(input)
% Coded bits deinterleaving. 3GPP TS 38.212 5.4.1.3
% Inputs:
%   input   : a vector representing the input sequence
% Outputs:
%   output  : a vector representing the bit deinterleaved output sequence

    % Get T off E
    E = length(input);
    T = ceil((-1+sqrt(1+8*E))/2);

    % Create the table with nulls (filled in row-wise)
    vTab = zeros(T, T, class(input));
    k = 0;
    for i = 0:T-1
        for j = 0:T-1-i
            if k < E
                vTab(i+1, j+1) = k+1;
            end
            k = k+1;
        end
    end

    % Write input to buffer column-wise, respecting vTab
    v = Inf*ones(T, T, class(input));
    k = 0;
    for j = 0:T-1
        for i = 0:T-1-j
            if k < E && vTab(i+1, j+1) ~= 0
                v(i+1, j+1) = input(k+1);
                k = k+1;
            end
        end
    end

    % Read output from buffer row-wise
    output = zeros(size(input),class(input));
    k = 0;
    for i = 0:T-1
        for j = 0:T-1-i
            if ~isinf(v(i+1, j+1))
                output(k+1) = v(i+1, j+1);
                k = k+1;
            end
        end
    end
end

function output = bitSelection(input, K, N)
% Bit selection. 3GPP 38.212 5.4.1.2.
% Inputs:
%   imput   : a vector representing the input selquence
%   K       : a number representing the number of information bits
%   N       : a number representing the polar encoded block length
% Outputs:
%   output  : a vector representing the output sequence

    E = length(input);
    if E >= N
        output = input(1:N);
    else
        if K/E <= 7/16  % puncturing
            output = zeros(N, 1, class(input));
            output(end-E+1:end) = input;
        else            % shortening
            output = 1e20*ones(N, 1, class(input));
            output(1:E) = input;
        end
    end
end

function output = subBlockDeinterleaving(input)
% Sub-block deinterleaving. 3GPP TS 38.212 5.4.1.1
% Inputs:
%   input   : a vector representing the input sequence
% Outputs:
%   output  : a vector representing the sub-block deinterleaved

    % Sub-block interleaver pattern P(i): 3GPP TS 38.212 Table 5.4.1.1-1
    P = [0;1;2;4; 3;5;6;7; 8;16;9;17; 10;18;11;19;
          12;20;13;21; 14;22;15;23; 24;25;26;28; 27;29;30;31];

    N = length(input);
    jn = zeros(N,1);
    for n = 0:N-1
        i = floor(32*n/N);
        jn(n+1) = P(i+1)*(N/32)+ mod(n,N/32);
    end
    output = zeros(N, 1, class(input));
    output(jn+1) = input;
end