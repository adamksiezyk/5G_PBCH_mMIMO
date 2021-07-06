%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_pr_sequence(c_init, seq_len)
%   generate pseudo-random sequence from 38.211 5.2.1
% Input:
%   c_init - initial value for sequence generation
%   seq_len - length of generated sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c_seq = generate_pr_sequence(c_init, seq_len)
    %pseudo-random sequence [38.211 5.2.1]
    Nc = 1600;
    x1_table = zeros(1, seq_len+Nc);
    x2_table = x1_table;
    x1_table(1) = 1;
    for n = 30:-1:0
        x2_table(n+1) = fix(c_init/2^n);
        c_init = mod(c_init, 2^n);
    end
    for n = 32:seq_len+Nc
        x1_table(n) = mod(x1_table(n-28) + x1_table(n-31), 2);
        x2_table(n) = mod(x2_table(n-28) + x2_table(n-29) + x2_table(n-30) + x2_table(n-31), 2);
    end
    c_seq = zeros(1, seq_len);
    for k = 1:seq_len
        c_seq(k) = mod(x1_table(k+Nc) + x2_table(k+Nc), 2);
    end
end