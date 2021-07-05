function NID1 = sss_decode(SSS, NID2)
%SSS_DECODE Decode the input SSS
% Inputs:
%   SSS     : a vector representing the SSS
%   NID2    : a number representing the decoded NID2
% Outputs:
%   NID2    : a number representing the decoded NID1

    N_SSS = 336; % Amout of all possible SSS
    
    % Correlation between detected SSS and a reference SSS
    corrSSS = zeros(1, N_SSS);
    for n = 1:N_SSS
        sss_ref = sss_generate(n-1, NID2);
        corrSSS(n) = sqrt(abs(sum(SSS .* conj(sss_ref(1:end))) .^2));
    end

    figure;
    stem(0:335, corrSSS, 'o');
    title('Detected SSS sequence');
    xlabel('Sequence number NID1');
    
    [~, NID1] = max(corrSSS);
    NID1 = NID1-1;
end