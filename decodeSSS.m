function NID1 = decodeSSS(SSS, NID2)
%DECODESSS Decode the input SSS
% Inputs:
%   SSS     : a vector representing the SSS
%   NID2    : a number representing the decoded NID2
% Outputs:
%   NID2    : a number representing the decoded NID1

    N_SSS = 336; % Amout of all possible SSS
    
    % Correlation between detected SSS and a reference SSS
    corrSSS = zeros(1, N_SSS);
    for n = 1:N_SSS
        sss_ref = generateSSS(n-1, NID2);
        corrSSS(n) = sqrt(abs(sum(SSS .* conj(sss_ref(1:end))) .^2));
    end

    figure;
    stem(0:335, corrSSS, 'o');
    title('Detected SSS sequence');
    xlabel('Sequence number N_{ID}^{(1)}');
    
    [max_val, NID1] = max(corrSSS);
    NID1 = NID1-1;
    
    hold on;
    plot(NID1, max_val, 'kx', 'LineWidth', 2, 'MarkerSize', 8);
end