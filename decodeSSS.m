function NID1 = decodeSSS(SSS, NID2, show_plots_)
%DECODESSS Decode the input SSS
% Inputs:
%   SSS             : a vector representing the SSS
%   NID2            : a number representing the decoded NID2
%   show_plots_     : a boolean representing the show plots flag
% Outputs:
%   NID2    : a number representing the decoded NID1

    if nargin<3
        show_plots = true;
    else
        show_plots = show_plots_;
    end

    N_SSS = 336; % Amout of all possible SSS
    
    % Correlation between detected SSS and a reference SSS
    corrSSS = zeros(1, N_SSS);
    for n = 1:N_SSS
        sss_ref = generateSSS(n-1, NID2);
        corrSSS(n) = sqrt(abs(sum(SSS .* conj(sss_ref(1:end))) .^2));
    end
    
    [max_val, NID1] = max(corrSSS);
    NID1 = NID1-1;

    if show_plots
        figure;
        hold on;
        stem(0:335, corrSSS, 'o');
        title('Detected SSS sequence');
        xlabel('Sequence number N_{ID}^{(1)}');
        plot(NID1, max_val, 'kx', 'LineWidth', 2, 'MarkerSize', 8);
    end
end