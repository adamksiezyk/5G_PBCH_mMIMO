%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pss_generate(NID1)
% generating PSS sequence for given NID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq = pss_generate(NID1)
    %[38.211 7.4.2.2.1]
    x_table = zeros(1, 127);
    x_table(1, 1:7) = [0, 1, 1, 0, 1, 1, 1];
    for n = 8:127
        x_table(1, n) = mod((x_table(1, n-3) - x_table(1, n-7)), 2);
    end
    n = (0:126);

    seq = 1 - 2 * x_table(1, (mod((n+43*NID1), 127)+1));
end