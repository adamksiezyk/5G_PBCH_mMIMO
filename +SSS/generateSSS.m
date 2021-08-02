function seq = generateSSS(NID1, NID2)
%GENERATESSS Generates the SSS squence besed on the provided NID1 and NID2
%and 3GPP TS 38.211 7.4.2.3.1
    x0_table = zeros(1, 127);
    x0_table(1:7) = [1 0 0 0 0 0 0];
    x1_table = zeros(1, 127);
    x1_table(1:7) = [1 0 0 0 0 0 0];
    for n = 8:127
        x0_table(n) = mod((x0_table(n-3) + x0_table(n-7)), 2);
        x1_table(n) = mod((x1_table(n-6) + x1_table(n-7)), 2);
    end
    m0 = 15*fix(NID1/112) + 5*NID2;
    m1 = mod(NID1, 112);
    n = (0:126);
    seq = (1 - 2*x0_table(mod(n+m0, 127)+1)) .* (1 - 2*x1_table(mod(n+m1, 127)+1));
end