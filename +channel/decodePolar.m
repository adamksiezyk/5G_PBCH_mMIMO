function output = decodePolar(input, K, E, n_max, n_PC, L, iIL, N_CRC)
%DECODEPOLAR Decodes the polar encoded sequenve INPUT. This process is the
%inverse of polar coding described in 3GPP TS 38.212 5.3.1.
% Inputs:
%   input   :
% Outputs:
%   output  :

    % Mother code block length: 5.3.1
    N = channel.getNPolarEncoded(K, E, n_max);
    
    % Frozen and information bits pattern: 5.3.1.2
    [F, ~] = getFrozeAndInformationIndices(K, E, N, n_PC);
    
    % Decode polar
    output = SCLDecode(input, F, L, iIL, N_CRC);
end

function Q = getPolarSequence()
%GETPOLARSEQUENCE Returns the polar sequence defined in 3GPP TS 38.212
%Table 5.3.1.2-1.
% Outputs:
%   Q   : a vector representing the polar sequence

    Q =   [  0    518     94    214    364    414    819     966
             1     54    204    309    654    223    814     755
             2     83    298    188    659    663    439     859
             4     57    400    449    335    692    929     940
             8    521    608    217    480    835    490     830
            16    112    352    408    315    619    623     911
            32    135    325    609    221    472    671     871
             3     78    533    596    370    455    739     639
             5    289    155    551    613    796    916     888
            64    194    210    650    422    809    463     479
             9     85    305    229    425    714    843     946
             6    276    547    159    451    721    381     750
            17    522    300    420    614    837    497     969
            10     58    109    310    543    716    930     508
            18    168    184    541    235    864    821     861
           128    139    534    773    412    810    726     757
            12     99    537    610    343    606    961     970
            33     86    115    657    372    912    872     919
            65     60    167    333    775    722    492     875
            20    280    225    119    317    696    631     862
           256     89    326    600    222    377    729     758
            34    290    306    339    426    435    700     948
            24    529    772    218    453    817    443     977
            36    524    157    368    237    319    741     923
             7    196    656    652    559    621    845     972
           129    141    329    230    833    812    920     761
            66    101    110    391    804    484    382     877
           512    147    117    313    712    430    822     952
            11    176    212    450    834    838    851     495
            40    142    171    542    661    667    730     703
            68    530    776    334    808    488    498     935
           130    321    330    233    779    239    880     978
            19     31    226    555    617    378    742     883
            13    200    549    774    604    459    445     762
            48     90    538    175    433    622    471     503
            14    545    387    123    720    627    635     925
            72    292    308    658    816    437    932     878
           257    322    216    612    836    380    687     735
            21    532    416    341    347    818    903     993
           132    263    271    777    897    461    825     885
            35    149    279    220    243    496    500     939
           258    102    158    314    662    669    846     994
            26    105    337    424    454    679    745     980
           513    304    550    395    318    724    826     926
            80    296    672    673    675    841    732     764
            37    163    118    583    618    629    446     941
            25     92    332    355    898    351    962     967
            22     47    579    287    781    467    936     886
           136    267    540    183    376    438    475     831
           260    385    389    234    428    737    853     947
           264    546    173    125    665    251    867     507
            38    324    121    557    736    462    637     889
           514    208    553    660    567    442    907     984
            96    386    199    616    840    441    487     751
            67    150    784    342    625    469    695     942
            41    153    179    316    238    247    746     996
           144    165    228    241    359    683    828     971
            28    106    338    778    457    842    753     890
            69     55    312    563    399    738    854     509
            42    328    704    345    787    899    857     949
           516    536    390    452    591    670    504     973
            49    577    174    397    678    783    799    1000
            74    548    554    403    434    849    255     892
           272    113    581    207    677    820    964     950
           160    154    393    674    349    728    909     863
           520     79    283    558    245    928    719     759
           288    269    122    785    458    791    477    1008
           528    108    448    432    666    367    915     510
           192    578    353    357    620    901    638     979
           544    224    561    187    363    630    748     953
            70    166    203    236    127    685    944     763
            44    519     63    664    191    844    869     974
           131    552    340    624    782    633    491     954
            81    195    394    587    407    711    699     879
            50    270    527    780    436    253    754     981
            73    641    582    705    626    691    858     982
            15    523    556    126    571    824    478     927
           320    275    181    242    465    902    968     995
           133    580    295    565    681    686    383     765
            52    291    285    398    246    740    910     956
            23     59    232    346    707    850    815     887
           134    169    124    456    350    375    976     985
           384    560    205    358    599    444    870     997
            76    114    182    405    668    470    917     986
           137    277    643    303    790    483    727     943
            82    156    562    569    460    415    493     891
            56     87    286    244    249    485    873     998
            27    197    585    595    682    905    701     766
            97    116    299    189    573    795    931     511
            39    170    354    566    411    473    756     988
           259     61    211    676    803    634    860    1001
            84    531    401    361    789    744    499     951
           138    525    185    706    709    852    731    1002
           145    642    396    589    365    960    823     893
           261    281    344    215    440    865    922     975
            29    278    586    786    628    693    874     894
            43    526    645    647    689    797    918    1009
            98    177    593    348    374    906    502     955
           515    293    535    419    423    715    933    1004
            88    388    240    406    466    807    743    1010
           140     91    206    464    793    474    760     957
            30    584     95    680    250    636    881     983
           146    769    327    801    371    694    494     958
            71    198    564    362    481    254    702     987
           262    172    800    590    574    717    921    1012
           265    120    402    409    413    575    501     999
           161    201    356    570    603    913    876    1016
           576    336    307    788    366    798    847     767
            45     62    301    597    468    811    992     989
           100    282    417    572    655    379    447    1003
           640    143    213    219    900    697    733     990
            51    103    568    311    805    431    827    1005
           148    178    832    708    615    607    934     959
            46    294    588    598    684    489    882    1011
            75     93    186    601    710    866    937    1013
           266    644    646    651    429    723    963     895
           273    202    404    421    794    486    747    1006
           517    592    227    792    252    908    505    1014
           104    323    896    802    373    718    855    1017
           162    392    594    611    605    813    924    1018
            53    297    418    602    848    476    734     991
           193    770    302    410    690    856    829    1020
           152    107    649    231    713    839    965    1007
            77    180    771    688    632    725    938    1015
           164    151    360    653    482    698    884    1019
           768    209    539    248    806    914    506    1021
           268    284    111    369    427    752    749    1022
           274    648    331    190    904    868    945    1023 ];
       
    Q = Q(:);
end

function [F, I] = getFrozeAndInformationIndices(K, E, N, n_PC)
%GETFROZENANDINFORMATIONINDICES Returns the frozen and information bit
%indices based on 3GPP TS 38.212 5.3.1.2 and 5.4.1.1

    % The polar sequence Q_N: 5.3.1.2
    Q_0_Nmax = getPolarSequence();
    Q_0_N = Q_0_Nmax(Q_0_Nmax<N);

    % Frozen and information bits indices
    % Subblock interleave map: 5.4.1.1
    P = [0 1 2 4 3 5 6 7 8 16 9 17 10 18 11 19 12 20 13 21 14 22 15 23 ...
         24 25 26 28 27 29 30 31];

    J = zeros(1,N);
    for n=0:N-1
        i=floor(32*n/N);
        J(n+1) = P(i+1)*(N/32)+mod(n,N/32);
    end

    Q_Ftmp_N = [];
    if E < N
        if K/E <= 7/16 % puncturing
            for n=0:N-E-1
                Q_Ftmp_N = [Q_Ftmp_N, J(n+1)];
            end
            if E >= 3*N/4
                Q_Ftmp_N = [Q_Ftmp_N,0:ceil(3*N/4-E/2)-1];
            else
                Q_Ftmp_N = [Q_Ftmp_N,0:ceil(9*N/16-E/4)-1];
            end
            Q_Ftmp_N = unique(Q_Ftmp_N);
        else % shortening
            for n=E:N-1
                Q_Ftmp_N = [Q_Ftmp_N, J(n+1)];
            end
        end
    end
    
    % Q_I_N sequence
    Q_I_N = zeros(K+n_PC, 1);
    j = 0;
    for i = 1:N
        ind = Q_0_N(N-i+1);     % Flip for most reliable
        if any(ind==Q_Ftmp_N)
            continue;
        end
        j = j+1;
        Q_I_N(j) = ind;
        if j==(K+n_PC)
            break;
        end
    end
    
    % Q_F_N sequence
    Q_F_N = setdiff(Q_0_N, Q_I_N);
    
    % Frozen and information bits indices
    F = zeros(N, 1);
    F(Q_F_N+1) = 1;
    I = zeros(N, 1);
    I(Q_I_N) = 1;
end

function output = SCLDecode(input, F, L, iIL, N_CRC)
%SCLDECODE A Successive Cancellation List Decoder for polar code decoding
%based on "List decoding of Polar Codes" by Tal and Vardy and "LLR-Based
%Successive Cancellation List Decoding of Polar Codes" by Stimming, Parizi
%and Burg.
% Inputs:
%   input   : a vector representing the soft demodulated symbols sequence
%   F       : a vector representing the frozen bits indices
%   L       : a number representing the number of decoding lists
%   iIL     : 
%   N_CRC   : a number representing the CRC length
%   padCRC  : a boolean, if true the CRC is added
%   RNTI    : a number representing the Radio Network Temporary Identifier.
%You can use this syntax when the value of rnti masks the CRC parity bits at the transmit end
% Outputs:
%   output  :

    % Setup
    N = length(F);
    m = log2(N);
    K = sum(F==0);      % includes nPC bits as well

    % CRCs as per TS 38.212, Section 5.1
    if N_CRC == 24         % '24C', downlink
        polyStr = '24C';
    elseif N_CRC == 11     % '11', uplink
        polyStr = '11';
    else % crcLen == 6      % '6', uplink
        polyStr = '6';
    end

    br = zeros(N,1);
    for idxBit = 0:N-1
        % 0-based indexing
        br(idxBit+1) = bitReverse(idxBit,m);
    end

    if iIL
        piInterl = nr5g.internal.polar.interleaveMap(K);
    else
        piInterl = (0:K-1).';
    end

    % Initialize core
    [sttStr, arrayPtrLLR, arrayPtrC] = initializeDataStructures(N,L,m);
    [iniPathIdx, sttStr] = assignInitialPath(sttStr);
    [sp, sttStr, arrayPtrLLR, arrayPtrC] = getArrayPtrP(sttStr, ...
        arrayPtrLLR, arrayPtrC, 1, iniPathIdx);
    arrayPtrLLR{1,sp}(:,1) = input(br+1);  % LLRs
    mplus1 = m+1;

    % Main loop
    for phase = 1:N
        [sttStr, arrayPtrLLR, arrayPtrC] = recursivelyCalcP(sttStr, ...
            arrayPtrLLR, arrayPtrC, mplus1, phase);

        pm2 = mod(phase-1,2);
        if F(phase)==1
            % Path for frozen (and punctured) bits
            for pathIdx = 1:L
                if ~sttStr.activePath(pathIdx)
                    continue;
                end
                [sc, sttStr, arrayPtrLLR, arrayPtrC] = getArrayPtrC(sttStr, ...
                    arrayPtrLLR, arrayPtrC, mplus1, pathIdx);
                arrayPtrC{mplus1,sc}(1,pm2+1) = 0; % set to 0

                tmp = arrayPtrLLR{mplus1,sc}(1,1);
                if tmp < 0
                    sttStr.llrPathMetric(pathIdx) = ...
                        sttStr.llrPathMetric(pathIdx) + abs(tmp);
                    % Else branch doesnt need an update
                end
            end
        else % Path for info bits
            [sttStr, arrayPtrLLR, arrayPtrC] = continuePathsUnfrozenBit(sttStr, ...
                arrayPtrLLR, arrayPtrC, phase);
        end

        if pm2==1
            [sttStr, arrayPtrLLR, arrayPtrC] = recursivelyUpdateC(sttStr, ...
                arrayPtrLLR, arrayPtrC, mplus1, phase);
        end
    end

    % Return the best codeword in the list. Use CRC checks, if enabled
    pathIdx1 = 1;
    p1 = realmax;
    crcCW = false;
    for pathIdx = 1:L
        if ~sttStr.activePath(pathIdx)
            continue;
        end

        [sc, sttStr, arrayPtrLLR, arrayPtrC] = getArrayPtrC(sttStr, ...
            arrayPtrLLR, arrayPtrC, mplus1, pathIdx);
        if N_CRC>0
            canCW = sttStr.savedCWs(:,sc);  % N, with frozen bits
            canMsg = canCW(F==0,1);         % K bits only (with nPC)
            canMsg(piInterl+1) = canMsg;    % deinterleave (for k+nPC)

            out = canMsg;
            
            % Check CRC: errFlag is 1 for error, 0 for no errors
            [~, errFlag] = utils.decodeCRC(out, polyStr);
            if errFlag      % ~0 => fail
                continue;   % move to next path
            end
        end
        crcCW = true;
        if p1 > sttStr.llrPathMetric(pathIdx)
            p1 = sttStr.llrPathMetric(pathIdx);
            pathIdx1 = pathIdx;
        end
    end

    if ~crcCW   % no codeword found which passes crcCheck
        pathIdx1 = 1;
        p1 = realmax;
        for pathIdx = 1:L
            if ~sttStr.activePath(pathIdx)
                continue;
            end

            if p1 > sttStr.llrPathMetric(pathIdx)
                p1 = sttStr.llrPathMetric(pathIdx);
                pathIdx1 = pathIdx;
            end
        end
    end

    % Get decoded bits
    [sc, sttStr] = getArrayPtrC(sttStr,arrayPtrLLR,arrayPtrC,mplus1, ...
                                pathIdx1);
    decCW = sttStr.savedCWs(:,sc);  % N, with frozen bits
    output = decCW(F==0,1);         % K, info + nPC bits only
    output(piInterl+1) = output;    % Deinterleave output, K+nPC
    output = logical(output);
end

function [struct, arrayPtrLLR, arrayPtrC] = initializeDataStructures(N,L,m)
% Algorithm 8: Initialize helper data structures

    struct.m = m;
    struct.L = L;

    parrayPtrLLR = cell(m+1,L);  % store arrays
    parrayPtrC = cell(m+1,L);    % store arrays
    coder.unroll(false);            % Allows maxMplus1, m+1 to differ
    for layer = 1:m+1
        expm = 2^(m+1-layer);
        coder.unroll(false);
        for s = 1:L
            parrayPtrLLR{layer,s} = zeros(expm,1);
            parrayPtrC{layer,s} = zeros(expm,2); % binary-valued: 0,1
        end
    end
    % An extra layer of in-direction is needed for codegen
    arrayPtrLLR = parrayPtrLLR;           % (m+1)-by-L
    arrayPtrC = parrayPtrC;               % (m+1)-by-L

    struct.llrPathMetric = zeros(L,1);

    struct.pathIdxToArrayIdx = ones(m+1,L);   % (m+1)-by-L

    struct.inactiveArrayIndices = zeros(m+1,L);
    struct.inactiveArrayIndicesLen = zeros(m+1,1);
    for layer = 1:m+1
        struct.inactiveArrayIndices(layer,:) = 1:L;
        struct.inactiveArrayIndicesLen(layer,1) = L;
    end
    struct.arrayReferenceCount = zeros(m+1,L);

    struct.inactivePathIndices = (1:L).';     % all paths are inactive
    struct.inactivePathIndicesLen = L;        % 1-by-1, stack depth
    struct.activePath = zeros(L,1,'logical'); % no paths are active

    struct.savedCWs = zeros(N,L);             % saved codewords

end

function [pathIdx, struct] = assignInitialPath(struct)
% Algorithm 9: Assigns initial path and returns the path index

    pathIdx = struct.inactivePathIndices(struct.inactivePathIndicesLen, 1);
    struct.inactivePathIndicesLen = struct.inactivePathIndicesLen-1;
    struct.activePath(pathIdx) = true;

    % Associate arrays with path index
    for layer = 1:struct.m+1
        s = struct.inactiveArrayIndices(layer, ...
            struct.inactiveArrayIndicesLen(layer));
        struct.inactiveArrayIndicesLen(layer) = ...
            struct.inactiveArrayIndicesLen(layer)-1;

        struct.pathIdxToArrayIdx(layer,pathIdx) = s;
        struct.arrayReferenceCount(layer,pathIdx) = 1;
    end
end

function [clPathIdx, struct] = clonePath(struct, pathIdx)
% Algorithm 10: Clones the given path
% Inpust:
%   pathIdx     : path index to clone
% Outputs:
%   clPathIdx   : cloned path index

    clPathIdx = struct.inactivePathIndices(struct.inactivePathIndicesLen,1);
    struct.inactivePathIndicesLen = struct.inactivePathIndicesLen-1;
    struct.activePath(clPathIdx) = true;
    struct.llrPathMetric(clPathIdx) = struct.llrPathMetric(pathIdx);

    % Make clPathIdx (l') reference same arrays as pathIdx (l)
    for layer = 1:struct.m+1
        s = struct.pathIdxToArrayIdx(layer,pathIdx);

        struct.pathIdxToArrayIdx(layer,clPathIdx) = s;
        struct.arrayReferenceCount(layer,s) = ...
            struct.arrayReferenceCount(layer,s)+1;
    end
end

function struct = killPath(struct, pathIdx)
% Algorithm 11: Kills the given path
% Inputs:
%   pathIdx     : path index to kill

    % Mark path pathIdx as inactive
    struct.activePath(pathIdx) = false;
    struct.inactivePathIndicesLen = struct.inactivePathIndicesLen+1;
    struct.inactivePathIndices(struct.inactivePathIndicesLen,1) = pathIdx;
    struct.llrPathMetric(pathIdx) = 0;

    % Disassociate arrays with path Idx
    for layer = 1:struct.m+1
        s = struct.pathIdxToArrayIdx(layer,pathIdx);
        struct.arrayReferenceCount(layer,s) = ...
            struct.arrayReferenceCount(layer,s)-1;

        if struct.arrayReferenceCount(layer,s)==0
            if struct.inactiveArrayIndicesLen(layer,1) < struct.L
                struct.inactiveArrayIndicesLen(layer,1) = ...
                    struct.inactiveArrayIndicesLen(layer,1)+1;
            end
            struct.inactiveArrayIndices(layer, ...
                struct.inactiveArrayIndicesLen(layer,1)) = s;
        end
    end
end

function [s2, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrP(struct, ...
    arrayPtrLLR, arrayPtrC, layer, pathIdx)
% Algorithm 12: Returns the pointer to corresponding probability array
% Inputs:
%   layer       : layer lambda
%   pathIdx     : path index
% Outputs:
%   s2  : corresponding pathIdx for same layer

    s = struct.pathIdxToArrayIdx(layer,pathIdx);
    if struct.arrayReferenceCount(layer,s)==1
        s2 = s;
    else
        s2 = struct.inactiveArrayIndices(layer, ...
            struct.inactiveArrayIndicesLen(layer,1));
        if struct.inactiveArrayIndicesLen(layer,1) > 1
            struct.inactiveArrayIndicesLen(layer,1) = ...
                struct.inactiveArrayIndicesLen(layer,1)-1;
        end

        % deep copy
        arrayPtrLLR{layer,s2} = arrayPtrLLR{layer,s};
        arrayPtrC{layer,s2} = arrayPtrC{layer,s};

        struct.arrayReferenceCount(layer,s) = ...
            struct.arrayReferenceCount(layer,s)-1;
        struct.arrayReferenceCount(layer,s2) = 1;
        struct.pathIdxToArrayIdx(layer,pathIdx) = s2;
    end
end

function [s2, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrC(struct, ...
    arrayPtrLLR, arrayPtrC, layer, pathIdx)
% Algorithm 12:  Returns the pointer to corresponding probability array
% Inputs:
%   layer       : layer lambda
%   pathIdx     : path index l
% Outputs:
%   s2  : corresponding pathIdx for same layer

    s = struct.pathIdxToArrayIdx(layer,pathIdx);
    if struct.arrayReferenceCount(layer,s)==1
        s2 = s;
    else
        s2 = struct.inactiveArrayIndices(layer, ...
            struct.inactiveArrayIndicesLen(layer,1));
        if struct.inactiveArrayIndicesLen(layer,1) > 1
            struct.inactiveArrayIndicesLen(layer,1) = ...
                struct.inactiveArrayIndicesLen(layer,1)-1;
        end

        % deep copy
        arrayPtrC{layer,s2} = arrayPtrC{layer,s};
        arrayPtrLLR{layer,s2} = arrayPtrLLR{layer,s};

        struct.arrayReferenceCount(layer,s) = ...
            struct.arrayReferenceCount(layer,s)-1;
        struct.arrayReferenceCount(layer,s2) = 1;
        struct.pathIdxToArrayIdx(layer,pathIdx) = s2;
    end
end

function br = bitReverse(b, n)
%BITREVERSE Bit-wise reverse input value

    br = comm.internal.utilities.bi2deRightMSB( ... % for fliplr
         comm.internal.utilities.de2biBase2LeftMSB(b,n),2);
end

function [struct, arrayPtrLLR, arrayPtrC] = recursivelyCalcP(struct, ...
    arrayPtrLLR, arrayPtrC, layer, phase)
% Algorithm 3
% Input:
%   layer: layer lambda
%   phase: phase phi

    if layer==1
        return;
    end
    psi = floor((phase-1)/2)+1;
    pm2 = mod(phase-1,2);
    if pm2==0
        [struct, arrayPtrLLR, arrayPtrC] = recursivelyCalcP(struct, ...
            arrayPtrLLR, arrayPtrC, layer-1, psi);
    end

    expm = 2^(struct.m-layer+1);
    for pathIdx = 1:struct.L
        if ~struct.activePath(pathIdx)
            continue;
        end

        [sp, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrP(struct, ...
            arrayPtrLLR, arrayPtrC, layer, pathIdx);
        [spminus1, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrP( ...
            struct, arrayPtrLLR, arrayPtrC, layer-1, pathIdx);
        [sc, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrC(struct, ...
            arrayPtrLLR, arrayPtrC, layer, pathIdx);
        for beta = 0:expm-1
            % LLR
            aa = arrayPtrLLR{layer-1,spminus1}( 2*beta+1,1 );
            bb = arrayPtrLLR{layer-1,spminus1}( 2*beta+2,1 );
            if pm2==0
                % Equation 4, Stimming
                arrayPtrLLR{layer,sp}(beta+1,1) = ...
                    sign(aa)*sign(bb)*min(abs(aa),abs(bb));
            else
                u1 = arrayPtrC{layer,sc}(beta+1,1);
                % Equation 5, Stimming
                arrayPtrLLR{layer,sp}(beta+1,1) = (-1)^u1 * aa + bb;
            end
        end
    end
end

function [struct, arrayPtrLLR, arrayPtrC] = recursivelyUpdateC(struct, ...
    arrayPtrLLR, arrayPtrC, layer, phase)
% Algorithm 4
% Input:
%   layer: layer lambda
%   phase: phase phi, has to be odd

    psi = floor((phase-1)/2);
    pm2 = mod(psi,2);
    expm = 2^(struct.m-layer+1);
    for pathIdx = 1:struct.L
        if ~struct.activePath(pathIdx)
            continue;
        end
        [sc, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrC(struct, ...
            arrayPtrLLR, arrayPtrC, layer, pathIdx);
        [scminus1, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrC( ...
            struct, arrayPtrLLR, arrayPtrC, layer-1, pathIdx);
        for beta = 0:expm-1
            arrayPtrC{layer-1,scminus1}(2*beta+1,pm2+1) = ...
                xor(arrayPtrC{layer,sc}(beta+1,1), ...
                arrayPtrC{layer,sc}(beta+1,2));
            arrayPtrC{layer-1,scminus1}(2*beta+2,pm2+1) = ...
                arrayPtrC{layer,sc}(beta+1,2);
        end
    end

    if pm2==1
        [struct, arrayPtrLLR, arrayPtrC] = recursivelyUpdateC(struct, ...
            arrayPtrLLR, arrayPtrC, layer-1, psi+1);
    end
end

function [struct, arrayPtrLLR, arrayPtrC] = continuePathsUnfrozenBit(...
    struct, arrayPtrLLR, arrayPtrC, phase)
% Helper function. Continue path for unfrozen bits
% Input:
%   phase: phase phi

    % Populate probForks
    probForks = -realmax*ones(struct.L,2);
    i = 0;
    mplus1 = struct.m+1;
    for pathIdx = 1:struct.L
        if struct.activePath(pathIdx)
            [sp, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrP(struct, ...
                arrayPtrLLR, arrayPtrC, mplus1, pathIdx);

            tmp = arrayPtrLLR{mplus1,sp}(1,1);
            if tmp > 0
                probForks(pathIdx,1) = - struct.llrPathMetric(pathIdx);
                probForks(pathIdx,2) = -(struct.llrPathMetric(pathIdx) ...
                    + tmp);
            else
                probForks(pathIdx,1) = - (struct.llrPathMetric(pathIdx) ...
                    + abs(tmp));
                probForks(pathIdx,2) = - struct.llrPathMetric(pathIdx);
            end
            i = i+1;
        end
    end

    rho = min(2*i,struct.L);
    contForks = zeros(struct.L,2);
    % Populate contForks such that contForks(l,b) is true iff
    % probForks(l,b) is one of rho largest entries in probForks.
    prob = sort(probForks(:), 'descend');
    if rho>0
        threshold = prob(rho);
    else
        threshold = prob(1); % Largest
    end
    numPop = 0;
    for pathIdx = 1:struct.L
        for bIdx = 1:2
            if numPop==rho
                break;
            end
            if probForks(pathIdx,bIdx)>threshold
                contForks(pathIdx,bIdx) = 1;
                numPop = numPop+1;
            end
        end
    end

    if numPop<rho
        for pathIdx = 1:struct.L
            for bIdx = 1:2
                if numPop==rho
                    break;
                end
                if probForks(pathIdx,bIdx)==threshold
                    contForks(pathIdx,bIdx) = 1;
                    numPop = numPop+1;
                end
            end
        end
    end

    % First, kill-off non-continuing paths
    for pathIdx = 1:struct.L
        if ~struct.activePath(pathIdx)
            continue;
        end
        if contForks(pathIdx,1)==0 && contForks(pathIdx,2)==0
            struct = killPath(struct, pathIdx);
        end
    end

    % Continue relevant paths, duplicating if necessary
    pm2 = mod(phase-1,2);
    for pathIdx = 1:struct.L
        if contForks(pathIdx,1)==0 && contForks(pathIdx,2)==0
            % Both forks are bad
            continue;
        end

        [sc, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrC(struct, ...
            arrayPtrLLR, arrayPtrC, mplus1, pathIdx);
        if contForks(pathIdx,1)==1 && contForks(pathIdx,2)==1
            % Both forks are good
            arrayPtrC{mplus1,sc}(1,pm2+1) = 0;
            struct.savedCWs(phase,sc) = 0;

            [pathIdx1, struct] = clonePath(struct, pathIdx);
            [sc2, struct, arrayPtrLLR, arrayPtrC] = getArrayPtrC(struct, ...
                arrayPtrLLR, arrayPtrC, mplus1, pathIdx1);
            struct.savedCWs(1:phase-1,sc2) = struct.savedCWs(1:phase-1,sc);

            arrayPtrC{mplus1,sc2}(1,pm2+1) = 1;
            struct.savedCWs(phase,sc2) = 1;

            tmp = arrayPtrLLR{mplus1,sc}(1,1);
            if tmp < 0
                struct.llrPathMetric(pathIdx) = ...
                    struct.llrPathMetric(pathIdx) + abs(tmp);
                % Else branch doesnt need an update
            end

            tmp2 = arrayPtrLLR{mplus1,sc2}(1,1);
            if tmp2 > 0
                struct.llrPathMetric(pathIdx1) = ...
                    struct.llrPathMetric(pathIdx1) + tmp2;
                % Else branch doesnt need an update
            end
        else
            % Exactly one fork is good
            tmp = arrayPtrLLR{mplus1,sc}(1,1);
            if contForks(pathIdx,1)==1
                arrayPtrC{mplus1,sc}(1,pm2+1) = 0;
                struct.savedCWs(phase,sc) = 0;

                if tmp < 0
                    struct.llrPathMetric(pathIdx) = ...
                        struct.llrPathMetric(pathIdx) + abs(tmp);
                    % Else branch doesnt need an update
                end
            else
                arrayPtrC{mplus1,sc}(1,pm2+1) = 1;
                struct.savedCWs(phase,sc) = 1;

                if tmp > 0
                    struct.llrPathMetric(pathIdx) = ...
                        struct.llrPathMetric(pathIdx) + tmp;
                    % Else branch doesnt need an update
                end
            end
        end
    end
end