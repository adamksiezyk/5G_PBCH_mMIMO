function blkcrc = encodeCRC(blk, poly)
%ENCODECRC Add CRC to block
% Inputs:
%   blk     : a vector representing the code block
%   poly    : a string representing the CRC polynomial
% Outputs:
%   blkcrc  : a vector representing the CRC encoded block

    % Perform cyclic redundancy check
    [codeLen, numCodeBlocks] = size(blk);
    blkL = logical(blk);
    if isempty(blk)
        blkcrc = zeros([0,numCodeBlocks],'like',blk);
    else
        [gLen, gPoly] = getPoly(poly);
        blkcrcL = false(codeLen+gLen,numCodeBlocks);
        for i = 1:numCodeBlocks
            blkcrcL(:,i) = crcEncode(double(blkL(:,i)),gPoly,gLen);
        end
        blkcrc = [blk; cast(blkcrcL(end-gLen+1:end,:),class(blk))];
    end
end

function out = crcEncode(in,gPoly,gLen)
% CRC Encode with all doubles inputs, logical out.

    % Append zeros to the data block
    inPad = [in; zeros(gLen,1)];

    % Perform cyclic redundancy check
    remBits = [0; inPad(1:gLen,1)];
    for i = 1:length(inPad)-gLen
        dividendBlk = [remBits(2:end); inPad(i+gLen)];
        if dividendBlk(1) == 1
            remBits = rem(gPoly+dividendBlk,2);
        else
            remBits = dividendBlk;
        end
    end
    parityBits = remBits(2:end);

    out = logical([in; parityBits]);

end

function [len, gPoly] = getPoly(poly)
% Initialize CRC polynomial (gPoly)
    switch poly
        case '6'
            gPoly = [1 1 0 0 0 0 1]';
            len = 6;
        case '11'
            gPoly = [1 1 1 0 0 0 1 0 0 0 0 1]';
            len = 11;
        case '16'
            gPoly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1]';
            len = 16;
        case {'24a','24A'}
            gPoly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1]';
            len = 24;
        case {'24b','24B'}
            gPoly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1]';
            len = 24;
        otherwise % {'24c','24C'}
            gPoly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1]';
            len = 24;
    end
end