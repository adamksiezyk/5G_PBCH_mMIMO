function [block, err] = decodeCRC(crc_block, polynomial, mask)
%DECODECRC Calculate and remove the CRC according to 3GPP TS 38.212 5.1
% Inputs:
%   crc_block   : a vector representing the CRC encoded block
%   polynomial  : a string representing the CRC polynomial type
%   mask        : a number representing the scrambling mask (optional)
% Outputs:
%   block   : a vector representing the decoded block
%   err     : a number representing the CRC decoding error

    if nargin < 3
        mask = 0;
    end

    gLen(:) = getPolyLength(polynomial);
    
    % Perform cyclic redundancy check for data (only part of CRC block)
    reEncodedBlk = utils.encodeCRC(crc_block(1:end-gLen,:), polynomial, mask);
    [codeLen, numCodeBlocks] = size(crc_block);
    if isempty(crc_block)
        block = zeros(0, numCodeBlocks, class(crc_block));
        err = zeros(0, numCodeBlocks, 'uint32');
    else
        block = reEncodedBlk(1:end-gLen,:);
        if codeLen <= gLen
            % For input length less than parity bit length
            blkcrcL = [false(gLen-codeLen,numCodeBlocks); crc_block>0];
            if mask
                maskBits = comm.internal.utilities.de2biBase2LeftMSB( ...
                    double(mask), gLen)';
                errBits = xor(blkcrcL,repmat(maskBits>0,[1 numCodeBlocks]));
            else
                errBits = blkcrcL;
            end
        else
            errBits = xor(reEncodedBlk(end-gLen+1:end,:)>0, ...
                crc_block(end-gLen+1:end,:));
        end
        err = uint32(sum(double(errBits).*repmat((2.^(gLen-1:-1:0)'), ...
            [1 numCodeBlocks])));
    end
end

function len = getPolyLength(poly)
    switch poly
        case '6'
            len = 6;
        case '11'
            len = 11;
        case '16'
            len = 16;
        case {'24a','24A'}
            len = 24;
        case {'24b','24B'}
            len = 24;
        otherwise % {'24c','24C'}
            len = 24;
    end
end