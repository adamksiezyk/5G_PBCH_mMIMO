function [subcarriers, symbols] = getPDCCHResources(carrier, pdcch)
%GETPDCCHRESOURCES Summary of this function goes here
%   Detailed explanation goes here
    
    [is_occasion, candidates] = getPDCCHCandidates(carrier, pdcch);
end

function [is_occasion, candidates] = getPDCCHCandidates(carrier, pdcch)
%
    
    is_occasion = isOccasion(carrier, pdcch);
    if is_occasion
        % Compute candidates for all aggregation levels
        candidates = getAllCCEIndexes(pdcch.SearchSpace, pdcch.CORESET, ...
                                   double(pdcch.RNTI), slotNum);
    else
        % not monitored, bail out
        candidates = {zeros(0, 1, 'uint32')};
    end
end

function is_occasion = isOccasion(carrier, pdcch)
% Returns true if a PDCCH monitoring occasion exists in current slot
    slots
    ss = pdcch.SearchSpace;
    is_occasion = (mod(carrier.NFrame * slotsPerFrame + ...
        carrer.NSlot - ss.SlotOffset, ss.SlotPeriod) == 0);
end

function candidates = getAllCCEIndexes(search_space, coreset, rnti, N_slot)
% For all candidates, return the CCE indexes (L-sized).
% Slot based, not per monitored occasion in a slot (assume only one per SS)

    nCI = 0;   % Assumes nCI = 0 (carrier indicator field)
    numCCEs = double(coreset.Duration) * sum(coreset.FrequencyResources);

    Yp = nr5g.internal.pdcch.getYp(search_space.SearchSpaceType,search_space.CORESETID,...
            rnti,N_slot);

    % Determine candidates for each aggregation level
    aggLvls = [1 2 4 8 16];
    candidates = cell(5,1);
    for i = 1:5 % for AL {1,2,4,8,16}
        MsAL = search_space.NumCandidates(i);
        L = aggLvls(i);

        % Store column-wise CCEIndices for each candidate
        cceIdx = zeros(L,MsAL,'uint32');
        for ms = 0:MsAL-1
            for idx = 0:L-1
                cceIdx(idx+1,ms+1) = L*( mod(Yp + ...
                    floor(double(ms*numCCEs)/double(L*MsAL)) + nCI, ...
                    floor(numCCEs/L)) ) + idx;
            end
        end
        candidates{i} = cceIdx;
    end
end