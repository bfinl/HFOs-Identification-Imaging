%% Time-frequency analysis
% Time window sieving
% This function combines events happening simultaneously from multiple
% channels for later extraction of multichannel epochs for the identified
% events.
%
%%% License
% We provide our code and data under a CC-BY-NC-SA-4.0 license, "As is" and
% without any guarantee to the scientific community for academic and
% research purposes primarily, not commercial use.
% 
% You should have received a copy of the CC-BY-NC-SA-4.0 license along with
% this program. If not, see https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [sig1,sig2,sig_idx] = jc_tfa_timingSieve(sig_idx,nWin,sig1,sig2)

% sample index information
idxEventChan = find(cellfun(@isempty,sig1)==0);
nEventChan = length(idxEventChan);

% find out event index and number, channel index
eventInfo = cell(nEventChan,1); % channel index, event sample index
for i = 1:nEventChan
    idxEvent_tmp = sig_idx{idxEventChan(i)};
    nEvent_tmp   = length(idxEvent_tmp);
    eventInfo{i} = zeros(nEvent_tmp,3);
    eventInfo{i}(:,1) = idxEventChan(i);
    eventInfo{i}(:,2) = 1:nEvent_tmp;
    eventInfo{i}(:,3) = idxEvent_tmp(:);
end
eventInfo = cat(1,eventInfo{:});

if size(eventInfo,1) < 2
    % no simultaneous events
    return
else
    % loop through all events if there exist multiple events
    %==========================================================================
    % find all possible groups according to how close these events are
    eventInfo_tmp = eventInfo;
    nEvent = size(eventInfo_tmp,1);
    eventGrp = cell(nEvent,1);
    nEvent_left = nEvent;
    j = 1;
    while nEvent_left > 0

        if nEvent_left == 1
            eventGrp{j} = eventInfo_tmp(1,:);
            eventInfo_tmp(1,:) = [];
        else        
            eventCurrentInfo = eventInfo_tmp(1,:);
            eventPoolInfo = eventInfo_tmp(2:end,:);
            idxEventMatch = find( eventPoolInfo(:,3) > (eventCurrentInfo(3)-nWin/2) &...
                eventPoolInfo(:,3) < (eventCurrentInfo(3)+nWin/2) );
            eventGrp{j} = [eventInfo_tmp(1,:); eventPoolInfo(idxEventMatch,:)];
            eventInfo_tmp([1;idxEventMatch+1],:) = [];
        end
        j = j+1;
        nEvent_left = size(eventInfo_tmp,1);

    end

    % for each group, only keep the event with maximum MGFP or power of
    % the filtered signal
    idxEventGrp = find(cellfun(@isempty,eventGrp)==0);
    nEventGrp = length(idxEventGrp);
    nSamp = size(sig2{idxEventChan(1)},2);
    eventGrp_discard_idx_set = [];
    for i = 1:nEventGrp
        eventGrp_tmp = eventGrp{i};
        nEvent_tmp = size(eventGrp_tmp,1);
        if nEvent_tmp > 1
            sigFilt_tmp = zeros(nEvent_tmp,nSamp);
            for j = 1:nEvent_tmp
                sigFilt_tmp(j,:) = sig2{eventGrp_tmp(j,1)}(eventGrp_tmp(j,2),:);
            end
            % get event index to keep and to discard
            sigFiltPwr = rms(sigFilt_tmp(:,floor(end/2-nWin/2):floor(end/2+nWin/2)),2);
            sigFiltPwr_max = max(sigFiltPwr);
            eventGrp_tmp_discard_idx = eventGrp_tmp(sigFiltPwr~=sigFiltPwr_max,:);
            % store discard event index
            eventGrp_discard_idx_set = [eventGrp_discard_idx_set; eventGrp_tmp_discard_idx]; %#ok<AGROW>
        end
    end
    if nEventGrp < nEvent
        % discard events
        eventChan_discard = unique(eventGrp_discard_idx_set(:,1));
        nChan_tmp = length(eventChan_discard);
        for j = 1:nChan_tmp
            eventIdx_discard = eventGrp_discard_idx_set(:,1)==eventChan_discard(j);
            eventIdx_discard = eventGrp_discard_idx_set(eventIdx_discard,2);

            sig1{eventChan_discard(j)}(eventIdx_discard,:) = [];
            sig2{eventChan_discard(j)}(eventIdx_discard,:) = [];
            sig_idx{eventChan_discard(j)}(eventIdx_discard,:) = [];
        end
    end
    %==========================================================================
end

end


