%% HFO detection
% This script is to conduct initial detection using a high sensititivity
% automatic detector, built based on amplitude and power information.
%
% The detector generally follows the following steps.
% (1) Raw data is exported and imported into workspace.
% (2) A 64-order FIR band-pass filter is used in 80-240 Hz range.
% (3) Moving SD is computed with 100 ms windows.
% (4) An baseline threshold (power) is set at 2x median SDs.
% (5) An epoch of 128 ms before and after the sample exceeding the
% threshold is extracted on the raw data.
% (6) HFO sieving procedure:
%   - (a) Raw data with zero-crossing over 10 times are rejected as
%       artifacts.
%   - (b) Signal envelope is computed using Hilbert ransform for each
%       extracted segment after filtering.
%   - (c) A threshold is defined at 2x median of the background envelope,
%       which is obtained from the first and last 80 ms segment.
%   - (d) Filtered data with threshold-crossing over 8 times, sitting at the
%       center are considered potential HFO events.
%   - (e) Raw data with zero-crossing over 10 would be rejected as noisy
%       activities.
%   - (f) Overlapped events are combined for later multichannel event
%       extraction.
% (7) The remaining HFO candidates are then pooled for later procedures.
% (8) The results are saved and a brief summary of the detected events is
%   saved in the file 3_HFOEvents_compact.mat.
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


%% zero-phase FIR filter: 80-240 Hz (64-order)
%%% signal parameters
nChan     = segmentInfo{3,1}.nChan;
nEpochLen = segmentInfo{3,1}.nSample;
fs        = segmentInfo{3,1}.Fs;
%%% band-pass filtering
[sigRawEpoch,nEpoch] = jc_io_segmentEpoch(sigRaw,nChan,nEpochLen);
sigFiltEpoch         = jc_tfa_filtFIREpoch(sigRawEpoch,fs,80,240,64,'band');

%% preallocate detection epoch and info. structure
%%% sigEpochRawSieve_set and sigEpochFiltSieve_set contain temporary events
%%% sigEpochTime_set and sigEpochSamp_set are timing of events
%%% sigEpochRawWhole_set and sigEpochFiltWhole_set are multichannel data
%%% sigEpochIndex_set contains index of spike epoch for each event
sigEpochRawSieve_set  = cell(nChan,1);
sigEpochFiltSieve_set = cell(nChan,1);
sigEpochTime_set      = cell(nChan,1);
sigEpochSamp_set      = cell(nChan,1);
sigEpochRawWhole_set  = cell(nChan,1);
sigEpochFiltWhole_set = cell(nChan,1);
sigEpochIndex_set     = cell(nChan,1);

% main loop for detection
% prepare channel information
chan_list = jc_io_readChannelList(pathChannel);
for iEpoch = 1:nEpoch

    % get one epoch at each loop
    sigRaw_tmp  = sigRawEpoch(:,:,iEpoch);
    sigFilt_tmp = sigFiltEpoch(:,:,iEpoch);
    
    %%% Baseline detector
    % compute moving std along each segment
    % moving window: 100 ms
    nWin = 100/1000*fs;
    sigFiltStd = movstd(sigFilt_tmp,nWin,0,2);
    % window size for combining events, if found multiple: 100 ms
    % window size for extraction: 128 ms before and after
    % baseline threshold: 2x median of the moving std
    nWin = 100/1000*fs;
    nWin_epoch = 128/1000*fs;
    thrd = 2*median(sigFiltStd,2);
    sigEpochIdx = jc_tfa_ampDetector(sigFilt_tmp,thrd,nWin);
    [sigEpochRaw,sigEpochIdx] = jc_io_extractEpoch(sigRaw_tmp,sigEpochIdx,nWin_epoch);
    sigEpochFilt              = jc_io_extractEpoch(sigFilt_tmp,sigEpochIdx,nWin_epoch);

    %%% HFO sieving procedure
    % (a) raw data zero-crossing below 10 times
    % suppose the activity-of-interest is sitting at the middle of each epoch
    % count zero-crossing times within the defined signal range
    % here, we just use the full length of the signal
    nWin = [];
    [sigEpochRawSieve_1,sigEpochFiltSieve_1,sigEpochIdx] = ...
        jc_tfa_zeroxSieve(sigEpochRaw,0,10,'<',nWin,sigEpochFilt,sigEpochIdx);

    % (b) filtered envelope threshold with counted crossing times
    % threshold is defined based on the filtered background activity
    % threshold: 2x median of the background envelope
    % window size for background signal: first and last 80 ms of the signal
    % window size for target event: 80 ms at the middle
    nWin = 80/1000*fs;
    [sigEpochFiltSieve,sigEpochRawSieve,sigEpochIdx] = ...
        jc_tfa_envelopeSieve(sigEpochFiltSieve_1,2,8,'>=',nWin,sigEpochRawSieve_1,sigEpochIdx);

    %%% HFO event timing sieving
    % HFO events are detected in multi-channel signals
    % thus, multiple events could exist at different channels within a
    % small time window, say 50 ms before and after the primary event
    % this procedure compares the timing of simultaneous events and sieve
    % out those smaller events, and only keep the one with the prominant MGFP
    % as the whole multichannel epoch would be extracted later which 
    % includes these simultaneous activities
    % or, as an alternative, this step can be simply commented out, then
    % every events from smultiple channels would be retained independently
    nWin = 100/1000*fs;
    [sigEpochRawSieve,sigEpochFiltSieve,sigEpochIdx] = ...
        jc_tfa_timingSieve(sigEpochIdx,nWin,sigEpochRawSieve,sigEpochFiltSieve);

    % print results
    nEvent = size(cell2mat(sigEpochRawSieve),1);
    nEventChan = length(find(cellfun(@isempty,sigEpochRawSieve)==0));
    jc_print_block(sprintf('Epoch %d: %d / %d event(s) / electrode(s).',iEpoch,nEvent,nEventChan));
    
    % log results
    % get channel information
    tmp_chanEvent_idx  = find(~cellfun(@isempty,sigEpochRawSieve));
    [tmp_iEvent_idx,~] = cellfun(@size,sigEpochRawSieve);
    tmp_iEvent_idx(tmp_iEvent_idx==0) = [];
    tmp_chanEvent_list = repelem(chan_list(tmp_chanEvent_idx),tmp_iEvent_idx);
    % event sample and time location
    % event sample index is the sample index when event occurrs
    % event time is the time when event occurrs
    sigEpochSamp = cell(nChan,1);
    sigEpochSamp(tmp_chanEvent_idx) = sigEpochIdx(tmp_chanEvent_idx);
    sigEpochTime = cellfun(@jc_bs_convertSamp2Time,sigEpochSamp,repelem({fs},nChan,1),'uniformoutput',false);
    % extract event epoch for whole channels
    sigEpochRawWhole  = jc_io_extractEpochWhole(sigRaw_tmp,sigEpochSamp,nWin_epoch);
    sigEpochFiltWhole = jc_io_extractEpochWhole(sigFilt_tmp,sigEpochSamp,nWin_epoch);

    %% construct event pool and event info.
    for iChan = tmp_chanEvent_idx'

        sigEpochRawSieve_set{iChan}     = [sigEpochRawSieve_set{iChan};
                                           sigEpochRawSieve{iChan}];
        sigEpochFiltSieve_set{iChan}    = [sigEpochFiltSieve_set{iChan};
                                           sigEpochFiltSieve{iChan}];
        sigEpochTime_set{iChan}         = [sigEpochTime_set{iChan};
                                           sigEpochTime{iChan}];
        sigEpochSamp_set{iChan}         = [sigEpochSamp_set{iChan};
                                           sigEpochSamp{iChan}];
        sigEpochRawWhole_set{iChan}     = cat(3,sigEpochRawWhole_set{iChan},...
                                           sigEpochRawWhole{iChan});
        sigEpochFiltWhole_set{iChan}    = cat(3,sigEpochFiltWhole_set{iChan},...
                                           sigEpochFiltWhole{iChan});
        sigEpochIndex_set{iChan}        = [sigEpochIndex_set{iChan};...
                                           repmat(iEpoch,length(sigEpochTime{iChan}),1)];
    end

end

% construct detected event pool and info. structure
% compose event pool
sigEventPoolRaw  = cell2mat(sigEpochRawSieve_set);
sigEventPoolFilt = cell2mat(sigEpochFiltSieve_set);
sigEventPoolRawWhole  = jc_bs_composeEventPoolWhole(sigEpochRawWhole_set);
sigEventPoolFiltWhole = jc_bs_composeEventPoolWhole(sigEpochFiltWhole_set);

% compose event info. structure
nEvent = size(sigEventPoolRaw,1);
sigEventInfo = struct('epochIdx',cell(nEvent,1),...
    'chanName',cell(nEvent,1),...
    'timeLoc',cell(nEvent,1),...
    'sampLoc',cell(nEvent,1),...
    'segIdx',cell(nEvent,1),...
    'sbjState',cell(nEvent,1),...
    'spkType',cell(nEvent,1));

% channel information
chan_list = jc_io_readChannelList(pathChannel);
chanEvent_idx  = find(~cellfun(@isempty,sigEpochRawSieve_set));
[iEvent_idx,~] = cellfun(@size,sigEpochRawSieve_set);
iEvent_idx(iEvent_idx==0) = [];
chanEvent_list = repelem(chan_list(chanEvent_idx),iEvent_idx);
[sigEventInfo.chanName] = deal(chanEvent_list{:});

% event timing and index
sigEventTime = cell2mat(sigEpochTime_set);
sigEventSamp = cell2mat(sigEpochSamp_set);
sigEventIndex = cell2mat(sigEpochIndex_set);
sigEventTime_tmp = num2cell(sigEventTime);
[sigEventInfo.timeLoc] = deal(sigEventTime_tmp{:});
sigEventSamp_tmp = num2cell(sigEventSamp);
[sigEventInfo.sampLoc] = deal(sigEventSamp_tmp{:});
sigEventIndex_tmp = num2cell(sigEventIndex);
[sigEventInfo.epochIdx] = deal(sigEventIndex_tmp{:});

% other event information
for iEvent = 1:nEvent
    iEpoch = sigEventInfo(iEvent).epochIdx;
    sigEventInfo(iEvent).segIdx   = sigEpochInfo(iEpoch).segNum;
    sigEventInfo(iEvent).sbjState = sigEpochInfo(iEpoch).sbjState;
    sigEventInfo(iEvent).spkType  = sigEpochInfo(iEpoch).spkType;
end

% print results
nEvent = size(sigEventPoolRaw,1);
nEventChan = length(find(cellfun(@isempty,sigEpochRawSieve_set)==0));
jc_print_single(sprintf('Total %d event(s) in %d electrode(s).',nEvent,nEventChan));

% save results for later labeling
jc_io_save([pathOutput,'\3_HFOEvents_compact'],...
    'sigEventPoolRaw','sigEventPoolFilt',...
    'sigEventPoolRawWhole','sigEventPoolFiltWhole',...
    'sigEventTime','chanEvent_list','sigEventIndex',...
    'sigEpochInfo','sigEventInfo','chan_list','fs');
% save all results
% save(fullfile(pathOutput,'2_HFODetection.mat'));
jc_print_block('HFO initial detection done.');


