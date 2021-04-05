%% High-frequency Oscillation Simulation - Evaluation
% This function calculates several metrics to evaluate the results of
% auto-labeled HFO results, comparing to the simulated ground-truth.
% The simulation ground-truth is saved in the structure strAct, with the
% temporal and spatial (channel) information where the events, both true
% and false HFOs, are inserted in the simulated EEG epochs.
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


%% select optimal cluster of events
% use temporal homogeneity of events winin each cluster as reference
% This measure describes the uniformity of HFOs with regard to spikes as a
% reference for selecting optimal clusters.
% It is noteworthy that the cluster with HFOs riding on spikes is generally
% clear to identify by plotted spatiotemporal patterns. This measure is
% only a reference for selecting such clusters. This can also be integrated
% for automated selection as well, but more evaluation is recommended for
% the robustness of such measurement and implementations. Please verify and
% confirm the identified cluster following the procedure.
[ClusterSPLAvg,ClusterSPL] = computeClsUniMeasure(sigEventTime,clusterX,segmentInfo{3});

% visualize cluster selection criteria
fig = figure;
subplot(121), hold all;
plot(ClusterSPLAvg,'o-');
xlabel('Cluster Index');
ylabel('Within-cluster homogeneity');
title('Identify pHFOs cluster');
set(gca,'xtick',1:nCls,'fontsize',12);
subplot(122), hold all;
gscatter3(feature_pc(:,1),feature_pc(:,2),feature_pc(:,3),clusterX,chanEvent_list);
legend('location','northeast');
title('HFO Feature Space');
xlabel('Feature 1');
ylabel('Feature 2');
zlabel('Feature 3');
set(gca,'fontsize',12);
axis square
set(gcf,'position',[352 174 1205 298]);

% select optimal cluster from minimum SPL as reference
[~,clsSig_idx] = sort((ClusterSPLAvg./max(ClusterSPLAvg)));
clsSig_idx = clsSig_idx(1);
fprintf('Suggested optimal cluster: %d, please confirm.\n',clsSig_idx);
% ask for confirmation
prompt = {'Please confirm or enter pHFOs cluster #'};
ttl = 'User input';
dims = [1 35];
definput = {num2str(clsSig_idx)};
answer = inputdlg(prompt,ttl,dims,definput);
clsSig_idx = str2num(answer{1});

% compose events from selected cluster
iCls_set = unique(clusterX);
iEvent = [];
for iCls = 1:length(clsSig_idx)
    iClsIdx = iCls_set(clsSig_idx(iCls));
    iEvent = [iEvent; find(clusterX==iClsIdx)];
end
sigEventLabel_auto = zeros(length(clusterX),1);
sigEventLabel_auto(iEvent) = 1;
sigEventPoolRaw_new = sigEventPoolRaw(iEvent,:);
sigEventPoolFilt_new = sigEventPoolFilt(iEvent,:);
sigEventPoolRawWhole_new = sigEventPoolRawWhole(:,:,iEvent);
sigEventPoolFiltWhole_new = sigEventPoolFiltWhole(:,:,iEvent);
chanEvent_list_new = chanEvent_list(iEvent);
jc_print_single(sprintf('Identified Cluster#: %d.',clsSig_idx));
jc_print_single(sprintf('Total event#: %d.',length(chanEvent_list_new)));
jc_print_single(sprintf('Channels included: '));
display(chanEvent_list_new);

% store cluster results
jc_io_save([pathOutput,'\3_HFOEvents_compact'],...
    'sigEventLabel_auto','-append');
close(fig);

%% loop over all events
% determine the reference label for all identified events
% (1) thfo on spike  - 1
% (2) fhfo on spike  - 2
% (3) spike or noise - 0

jc_print_single('Evaluating HFOs identification results...');
% determine detection results one by one
nEvent = length(sigEventLabel_auto);
sigEventLabel_ref = zeros(nEvent,1);
for iEvent = 1:nEvent
    % get epoch index of the current event
    iEvent_idx = sigEventIndex(iEvent);
    % get channel index and temporal location of the current event
    iEvent_ichan = find(strcmpi(chan_list,chanEvent_list{iEvent}));
    iEvent_itpos = sigEventTime(iEvent)/1000*fs;
    % find matched simulated epoch
    strAct_tmp = strAct{iEvent_idx};
    iEpoch_thfo_ichan = strAct_tmp.iChanTHFO;
    iEpoch_fhfo_ichan = strAct_tmp.iChanFHFO;
    % determine if the event matches THFO
    if sum(iEpoch_thfo_ichan == iEvent_ichan) > 0
        % determine if temporal location matches
        nTHFO = length(iEpoch_thfo_ichan);
        for iTHFO = 1:nTHFO
            iEpoch_thfo_itpos = strAct_tmp.iPosTHFO(iTHFO,:);
            if iEvent_itpos >= min(iEpoch_thfo_itpos) &&...
                    iEvent_itpos <= max(iEpoch_thfo_itpos)
                sigEventLabel_ref(iEvent) = 1;
                break
            end
        end
    end
    % determine if the event matches FHFO
    if sum(iEpoch_fhfo_ichan == iEvent_ichan) > 0
        % determine if temporal location matches
        nFHFO = length(iEpoch_fhfo_ichan);
        for iFHFO = 1:nFHFO
            iEpoch_fhfo_itpos = strAct_tmp.iPosFHFO{iFHFO};
            if iEvent_itpos >= min(iEpoch_fhfo_itpos) &&...
                    iEvent_itpos <= max(iEpoch_fhfo_itpos)
                if sigEventLabel_ref(iEvent) == 0
                    sigEventLabel_ref(iEvent) = 2;
                    break
                else
                    warning('Already assigned label. Please check.');
                    keyboard
                end
            end
        end
    end
end

%% loop over all epochs
% determine the simulated events
% nSimTHFO, nSimFHFO

nEpoch = length(strAct);
nSimTHFO = 0;
nSimFHFO = 0;
simTHFO_idx = [];
simFHFO_idx = [];
for iEpoch = 1:nEpoch
    strAct_tmp = strAct{iEpoch};
    % get epoch index containing thfo
    if ~isempty(strAct_tmp.iChanTHFO)
        nTHFO = length(strAct_tmp.iChanTHFO);
        nSimTHFO = nSimTHFO + nTHFO;
        simTHFO_idx = [simTHFO_idx; repmat(iEpoch,nTHFO,1)];
    end
    % get epoch index containing fhfo
    if ~isempty(strAct_tmp.iChanFHFO)
        nFHFO = length(strAct_tmp.iChanFHFO);
        nSimFHFO = nSimFHFO + nFHFO;
        simFHFO_idx = [simFHFO_idx; repmat(iEpoch,nFHFO,1)];
    end
end

%% evaluate identification results
% definition:
%%% TP, true positive, identified and true HFO
%%% FP, false positive, identified but false HFO
%%% FN, false negative, missed (detected or not) but true HFO
%%% TN, true negative, negative samples would be difficult to define
% implement sensitivity, precision (PPV), F1 score

% setup event labels
% simulated events
%%% nSimTHFO and nSimFHFO
%%% determined labels
tagTHFO = find(sigEventLabel_ref==1);
tagFHFO = find(sigEventLabel_ref~=1);
% identified events
tagPHFO = find(sigEventLabel_auto == 1);
tagFHFO = find(sigEventLabel_auto ~= 1);
% calculate condition numbers
if length(tagTHFO) > 1
    nTP = length(intersect(tagTHFO,tagPHFO));
    nFP = length(setdiff(tagPHFO,tagTHFO));
    nFN = nSimTHFO - nTP;
    % calculate evaluation metrics
    rTPR = nTP/(nTP+nFN);
    rPPV = nTP/(nTP+nFP);
    rF1 = 2*(rPPV*rTPR)/(rPPV+rTPR);
else
    nTP = 0;
    nFP = 0;
    nFN = 0;
    rTPR = nan;
    rPPV = nan;
    rF1 = nan;
end
% store the results
evalEventResults = [rTPR rPPV rF1 nTP nFP nFN];

%% save results
% store labeling results
evalEventResults_ref = evalEventResults;
jc_io_save([pathOutput,'\3_HFOEvents_compact'],...
    'sigEventLabel_ref','evalEventResults_ref','-append');
fprintf('Evaluation recall: %.f%%, precision: %.f%%, F1: %.f%%.\n',...
    rTPR*100,rPPV*100,rF1*100);
jc_print_block('Saved HFO event results.');

