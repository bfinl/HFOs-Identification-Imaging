%% HFO detection
% This function labels clustered HFO events and extracts HFO activities.
%
% Overall steps of this algorithm are as follows.
% 1) construct hfo clusters
% 2) extract hfo events
% 3) save events and maps
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
% the robustness of such measurement and implementations.
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

%% estimate temporal-spectral range
event_idx = 1:length(chanEvent_list_new);
[nEvent,nEventLen] = size(sigEventPoolRaw_new(event_idx,:));
% determine hfo event spectra range
[sigEventPoolSampLoc,sigEventPoolSpecRange] = findSpectraMax(sigEventPoolFilt_new,fs);
fprintf('Event frequency at [%.2f %.2f] Hz.\n',sigEventPoolSpecRange);
fprintf('Event location at [%.0f %.0f] Sample.\n',sigEventPoolSampLoc);

% raw events - unfiltered
sigEvent = reshape(sigEventPoolRawWhole_new(:,:,event_idx),nChan,[])';
% filtered events
% sigEvent = reshape(sigEventPoolFiltWhole_new(:,:,event_idx),nChan,[])';
clsEventMap = sum(cat(2,eventRate_set{clsSig_idx}),2);

%% spatio-spectral decomposition, ssd

% config parameters
fCut = [sigEventPoolSpecRange;
        sigEventPoolSpecRange+[-3 3];
        sigEventPoolSpecRange+[-1 1]];
eventIdx = zeros(nEvent,2);
sigEventPoolSampLoc_set = repmat(sigEventPoolSampLoc,nEvent,1);
for iEvent = 1:nEvent
    eventIdx(iEvent,:) = sigEventPoolSampLoc_set(iEvent,:) + (iEvent-1)*nEventLen;
end
% perform ssd
[W, A, lambda, C_s, X_ssd] = ssd(sigEvent, fCut, fs, [], eventIdx);
% print results
jc_print_single(sprintf('Estimated total SSD comp.s#: %d.',size(A,2)));

%% find dominant components
%%% In this part, the SSD components are clustered to find the optimal ones
%%% that are concordant with the identified events. In general, most of the
%%% dominant components should be consistent with the identified events,
%%% since SSD is operated on the spatial and spectral information of the
%%% events. While some components could include various noise, the
%%% components could be selected by comparing the cosine distance of the
%%% components to the identified spatial distribution of the positive
%%% events.

jc_print_single(sprintf('Finding optimal components...'));
% elbow criterion to find dominant components
nComp = jc_ml_findElbow(lambda,0.95);
% % plot component selection result
% fig1 = figure; hold all;
% plot(lambda,'*-','linewidth',1);
% plot(1:nComp,lambda(1:nComp),'o','linewidth',1.5);
% xlabel('Component Index');
% ylabel('\lambda');
% title('Dominant SSD Components');
% box on;
% fprintf('Event Components to examine: %d.\n',nComp);

% cluster dominant ssd components
A_new = A(:,1:nComp);
[clsSig,nCls,clsIdxOpt,fig] = clusterEpoch(A_new,[],pathChannel,clsEventMap);
set(fig,'position',[680 50 560 420]);

% ask for confirmation
%%% Please confirm that the clustered components are concordant ones.
%%% In general, optimal cluster is class 1, which is recommended for
%%% subsequent analysis. The user could also manually identify components.
prompt = {'Please confirm or enter optimal components #'};
ttl = 'Input';
dims = [1 35];
definput = {num2str(clsIdxOpt')};
answer = inputdlg(prompt,ttl,dims,definput);
clsIdx = str2num(answer{1});
fprintf('Event Components to keep: %d.\n',length(clsIdx));
close all;

%% compose multichannel event data by optimal cluster

jc_print_single(sprintf('Extracting HFOs components...'));
% reconstruct dominent components
X_s_main = A(:,clsIdx)*X_ssd(:,clsIdx)';
nLenMax = max(sigEventPoolSampLoc(:,2)-sigEventPoolSampLoc(:,1)+1);
X_s_main_set = zeros(size(X_s_main,1),nLenMax,size(eventIdx,1));
for i = 1:size(eventIdx,1)
    X_s_main_set(:,1:length(eventIdx(i,1):eventIdx(i,2)),i) = X_s_main(:,eventIdx(i,1):eventIdx(i,2));
end
[~,alignIdx] = jc_tfa_alignEpoch(X_s_main_set(clsEventMap>0,:,:),sum(clsEventMap>0));
X_s_main_avg = mean(jc_tfa_alignEpochIdx(X_s_main_set,alignIdx,nChan,size(X_s_main_set,2),size(X_s_main_set,3)),3);

% plot topo-map
figure;
subplot(1,2,1);
plotTopoMap(var(X_s_main_avg,[],2),gcf,pathChannel,[]);
axis square;
title('Event topo-map');
subplot(1,2,2)
plot(X_s_main_avg')
axis square
title('Event activity');

%% save hfo event classes and maps
 
% compose results
hfoEvent = struct(...
    'hfoEventRaw',sigEventPoolRaw_new,...
    'hfoEventFilt',sigEventPoolFilt_new,...
    'hfoEventRawWhole',sigEventPoolRawWhole_new,...
    'hfoEventFiltWhole',sigEventPoolFiltWhole_new,...
    'hfoEventChan',chanEvent_list_new,...
    'nChan',nChan,...
    'nSamp',nEventLen,...
    'nEvent',nEvent);
hfoSig_set = X_s_main;
hfoSig_loc = eventIdx;
hfoSig = X_s_main_avg;
if size(hfoSig,2)>1
    hfoMap = var(hfoSig,[],2);
else
    hfoMap = hfoSig;
end

% save data file
fileNameEEG = fullfile(pathOutput,'\HFOResults');
jc_io_save(fileNameEEG,'hfoEvent','hfoSig','hfoMap','hfoSig_set','hfoSig_loc','chan_list');
% save averaged spike data
save(fullfile(pathOutput,'SpikeAvg.mat'),'sigEpochSpike','sigEpochSpikeAvg')

% save all results
% save(fullfile(pathOutput,'6_HFOResults.mat'));
jc_print_block('Saved HFO event results.');


