%% Visualization
% This function plots piled events for each cluster.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [clusterY,eventRate_set,meanEventFreq_set] = plotClusterPile(...
    sigEventRaw,sigEventFilt,clusterX,chanEvent_list,...
    chan_list,fs,pathChannel)

nChan = length(chan_list);
iCls_set = unique(clusterX);
nCls = length(iCls_set);

% preallocate for event cells
sigEventPoolRaw_class  = cell(nCls,1);
sigEventPoolFilt_class = cell(nCls,1);
chanEvent_list_class = cell(nCls,1);
clusterY = zeros(size(clusterX));
for iCls = 1:nCls
    clsIsx = iCls_set(iCls);
    iEvent = find(clusterX==clsIsx);
    
    sigEventPoolRaw_class{iCls}  = sigEventRaw(iEvent,:);
    sigEventPoolFilt_class{iCls} = sigEventFilt(iEvent,:);
        
    chanEvent_list_class{iCls} = chanEvent_list(iEvent);

    clusterY(iEvent) = iCls;
end

% evaluate on clusters
% regular metrics
% mean frequency and channel rates
clsIdx_set = num2cell(1:nCls);
meanEventFreq_set = cell(size(clsIdx_set));
eventRate_set = cell(size(clsIdx_set));
for iClsIdx = 1:length(clsIdx_set)
    
    clsIsx = clsIdx_set{iClsIdx};
    % construct class-set data
    sigEventPoolRaw_class_set  = cell2mat(sigEventPoolRaw_class(clsIsx));
    sigEventPoolFilt_class_set = cell2mat(sigEventPoolFilt_class(clsIsx));
    chanEvent_list_class_set = cat(1,chanEvent_list_class{clsIsx});

    [nEvent,~] = size(sigEventPoolRaw_class_set);
    % event mean frequency
    meanEventFreq = zeros(nEvent,2);
    for iEvent = 1:nEvent
        meanEventFreq(iEvent,1) = meanfreq(sigEventPoolFilt_class_set(iEvent,:),fs);
        meanEventFreq(iEvent,2) = meanfreq(sigEventPoolRaw_class_set(iEvent,:),fs);
    end
    
    % event rates per channels
    eventRate = zeros(nChan,1);
    for iEvent = 1:nEvent
        eventChan_idx = find(strcmp(chan_list,chanEvent_list_class_set{iEvent}));
        eventRate(eventChan_idx) = eventRate(eventChan_idx)+1;
    end
    
    % restore results
    meanEventFreq_set{iClsIdx} = meanEventFreq;
    eventRate_set{iClsIdx} = eventRate;
    
end

% plot clustered events profile
figure, hold all;
meanEventFreq_hf = cat(1,meanEventFreq_set{:});
meanEventFreq_lf = meanEventFreq_hf(:,2);
meanEventFreq_hf = meanEventFreq_hf(:,1);
eventLabel = zeros(length(meanEventFreq_lf),1);
for iCls = 1:length(meanEventFreq_set)
    nEvent = size(meanEventFreq_set{iCls},1);
    eventLabel(find(eventLabel==0,1)+(0:nEvent-1),1) = iCls;
end
boxplot(meanEventFreq_hf,eventLabel);
xlabel('Class Index');
ylabel('Event Mean Frequency (Hz)');

% plot pile morophology
% plot event rates
figure, hold all;
tSig = (0:1:length(sigEventPoolRaw_class{1}(1,:))-1).*1000/fs;
ySigRaw_max = max(max(sigEventRaw)); ySigRaw_max = ySigRaw_max*(1+sign(ySigRaw_max)*0.05);
ySigRaw_min = min(min(sigEventRaw)); ySigRaw_min = ySigRaw_min*(1-sign(ySigRaw_min)*0.05);
ySigFilt_max = max(max(sigEventFilt))/2; ySigFilt_max = ySigFilt_max*(1+sign(ySigFilt_max)*0.05);
ySigFilt_min = min(min(sigEventFilt))/2; ySigFilt_min = ySigFilt_min*(1-sign(ySigFilt_min)*0.05);
for iClsIdx = 1:length(clsIdx_set)
    subplot(3,length(clsIdx_set),iClsIdx)
    plot(tSig,sigEventPoolRaw_class{iClsIdx}');
    title(['Cluster ',num2str(iClsIdx)],'fontname','Times New Roman');
    set(gca,'xlim',[min(tSig),max(tSig)],'ylim',[ySigRaw_min ySigRaw_max],...
    'ytick',sort(floor([ySigRaw_min/2 0 ySigRaw_max/2])),...
    'box','on','linewidth',1,'fontname','Times New Roman','fontsize',12);
    if iClsIdx == 1
        ylabel('Amp. (\muV)');
    end
    
    subplot(3,length(clsIdx_set),iClsIdx+length(clsIdx_set))
    plot(tSig,sigEventPoolFilt_class{iClsIdx}');
    set(gca,'xlim',[min(tSig),max(tSig)],'ylim',...
        [ySigFilt_min ySigFilt_max],'ytick',sort(round([ySigFilt_min/2 0 ySigFilt_max/2])),...
        'box','on','linewidth',1,'fontname','Times New Roman','fontsize',12);
    xlabel('t (ms)');
    if iClsIdx == 1
        ylabel('Amp. (\muV)');
    end
    
    subplot(3,length(clsIdx_set),iClsIdx+2*length(clsIdx_set))
    plotTopoMap(eventRate_set{iClsIdx},gcf,pathChannel,[],[]);
    set(gca,'box','on','linewidth',1,'fontname','Times New Roman','fontsize',12);
end
savefig(gcf,'5_cluster_pile_topo.fig','compact');

end

