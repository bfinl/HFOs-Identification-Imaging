%% HFO detection
% Cluster uniformity measure
% This function calculates within cluster temporal homogeneity measure
% which describes the uniformity of HFOs with regard to spikes as reference
% for selecting optimal clusters.
%
% It is noteworthy that the cluster with HFOs riding on spikes is generally
% clear to identify by plotted spatiotemporal patterns. This measure is
% only a reference for selecting such clusters. This can also be integrated
% for automated selection as well, but more evaluation is recommended for
% the robustness of such measurement and implementations.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [ClusterSPLAvg,ClusterSPL] = computeClsUniMeasure(...
    sigEventTime,clusterX,sigEpochInfo)

iCls_set = unique(clusterX);
nCls = length(iCls_set);
% temporal homogeneity
%%% use spike duration center as reference
sigEpochCentr = sigEpochInfo.Durations/2;
ClusterSPL    = zeros(nCls,1);
ClusterSPLAvg = zeros(nCls,1);
for iCls = 1:nCls
    iClsIdx = iCls_set(iCls);
    clsEventIdx = find(clusterX == iClsIdx);
    clsEventTime_set  = sigEventTime(clsEventIdx,:);
    % calculate, temporal homogeniety between low/high signals
    ClusterSPLAvg(iCls) = std(abs(clsEventTime_set - sigEpochCentr));
    ClusterSPL(iCls) = ClusterSPLAvg(iCls)*(length(clsEventIdx)-1);
end

% discard one/two-element clusters
iCls_set = unique(clusterX);
nCls = length(iCls_set);
for iCls = 1:nCls
    iClsIdx = iCls_set(iCls);
    clsEventIdx = find(clusterX == iClsIdx);
    if length(clsEventIdx) < 3
       ClusterSPLAvg(iCls) = max(ClusterSPLAvg)+1;
    end
end

end
