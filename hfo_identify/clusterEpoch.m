%% Clustering
% This function clusters epochs according to input maps using k-means method.
%
% Input:
%   (1) sig_map:
%       - sig(nSegemtns, nSamples)
%   (2) nClsMax:
%       - number of clusters
%   (3) pathChan:
%       - path of electrode digitizer
%   (4) clsEventMap:
%      - reference map for calculating relative relavance of clustered maps
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [clsSigSort,nCls,optCls,fig] = clusterEpoch(sig_map,nClsMax,pathChan,clsEventMap)

% signal dimension
[nChan,nComp] = size(sig_map);

if nargin < 2 || isempty(nClsMax)
    nClsMax = round(sqrt(nComp));
end

% prepare feature
sig_map_vect = normSig(sig_map);
sig_map_feat = sig_map_vect';

% For reproducibility
rng('default');
rng(3);
% k-means clustering
[clsSig,~,~,nCls] = kmeansElbow(sig_map_feat,nClsMax,0.95,20,'cosine');
% note: use 'cosine', since topo vector can be non-zero sum, we would like
% to have similarity between vectors without extracting the mean

% sort optimal cluster according to the reference map
% compare the cluster map to the scalp distribution of the cluster events
clsEventMap_norm = normSig(clsEventMap);
clsEventMap_feat = clsEventMap_norm';
clsLbd = zeros(nCls,1);
for iCls = 1:nCls
    clsCompMap = sig_map_feat(clsSig == iCls,:);
    clsLbd_tmp = simCosine(clsEventMap_feat,clsCompMap);
    clsLbd(iCls) = mean(clsLbd_tmp);
end
[~,clsIdx] = sort(clsLbd,'descend');

% plot maps to visualize clusters
% define class labels
clsSig_set = cell(nCls,1);
nClsSig_set = zeros(nCls,1);
clsSigSort = zeros(size(clsSig));
for iCls = 1:nCls
    clsSig_set{iCls} = find(clsSig==clsIdx(iCls));
    nClsSig_set(iCls) = length(clsSig_set{iCls});
    clsSigSort(clsSig_set{iCls}) = iCls;
end
optCls = clsSig_set{1}; % optimal cluster

nClsSig_max = max(nClsSig_set)+1;
% plot maps
fig = figure;
ha = tightSubplot(nCls,nClsSig_max,[.01 .01],[.01 .01],[.01 .01]);
axis(ha,'off');
for iCls = 1:nCls
    clsSub_idx = clsSig_set{iCls};
    nClsSub = length(clsSub_idx);
    cls_sub_idx_0 = (iCls-1)*nClsSig_max + 1;
    % plot class label
    axes(ha(cls_sub_idx_0));
    text(0,0,['Class',num2str(iCls)]);
    set(gca,'xlim',[0 1],'ylim',[-1 1],'fontweight','bold');
    axis off;
    for iFig = 1:nClsSub
        % plot cluster map
        axes(ha(iFig+cls_sub_idx_0));
        plotTopoMap(sig_map_vect(:,clsSub_idx(iFig)),fig,pathChan,[]);
        text(-0.5,-0.5,num2str(clsSub_idx(iFig)));
    end
end

end


%% sub-function
function sig_new = normSig(sig)
    sig_norm = sqrt(sum(sig.^2,1));
    sig_new  = abs(sig)./repmat(sig_norm,size(sig,1),1);
end

function sim = simCosine(a,B)
    nCol = size(B,1);
    sim = zeros(nCol,1);
    for iCol = 1:nCol
        b = B(iCol,:);
        sim(iCol,1) = (a(:)'*b(:)) / (sqrt(sum(a.^2))*sqrt(sum(b.^2)));
    end
end

    
    
    