%% HFO detection
% This function clusters detected hfo events.
%
% Overall steps of this algorithm are as follows.
% 1) find out optimal cluster number based on elbow curve
% 2) run k-means for initialization
% 3) run gmm based on the initialization
%       - appropriate gmm configuration, optimized on feature set
% 4) evaluate the results by various means
%       - raw/filtered piled plots
%       - spatial activation maps
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


% for reproducibility
rng('default');
rng(3);

% compose feature set
feature = zscore([feature_1 feature_2]);

% find optimal cluster number
nCls = findClusterNumber(feature);

%% clustering

%%% k-means clustering
clusterX = kmeans(feature,nCls,'replicate',10,...
    'MaxIter',1000,'OnlinePhase','on');

%%% GMM clustering
% dim. reduction on features
[feature_pc,mu,V] = reduceDimension(feature,[],true);
if size(feature_pc,2) < 4
    feature_pc = feature;
end
% model configuration
feature_gmm = feature_pc; clusterX_gmm = clusterX; nCls_gmm = nCls;
[nK_opt,Sigma_opt,ShrdCov_opt,gm_opt,allConverge] = jc_ml_optimizeGMM(...
    feature_gmm,nCls_gmm,clusterX_gmm);
% gmm clustering
[clusterX,nlogL,P,logpdf,d2] = cluster(gm_opt,feature_gmm);

%%% plot clusters with channel label
% feature space in 3D/2D space
if size(feature_pc,2) > 2
    figure;
    gscatter3(feature_pc(:,1),feature_pc(:,2),feature_pc(:,3),...
        clusterX,chanEvent_list);
    legend('location','northeast');
    title('HFO Feature Space');
    xlabel('Principal Feature 1');
    ylabel('Principal Feature 2');
    zlabel('Principal Feature 3');
else
    figure;
    gscatter2(feature_pc(:,1),feature_pc(:,2),feature_pc(:,3),...
        clusterX,chanEvent_list);
    legend('location','northeast');
    title('HFO Feature Space');
    xlabel('Principal Feature 1');
    ylabel('Principal Feature 2');
    zlabel('Principal Feature 3');
end

%% plot clusters with piled events
[clusterX,eventRate_set] = plotClusterPile(sigEventPoolRaw,sigEventPoolFilt,...
    clusterX,chanEvent_list,chan_list,fs,pathChannel);

%% save all results
% save(fullfile(pathOutput,'5_HFOClusters.mat'));
jc_print_block('HFO clustering done.');

