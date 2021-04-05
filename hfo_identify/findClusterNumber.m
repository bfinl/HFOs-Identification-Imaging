%% Clustering
% Evaluate and optimize the number of clusters
% This function runs through several evaluation criteria to compare and
% determine the optimal cluster number.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function nCls = findClusterNumber(feature,varargin)

% for reproducibility
rng('default');
rng(3);

% setup evaluation criteria
eva_criterion = {'CalinskiHarabasz','DaviesBouldin','gap','silhouette'};
n_eva_crt = length(eva_criterion);

% clustering function
nReplicate = 10;
options.UseParallel = false;
myfunc = @(X,K)(kmeans(X, K,'replicate',nReplicate,...
    'MaxIter',1000,'OnlinePhase','on','options',options,varargin{:}));
% estimation of the upper bound for cluster number
nCls_max = round(sqrt(size(feature,1)));

% run through all evaluation criteria
nCls_set = zeros(n_eva_crt,1);
for i_eva_crt = 1:n_eva_crt
    fprintf('Calculating optimal cluster #...%d%%\n',100*i_eva_crt/(n_eva_crt+1));
    eva = evalclusters(feature,myfunc,eva_criterion{i_eva_crt},...
        'klist',(1:nCls_max));
    nCls_set(i_eva_crt) = eva.OptimalK;
end
% elbow method
fprintf('Calculating optimal cluster #...%d%%\n',100);
[~,~,~,nCls_add] = kmeansElbow(feature,nCls_max,0.95,nReplicate);

% find optimal cluster number by distribution
nCls_set = [nCls_set;nCls_add];
nCls = round(median(nCls_set));

% output result
% nCls_set_lbl = [eva_criterion,{'Elbow'}];
% figure, bar(nCls_set); set(gca,'xticklabel',nCls_set_lbl);
fprintf('Optimal cluster # found: %d.\n', nCls);

end



