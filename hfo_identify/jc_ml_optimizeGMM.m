%% Clustering
% Perform GMM clustering with optimization.
% This function evaluates GMM clustering using AIC/BIC criteria to
% determine the optimal clustering model.
%
% Reference from Matlab instruction for Tuning GMM.
% Suppose k is the number of desired components or clusters, and  is the
% covariance structure for all components. Follow these steps to tune a
% GMM.
%   1. Choose a (k,Sigma) pair, and then fit a GMM using the chosen
%   parameter specification and the entire data set.
%   2. Estimate the AIC and BIC.
%   3. Repeat steps 1 and 2 until you exhaust all (k,Sigma) pairs of
%   interest.
%   4. Choose the fitted GMM that balances low AIC with simplicity.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [nK_opt,Sigma_opt,ShrdCov_opt,gm_opt,allConverge] = jc_ml_optimizeGMM(...
    X,nCls,initClsIdx,nReplicate)

if nargin < 2 || isempty(nCls)
    nCls = 1:size(X,1)-1;
end
if nargin < 3 || isempty(initClsIdx)
    initClsIdx = 'plus';
end
if nargin < 4 || isempty(nReplicate)
    if isvector(initClsIdx)
        nReplicate = 1;
    else
        nReplicate = 10;
    end
end

% setup models
k = nCls;
nK = numel(k);
Sigma = {'diagonal','full'};
nSigma = numel(Sigma);
SharedCovariance = {true,false};
SCtext = {'true','false'};
nSC = numel(SharedCovariance);
RegularizationValue = 1e-8; % ensure positive-definite covariance
options = statset('MaxIter',10000);

% Preallocation
gm = cell(nK,nSigma,nSC);
aic = zeros(nK,nSigma,nSC);
bic = zeros(nK,nSigma,nSC);
converged = false(nK,nSigma,nSC);

% For reproducibility
rng('default');
rng(3);

% Fit all models
for m = 1:nSC
    for j = 1:nSigma
        for i = 1:nK
            gm{i,j,m} = fitgmdist(X,k(i),...
                'CovarianceType',Sigma{j},...
                'SharedCovariance',SharedCovariance{m},...
                'RegularizationValue',RegularizationValue,...
                'Start',initClsIdx,...
                'Replicates',nReplicate,...
                'Options',options);
            aic(i,j,m) = gm{i,j,m}.AIC;
            bic(i,j,m) = gm{i,j,m}.BIC;
            converged(i,j,m) = gm{i,j,m}.Converged;
        end
    end
end

% plot, if needed
% figure, subplot(121);
% bar(reshape(aic,nK,nSigma*nSC));
% title('AIC For Various $k$ and $\Sigma$ Choices','Interpreter','latex');
% xlabel('$k$','Interpreter','Latex');
% ylabel('AIC');
% legend({'Diagonal-shared','Full-shared','Diagonal-unshared',...
%     'Full-unshared'});
% 
% subplot(122);
% bar(reshape(bic,nK,nSigma*nSC));
% title('BIC For Various $k$ and $\Sigma$ Choices','Interpreter','latex');
% xlabel('$c$','Interpreter','Latex');
% ylabel('BIC');
% legend({'Diagonal-shared','Full-shared','Diagonal-unshared',...
%     'Full-unshared'});

% evaluate criteria
allConverge = (sum(converged(:)) == nK*nSigma*nSC);
abic = zscore(aic(:)) + zscore(bic(:));
[~,abic_loc] = min(abic(:));
[abic_idx(1),abic_idx(2),abic_idx(3)] = ind2sub(size(aic),abic_loc);
nK_opt = k(abic_idx(1));
Sigma_opt = Sigma{abic_idx(2)};
ShrdCov_opt = SharedCovariance{abic_idx(3)};
gm_opt = gm{abic_idx(1),abic_idx(2),abic_idx(3)};

% display results
fprintf('Optimal GMM configuration: %d Cls with %s & %s shared COV.\n',nK_opt,Sigma_opt,SCtext{abic_idx(3)});

end



