%% Dimensionality reduction
% This function reduces dimension of input data using PCA, principal
% component analysis.
%
% The process of dimensionality reduction takes input data and trandforms
% it into a reduced dimensional space. The main purpose of this function is
% to perform dimensionality reduction via various methods, such as PCA.
% This process could also be integrated with scaling factors to standardize
% the data components to variance 1 after reducing the dimensionality.
% Here, we implemented a commonly used whitening technique, PCA, principal
% component analysis.
% 
% Inputs:
%   (1) data_raw: matrix, samples (observations) X features (predictors)
%           For instance, each row is one example or observation, and each
%           column is one feature or measure kind.
%   (2) n_comp (two options):
%           (1) int, number of components to keep
%           (2) real, [0,1], percentage of variance to keep, DEFAULT: 1
%   (3) f_whtn: boolean, Flag for whitening the data after dimensionality
%   reduction, DEFAULT: false, not whitening
%
% Outputs:
%   (1) data_rd: PCA transformed data
%   (2) mu, mean of original data
%   (3) V, PCA transformation matrix
%   (4) n_comp, number of components kept in output data
%   (5) S: singular values of the components in the output data
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [data_rd,mu,V,n_comp,S] = reduceDimension(data_raw,n_comp,f_whtn)

global EPSILON
EPSILON = 10e-8;

% default dimensions to keep
% two options for input
%   (1) 0 <= n_comp <= 1, variance to keep
%   (2) 1 < n_comp <= full rank, components to keep
if ~exist('n_comp','var') || isempty(n_comp)
    n_comp = 1; % keep all variance
end
% default whitening option
if ~exist('f_whtn','var') || isempty(f_whtn)
    f_whtn = false;
end

% standardization
% principal component analysis
%
% compute eigen basis for the feature space
[U,S,data_new,mu] = computeEigenBasis(data_raw);
% compute dimensions to keep, if not given
if (n_comp >= 0) && (n_comp <= 1)
    n_comp = findEigenBasisByVariance(S,n_comp);
end
% keep given components
U = U(:,1:n_comp); S = S(1:n_comp);
% compute the PCA whitened data
% epsilon = 10e-8;
% whitening option
if f_whtn
    V = bsxfun(@times,U,1./sqrt(S + EPSILON));
else
    V = U;
end
% compute the PCA whitened data
data_rd = data_new * V;

end


% sun-function
% compute eigen decomposition using SVD
function [U,S,dataZeroMean,mu] = computeEigenBasis(dataRaw)

	% get size of data
    % [N,P] = size(dataRaw);
    % assert(N >= P); % samples should exceed predictors
    
    % compute mean data value separately for each column or feature
    mu = mean(dataRaw,1);
    % subtract the mean from the data in each feature column
    % dataZeroMean = dataRaw - repmat(mu,size(dataRaw,1),1);
    dataZeroMean = bsxfun(@minus,dataRaw,mu);
    % compute the covariance matrix
    dataSigma = cov(dataZeroMean);
    
    % compute eigenvectors and eigenvalues using SVD
    %
    % Matrix U contains the eigenvectors of ? (one eigenvector per
    % column, sorted in order from top to bottom eigenvector), and the
    % diagonal entries of the matrix S contains the corresponding
    % eigenvalues (also sorted in decreasing order). The matrix V is
    % equal to U, and can be safely ignored.
    [U,S,~] = svd(dataSigma);
    S       = diag(S);
    S       = S(:)';
    
    % figure out which values of S are non-zero.
    tol = eps(class(dataRaw));
    idx = (S > max(S)*tol);
    assert(~all(idx == 0));
    % get the non-zero elements of Sig and corresponding columns of U.
    S = S(idx);
    U = U(:,idx);
    
    % compute prewhitened data.
    % Z = bsxfun(@times,Z*U,1./sqrt(S));
    
end

% find index for components that count for some threshold of variance
function nComp = findEigenBasisByVariance(lambda,thrd)

    global EPSILON

    if ~isvector(lambda)
        lambda = diag(lambda);
    else
        lambda = lambda(:);
    end
    [lambda_sort,lambda_idx] = sort(lambda,'descend');

    lambda_ratio = lambda_sort/sum(lambda_sort);
    lambda_ratio_cum = lambda_ratio'*triu(ones(length(lambda)));

    % find 90% and 99% points or defined by thrd
    lambda_ratio_cum_idx = find(lambda_ratio_cum > thrd-EPSILON,1);
    nComp = lambda_idx(1:lambda_ratio_cum_idx);
    nComp = max(nComp);

    % plot component selection result, if needed
    if false
        figure, hold all; %#ok<UNRCH>
        plot(lambda,'*-','linewidth',1);
        plot(1:nComp,lambda(1:nComp),'o','linewidth',1.5);
        xlabel('Component Index');
        ylabel('\lambda');
        box on;
        set(gca,'fontsize',12);
    end

end

