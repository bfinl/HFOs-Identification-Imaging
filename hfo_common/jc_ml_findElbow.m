%% Optimization
% This function finds the elbow point of the input L-shaped curve.
% The basic idea is to use the identified elbow point as an optimal
% parameter, such as number of components, in various of apoplications.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function K = jc_ml_findElbow(param,thrd)

if nargin < 2
    thrd = 0.95;
end

% define change gradient
paramGrd = abs(gradient(param));
% define the gradient thrd - significant descend
paramGrdThrd = quantile(paramGrd,thrd);
% find the best index
K = find(paramGrd>=paramGrdThrd,1)-1;
if K<1
    K = [];
    warning('Elbow point NOT found.')
end

