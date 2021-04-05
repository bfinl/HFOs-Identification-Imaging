%% Noise estimation
% This function estimates the noise covariance of the input signal using
% PCA, principal component analysis. The basic idea is to first extract the
% dominant components as the target signal, and treat the remaining
% components as the noisy signals. The level of dominant components is
% measured by variance explained by the components of the input data. The
% default level is set to 95% of variance.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [C_noise,SNR] = noiseEstimation(phi_noise,nPCADiscard)

if nargin < 2
    nPCADiscard = 0.95;
end

% PCA to remove main activity
[~,~,~,~,S_all] = reduceDimension(phi_noise,1,false);
[~,~,~,~,S_act] = reduceDimension(phi_noise,nPCADiscard,false);

% use data variance for estimation
% empirical noise estimation
C_noise = sum(S_all) - sum(S_act);
if C_noise > 0
    SNR = sum(S_act)/C_noise;
    C_noise = ones(size(phi_noise,1),1)*C_noise;
else
    SNR = [];
    C_noise = ones(size(phi_noise,1),1);
end

end




