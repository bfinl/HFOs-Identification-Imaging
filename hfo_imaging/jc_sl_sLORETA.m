%% Source imaging
% This function performs sLORETA, Standardized low resolution brain
% electromagnetic tomography, to solve electrophysiological source imaging
% estimation. This process is to estimate the cortical electrical source
% distribution from the scalp M/EEG measurements through leadfield matrix
% modeled from geometric brain model, such as BEM model from MRI images.
%
% References are listed as follows.
% Pascual-Marqui, R. D. (2002). "Standardized low-resolution brain
% electromagnetic tomography (sLORETA): technical details." Methods Find
% Exp Clin Pharmacol 24 Suppl D: 5-12.
%
% P. C. Hansen, Regularization Tools: A Matlab package for analysis and
% solution of discrete ill-posed problems, Numerical Algorithms, 6 (1994),
% pp. 1-35.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2017.01.16
% Initial implementation
%
% Rui Sun, Zhengxiang Cai
% 2018.03.02
% Added orientational estimation.
%
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [j_opt,cost_opt,lambda,S_j] = jc_sl_sLORETA(s,lambda,phi,S_j)

% Primary inputs ----------------------------------------------------------
% phi,leadfield,lambda
% leadfield
K = s.Lfd;
% number of sources
% not necessarily equals to size of leadfield matrix
nSource = s.nSource;
% mainly for paralle
% variables defined in paralle loop
% if not, defined in input structure
if nargin == 1
    lambda  = s.lambda;
    phi     = s.phi;
    S_j     = [];
elseif nargin == 2
    phi     = s.phi;
    S_j     = [];
elseif nargin == 3
    S_j     = [];
end
% -------------------------------------------------------------------------
% Additional inputs -------------------------------------------------------
% sensor space covariance, for normalization
% not used currently, needs more test
if isfield(s,'C_phi')
    C_phi = s.C_phi;
else
    C_phi = [];
end
% sensor space noise, for whitening
if isfield(s,'C_noise')
    C_noise = s.C_noise;
else
    C_noise = 1;
end
% sensor space SNR, alternative for noise correction
% usually not needed, set to 1
if isfield(s,'SNR')
    SNR = s.SNR;
else
    SNR = 1;
end
% method for hyperparameter (lambda) estimation
% usually use L-curve
if isfield(s,'paramMethod')
    strMethod = s.paramMethod;
else
    strMethod = 'lcurve';
end
% -------------------------------------------------------------------------
% Inverse method setup
% sensor and source size
nElectrode = size(phi,1);
nLfd = size(K,2);
% whiten sensor data
[phi_w,K_w] = jc_sl_whtnSensor(phi,K,C_noise);

% Main loop for solution ==================================================

% solve for lambda
% if solved, return lambda for further usage
if isempty(lambda)
    lambda = jc_sl_solveLambda(phi_w,K_w,ones(nLfd,1),strMethod);
    lambda = lambda/SNR;
    % normalization, for sLORETA
    if isempty(S_j)
        if ~isempty(C_phi)
            S_j_tmp = K_w'*((C_phi)\K_w);
        else
            S_j_tmp = K_w'/(K_w*K_w' + lambda*eye(nElectrode))*K_w;
        end
        % efficient code of S_j
        % diagonal if fixed orientation
        % block wise diagonal is rotational
        if nSource==nLfd
            S_j = diag(S_j_tmp).^(-0.5);
        else
            S_j = zeros(3,nSource);
            for i = 1:nSource
                j_idx = (i-1)*3+(1:3);
                S_j(:,j_idx) = S_j_tmp(j_idx,j_idx)^(-0.5);
            end
        end
    end
    % set solution to empty
    j_opt = []; cost_opt = [];
    return
end

% solve pure minimum norm problem
% whitened, lambda solved
[j_opt,cost_opt] = solver(phi_w,K_w,lambda);

% normalization, for sLORETA
if isempty(S_j)
    if ~isempty(C_phi)
        S_j_tmp = K_w'*((C_phi)\K_w);
    else
        S_j_tmp = K_w'/(K_w*K_w' + lambda*eye(nElectrode))*K_w;
    end
    % efficient code of S_j
    % diagonal if fixed orientation
    % block wise diagonal is rotational
    if nSource==nLfd
        S_j = diag(S_j_tmp).^(-0.5);
    else
        S_j = zeros(3,nSource);
        for i = 1:nSource
            j_idx = (i-1)*3+(1:3);
            S_j(:,j_idx) = S_j_tmp(j_idx,j_idx)^(-0.5);
        end
    end
end

if nSource==nLfd
    j_opt = j_opt.*S_j;
else
    j_amp = zeros(size(j_opt));
    for i = 1:nSource
        j_amp_idx = (i-1)*3+(1:3);
        j_amp(j_amp_idx) = S_j(:,j_amp_idx)*j_opt(j_amp_idx);
    end
    j_opt = j_amp;
end

% =========================================================================
end

%% Sub-functions
% solver
function [J,cost] = solver(phi,K,lambda)

    % this solver uses matrix formulation of sLORETA
    M = size(K,1);

    J = K'/(K*K'+lambda*eye(M))*phi;
    cost = [norm( K * J - phi, 2 ) norm( J, 2 )];

end

% whiten sensor space data
function [phi_w,K_w] = jc_sl_whtnSensor(phi,K,C_noise)

    % C_noise is the covariance matrix for sensor noise space
    if isvector(C_noise)
        if isscalar(C_noise)
            C_noise = repmat(C_noise,size(phi,1),1);
        end
        if ~iscolumn(C_noise)
            C_noise = C_noise(:);
        end
        phi_w = phi./repmat(sqrt(C_noise),1,size(phi,2));
        K_w   = K./repmat(sqrt(C_noise),1,size(K,2));
    else
        [U,S,~] = svd(C_noise,0);
        C_noise_2 = ((S + 0.001*S(1,1)*eye(size(S)))^(-0.5))*(U.');
        phi_w = C_noise_2*phi;
        K_w   = C_noise_2*K;
    end
    
end

% solve for lambda, imaging hyperparameter
function [lambda,lambda_norm] = jc_sl_solveLambda(phi_w,K_w,W,strMethod)

    % de-weight
    if isvector(W)
        K_w_dw = K_w.*repmat(1./W',size(K_w,1),1);
    else
        K_w_dw = K_w*W;
    end
    % normalize Lfd
    [K_w_dw,K_scale] = sl_normLfd(K_w_dw);

    % svd
    [U,s,V] = csvd(K_w_dw);

    % solve for lambda
    switch strMethod
        case 'lcurve'
            [lambda_norm_1,~,~,~] = l_curve(U,s,phi_w);
            [~,lambda_norm_2] = discrep(U,s,V,phi_w,sqrt(chi2inv(0.99,length(phi_w))));
            lambda_norm = max(lambda_norm_1,lambda_norm_2);
        case 'dcp'
            [~,lambda_norm] = discrep(U,s,V,phi_w,sqrt(chi2inv(0.99,length(phi_w))));
    end

    % lambda calculated from l-curve is for rho and eta, not rho^2 and eta^2
    lambda = (lambda_norm*K_scale)^2;

end

% scale leadfield
function [Lfd,scale] = sl_normLfd(Lfd0)

    % scale lead-field matrix with its maximum
    scale = max(max(Lfd0));
    Lfd = Lfd0/scale;

end


