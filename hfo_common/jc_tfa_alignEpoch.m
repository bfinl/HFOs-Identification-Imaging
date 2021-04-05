%% Time-frequency analysis
% Auto-Align Multi-Channel Segments
%
% This function aligns segments according to the cross-correlation.
% First, define correlation matrix across all signal samples, then find
% the maximum peaks in correlation for all signal samples, and align
% them regarding to their own alignment point.
% Finally, the signal segments are pruned to have the same length
%
% Input:
%   (1) sigOrigSeg:
%       - data segments in 3d dimension (nChan x nSample x nSegments)
%   (2) nChan:
%       - number of channels in sigOrigSeg
%   (3) nWin:
%       - sample window of interest, in sample points
% Output:
%   (1) sigAlignedSeg:
%       - data segments aligned in the same 3d dimension
%   (2) sig_peak_idx:
%       - relative sample points of alignment
%
% In general, such alignment could facilitate averaging across several event
% segements, which might potentially improve the SNR of signals of interest.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [sigAlignedSeg,sig_peak_idx] = jc_tfa_alignEpoch(sigOrigSeg,nChan,nWin)

if nargin < 3
    nWin = inf;
end

% check data size
[mSig,nSig] = size(sigOrigSeg(:,:,1));
if nChan ~= mSig
    if nSig == nChan
        sigOrigSeg = permute(sigOrigSeg,[2 1 3]);
    else
        error('Signal size is not correct.');
    end
end
% block samples outside window of interest
[~,nSamp,nSeg] = size(sigOrigSeg);
sigOrigSeg_tmp = sigOrigSeg;
if nWin < nSamp
    for iSeg = 1:nSeg
        nWinIdx = floor(nSamp/2 + (-nWin:nWin)/2);
        replaceSamp1 = sigOrigSeg_tmp(:,nWinIdx(1),iSeg);
        replaceSamp2 = sigOrigSeg_tmp(:,nWinIdx(end),iSeg);
        sigOrigSeg_tmp(:,1:nWinIdx(1),iSeg) = repmat(replaceSamp1,1,nWinIdx(1));
        sigOrigSeg_tmp(:,nWinIdx(end):end,iSeg) = repmat(replaceSamp2,1,nSamp-nWinIdx(end)+1);
    end
end

% compute cross-correlation
[nChan,nSamp,nSeg] = size(sigOrigSeg_tmp);
sigXCORR_peak = zeros(nSeg);
sigXCORR_peak_idx = zeros(nSeg);
maxLag = min([floor(nSamp/10),30,nWin]);
sigXCORR_lag = -(maxLag):(maxLag);

% loop through all segments
for iSeg = 1:nSeg
    jSeg_set = 1:nSeg; jSeg_set(jSeg_set==iSeg) = [];
    % loop through all other segemnts to define best aligned
    % signals using XCORR
    for jSeg = jSeg_set
        % calculate XCORR between the current two segements
        sigXCORR = zeros(1,2*maxLag+1);
        for iChan = 1:nChan
            % it is possible that the epochs would align at large
            % lags, which brings the aligned segmnent very short (truncated)
            % as we assume that the target signal is relatively around the
            % center, we limit the maximum lag to be 1/2 of the length
            % from -maxLag to +maxLag
            sigXCORR = sigXCORR +...
                xcorr(sigOrigSeg_tmp(iChan,:,iSeg),sigOrigSeg_tmp(iChan,:,jSeg),...
                floor(maxLag),'coef');
        end
        [sigXCORR_peak(jSeg,iSeg), sigXCORR_peak_idx(jSeg,iSeg)] = max(sigXCORR);
    end
end
% find best set of alignment using sum of XCORR
[~,sigXCORR_peak_max_idx] = max(sum(sigXCORR_peak,1));
sigXCORR_peak_idx_max = sigXCORR_peak_idx(:,sigXCORR_peak_max_idx);
sig_peak_idx = zeros(size(sigXCORR_peak_idx_max));
sig_peak_idx(sigXCORR_peak_idx_max~=0)...
    = -sigXCORR_lag(sigXCORR_peak_idx_max(sigXCORR_peak_idx_max~=0));

% convert to index (sample point)
sig_peak_idx = sig_peak_idx + floor(nSamp/2);

% align signals
sigAlignedSeg = jc_tfa_alignEpochIdx(sigOrigSeg,sig_peak_idx,nChan,nSamp,nSeg,'zp');

end






