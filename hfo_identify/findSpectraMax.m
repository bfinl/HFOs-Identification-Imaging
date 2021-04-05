%% Time-frequency analysis
% This function uses STFT to find maximum spectrum range and time location.
%
% This function allows signal to be multi-events
% in this case, the spectrum is estimated in an averaged sense
% temporal location will be find individually for each events
% 
% This function is also used to prepare parameters for SSD analysis.
% Some recomendations for optimizing SSD parameters, by Cohen, 2017
% Cohen, M. X. (2017). "Comparison of linear spatial filters for
% identifying oscillatory activity in multichannel data." J Neurosci
% Methods 278: 1-12.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [tempRange,specRange,tempPeak] = findSpectraMax(sigEvent,fs,tWin)

if nargin < 3
    % target signal of interest
    tWin = 100/1000; % ms
end

[nEvent,nSamp] = size(sigEvent);
% window size for default event location
% recomended 2Hz as flanking range, 6 Hz away from center frequency of interest
nWin = tWin*fs; % pts
fWin = 12; % Hz
nInitSpot = floor(nSamp/2);
% configurations for STFT
nWinSTFT = 32;
nShiftSTFT = 1;
nfftSTFT = 1024;
% perform STFT
specSum = 0;
for iEvent = 1:nEvent
    sig_tmp = sigEvent(iEvent,:);
    [tSpec,fSpec,spec] = stftSpectrum(sig_tmp,fs,'Ham',nWinSTFT,nShiftSTFT,nfftSTFT,false);
    specSum = specSum + spec;
end
specAvg = specSum/nEvent;
% zero out non-interested-regions
specAvg(:,[1:floor(nInitSpot-nWin/2),floor(nInitSpot+nWin/2):end]) = 0;

% temporal and spectral peak
[maxSpec,~] = max(specAvg);
[~,tMaxSpec] = max(maxSpec);
[maxSpec,~] = max(specAvg,[],2);
[maxSpec,fMaxSpec] = max(maxSpec);

% find spectral range
specFMax_loc = zeros(1,2);
specFMax = specAvg(:,tMaxSpec)';
[pks,~,~,~,wx] = findpeaks(specFMax);
pks_idx = find(pks == maxSpec);
if ~isempty(pks_idx)
    specFMax_loc(1,:) = wx(pks_idx,:);
    specFMax_loc(1) = max(specFMax_loc(1),find(fSpec>(fSpec(fMaxSpec)-fWin/2),1));
    specFMax_loc(2) = min(specFMax_loc(2),find(fSpec<(fSpec(fMaxSpec)+fWin/2),1,'last'));
else
    error('Could not find spectral peak, please check.')
end
% convert to frequency in Hz
specRange = fSpec(round(specFMax_loc));

% find temporal range
% specTMax_PkLoc = zeros(1,1);
specTMax_loc = zeros(1,2);
specTMax = specAvg(fMaxSpec,:);
[pks,locs,~,~,wx] = findpeaks(specTMax);
pks_idx = find(pks == maxSpec);
if ~isempty(pks_idx)
    specTMax_PkLoc = locs(pks_idx);
    specTMax_loc(1,:) = wx(pks_idx,:);
    % get left drop bound
    specTMax_loc_l = tMaxSpec+1 - find(specTMax(tMaxSpec:-1:1)<0.71*max(specTMax),1);
    specTMax_loc_r = tMaxSpec-1 + find(specTMax(tMaxSpec:end)<0.71*max(specTMax),1);
    specTMax_loc(1) = max(specTMax_loc(1),specTMax_loc_l);
    specTMax_loc(2) = min(specTMax_loc(2),specTMax_loc_r);
else
    error('Could not find temporal peak, please check.')
end
% adjust with half window size
tempRange = round(specTMax_loc*(nSamp/length(tSpec))*nShiftSTFT + nWinSTFT/2);
tempPeak  = round(specTMax_PkLoc*(nSamp/length(tSpec))*nShiftSTFT + nWinSTFT/2);

end





