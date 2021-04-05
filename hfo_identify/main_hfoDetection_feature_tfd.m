%% HFO detection
% Feature extraction sub-function
% Regular time-frequency domain features, defined at local events
%
% Implemented measure:
%   (1) time domain
%       mean, var, skewness
%   (2) frequency domain:
%       median frequency, spectral centroid, power band ratio
% Features extratced from both raw and cleaned signals
% In total, 15d x 1, per event sample
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


function TFD = main_hfoDetection_feature_tfd(sigRaw,sigFil,fs,nWin)

% extract wimt window, if needed
if nargin < 4
    sigRawWin = sigRaw;
    sigFilWin = sigFil;
else
    sigRawWin = jc_io_extractEpoch(sigRaw,{round(length(sigRaw)/2)},nWin);
    sigRawWin = sigRawWin{1};
    sigFilWin = jc_io_extractEpoch(sigFil,{round(length(sigFil)/2)},nWin);
    sigFilWin = sigFilWin{1};
end
% psd estimate
[sigRawSpect,fxx] = pwelch(sigRawWin, hann(64), 63, 2*(length(sigRaw)-1), fs); %#ok<ASGLU>
[sigFilSpect,fxx] = pwelch(sigFilWin, hann(64), 63, 2*(length(sigFil)-1), fs);

% time domain
% mean, var, skewness
tf1 = [mean(sigRawWin) mean(sigFilWin)] ;
tf2 = [var(sigRawWin) var(sigFilWin)];
tf3 = [skewness(sigRawWin,0) skewness(sigFilWin,0)];
% line length
tf4 = [sum(abs(diff(sigRawWin)))/length(sigRawWin)...
    sum(abs(diff(sigFilWin)))/length(sigFilWin)];
% global/average-local peak ratio
sigSmooth = smooth(sigRawWin,3);
[pkMax,~] = max(sigSmooth);
[pks,~] = findpeaks(sigSmooth);
tf5 = pkMax/(sum(pks(pks<pkMax))+eps)*(length(pks)-1);
sigSmooth = smooth(sigFilWin,3);
[pkMax,~] = max(sigSmooth);
[pks,~] = findpeaks(sigSmooth);
tf5 = [tf5 pkMax/sum(pks(pks<pkMax))*(length(pks)-1)];

% frequency domain
% median frequency
ff1 = [medfreq(sigRawWin) medfreq(sigFilWin)];
% spectral centroid
ff2 = sum((fxx./max(fxx)).*(sigRawSpect.^2))/sum(sigRawSpect.^2);
ff2 = [ff2 sum((fxx./max(fxx)).*(sigFilSpect.^2))/sum(sigFilSpect.^2)];
% power band ratio
lowBand_idx  = find( (fxx>1) & (fxx<80) );
highBand_idx = find( (fxx>80) & (fxx<240) );
ff3 = (sum(sigRawSpect(highBand_idx))) / (sum(sigRawSpect(lowBand_idx))); %#ok<FNDSB>

% combine all features
TFD = [tf1 tf2 tf3 tf4 tf5 ff1 ff2 ff3];



