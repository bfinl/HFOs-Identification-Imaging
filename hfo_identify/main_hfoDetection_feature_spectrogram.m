%% HFO detection
% Feature extraction sub-function
%
% Spectrogram features extrcted from spectrogram images
%   (1) image entropy
%   (2) third-order moment of statistical histogram
%   (3) directionality
%   (4) block-wise power spectrum density (PSD)
%
% The feature design was mainly based on image processing and inspired by
% several works as follows.
% Shi, X., et al. (2015). "Textural feature extraction based on
% time–frequency spectrograms of humans and vehicles." IET Radar, Sonar &
% Navigation 9(9): 1251-1259.
%
% Tzallas, A. T., et al. (2009). "Epileptic seizure detection in EEGs using
% time-frequency analysis." IEEE Trans Inf Technol Biomed 13(5): 703-710. 
%
% To run this function, the Wavelet and Imaging Processing Toolbox in
% Matlab is required for calculation some of the spectrogram features.
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
% 2018.06.05
% Zhengxiang Cai
% Initial implementation.
%
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function specFeature = main_hfoDetection_feature_spectrogram(sig,fs,nWin)

% extract signal windows, if needed
if nargin < 3
    sigWin = sig;
else
    sigWin = jc_io_extractEpoch(sig,{round(length(sig)/2)},nWin);
    sigWin = sigWin{1};
end

% compute spectrogram using wavelet transform
[~,specAmp,t,f,~] = jc_tfa_WaveletSpectrum(sigWin,fs,'morse',0);

% pre-processing of TF spectrogram
% compute energy sum over time
specEnrg = sum(specAmp,2);

% extract features from spectrogram
specFeature = zeros(1,18);
% (1) image entropy
specEnrg_p = jc_tfa_normPDF(specEnrg);
specFeature(1) = wentropy(specEnrg_p,'shannon');

% (2) third-order moment
[pixelCount, grayLevels] = imhist(specAmp);
pixelCount_mean = (pixelCount'*grayLevels)/sum(pixelCount);
specFeature(2) = sum( (grayLevels-pixelCount_mean).^3.*pixelCount );

% (3) directionality
% compute gradients
[Gx,Gy] = imgradientxy(specAmp);
% amplitude of gradients
specGAmp = sqrt(Gx.^2+Gy.^2)./2;
% find mask of gradient-map using threshold
specG_mask = specGAmp > mean(mean(specGAmp));
% compute gradient angle
specTheta = atan(Gy./Gx);
% zero out nan angles
specG_mask(isnan(specTheta)) = 0;
specTheta(isnan(specTheta)) = 0;
specTheta_effect = specTheta(specG_mask>0);
% compute histogram
[pixelCount, grayLevels] = imhist(specTheta_effect);
pixelCount_p = jc_tfa_normPDF(pixelCount);
% compute directionality
[~,pixelCount_peak_loc] = findpeaks(pixelCount,...
    'MinPeakDistance',10);
pixelCount_peak_thrd = median(pixelCount(pixelCount_peak_loc));
[~,pixelCount_peak_loc] = findpeaks(pixelCount,...
    'MinPeakDistance',10,'minpeakheight',pixelCount_peak_thrd);
nPeak = length(pixelCount_peak_loc); nPixel = length(pixelCount);
specFeature(3) = sum(sum(((repmat(grayLevels,1,nPeak) -...
    repmat(grayLevels(pixelCount_peak_loc),1,nPixel)').^2).*repmat(pixelCount_p,1,nPeak)));
if specFeature(3)==0
    keyboard
end

% (4) block-wise power spectrum density (PSD)
% the 2d spectrogram is divided into 12 blocks
% 3 equal interval in time
% 4 bands in frequency: 0-15 Hz, 15-30 Hz, 30-80 Hz, 80-150 Hz, 150-240 Hz
% all energy bins are normalized by the total energy of tghe whole
% spectrogram, as a ratio

% total spectrum energy
specEnrg_sum = sum(specEnrg);
% define time and frequency bins
binTime = [0, max(t)/3, 2*max(t)/3, max(t)];
binFreq = [0, 15, 30, 80, 150, 240];
% define block-wise energy bins
binEnrg = zeros(1,(length(binTime)-1)*(length(binFreq)-1));
for iBinFreq = 1:length(binFreq)-1
    idxFreq = f>binFreq(iBinFreq) & f<binFreq(iBinFreq+1);
    for iBinTime = 1:length(binTime)-1
        idxTime = t>binTime(iBinTime) & t<binTime(iBinTime+1);
        binEnrg( (iBinFreq-1)*(length(binTime)-1) + iBinTime ) = ...
            sum(sum(specAmp(idxFreq,idxTime)))/specEnrg_sum;
    end
end
specFeature(4:end) = binEnrg;

end










