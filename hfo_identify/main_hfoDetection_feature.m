%% HFO detection
% Feature extraction main function
%
% Two set of features are extracted for clustering
% 1) time-frequency domain features
% 2) spectrogram features
% Features extratced from both raw and cleaned signals
% In total, 30d x 1, per event sample
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


fprintf('Extracting features...Please wait...\n');

fprintf('Extracting features...Please wait...TFD\n');
% time-frequency domain features
%   (1) time domain:
%       mean, var, skewness
%   (2) frequency domain:
%       median frequency, spectral centroid, power band ratio
% Features extratced from both raw and cleaned signals
% 15d x 1, named feature_1
nWin = round(100/1000*fs);
nEvent = size(sigEventPoolRaw,1);
feature_1 = zeros(nEvent,15);
for iEvent = 1:nEvent
    sigEventPoolRaw_tmp = sigEventPoolRaw(iEvent,:);
    sigEventPoolFilt_tmp = sigEventPoolFilt(iEvent,:);
    feature_1(iEvent,:) = main_hfoDetection_feature_tfd(...
        sigEventPoolRaw_tmp,sigEventPoolFilt_tmp,fs,nWin);
end

fprintf('Extracting features...Please wait...Spectrogram\n');
% spectrogram features extrcted from spectrogram images
%   (1) image entropy
%   (2) third-order moment of statistical histogram
%   (3) directionality
%   (4) block-wise power spectrum density (PSD)
% Features extratced from both raw and cleaned signals
% 15d x 1, named feature_2
nWin = min(round(200/1000*fs),floor(size(sigEventPoolRaw,2)/2));
nEvent = size(sigEventPoolRaw,1);
feature_tmp1 = zeros(nEvent,3+15);
feature_tmp2 = zeros(nEvent,3+15);
for iEvent = 1:nEvent
    
    sigEventPoolRaw_tmp  = sigEventPoolRaw(iEvent,:);
    sigEventPoolFilt_tmp = sigEventPoolFilt(iEvent,:);
    
    feature_tmp1(iEvent,:) = main_hfoDetection_feature_spectrogram(sigEventPoolRaw_tmp,fs,nWin);   
    feature_tmp2(iEvent,:) = main_hfoDetection_feature_spectrogram(sigEventPoolFilt_tmp,fs,nWin);
    
end
feature_2 = [feature_tmp1(:,4:(4+3*2-1)) feature_tmp2(:,(4+3*2):end)];

% save all results
% save(fullfile(pathOutput,'4_HFOFeatures.mat'));
jc_print_block('Feature extraction done.');



