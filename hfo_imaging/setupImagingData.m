%% HFO imaging
% setup imaging data, extracted signals, and parameters
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


% prepare for source estimation
% regularization
paramMethod = 'lcurve';
% electrode file path
pathChannel = fullfile(pathInput,'electrodes.xyz');

% data file setup
% get folder structure
fileInputDataInfo_set = struct2cell(dir);
fileInputDataName_set = fileInputDataInfo_set(1,:);
% get hfo data
fileHFOLocation = fileInputDataName_set(...
    strcmpi(fileInputDataName_set,'HFOResults.mat'));
% get spike data
fileSpikeLocation = fileInputDataName_set(...
    strcmpi(fileInputDataName_set,'SpikeAvg.mat'));
% get lead field
fileLeadFieldLocation = fullfile(pathRaw,'leadField.mat');

% load hfo data
load(fileHFOLocation{1},'hfoMap','hfoSig','chan_list');
sig_hfo = hfoMap;
% load spike data
load(fileSpikeLocation{1},'sigEpochSpikeAvg');
sig_spk = sigEpochSpikeAvg{1};
% centralize data
sig_hfo = jc_sl_carSensor(sig_hfo);
sig_spk = jc_sl_carSensor(sig_spk);
% signal content location
sigSpot_hfo = [1 1 1];
sigSpot_spk = round([1 size(sig_spk,2) size(sig_spk,2)/2]);
% empirical noise data
noiseCov_hfo = noiseEstimation(hfoSig,0.95);
noiseCov_spk = noiseEstimation(sig_spk,0.95);

% check signal morphology
% figure;
% subplot(1,2,1), jc_vsl_topoMap(sig,gcf,pathChannel,[]);
% subplot(1,2,2), plot(sig');
% xtickIdx = find(sig > median(sig));
% set(gca,'xtick',xtickIdx,'xticklabel',chan_list(xtickIdx));


