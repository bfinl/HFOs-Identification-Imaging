%% HFO detection
% This function prepares data file from simulation.
% (1) load dataset with basic pre-processing
% (2) setup electrode file (containing channel information)
%
% Note that the dataset contains raw spike epochs simulated from realistic
% clinical dataset. The attached simulation data consists of epoched spikes
% with additive HFOs and artifacts, embeded in white or realistic EEG background
% noise of SNR at 5~20dB. The noise data is modeled from white Gaussian
% process and real clinical recordings using autoregressive model.
% The electrode file contains the digitizer of the surface channels, which
% defines the sequence, the name, and the 3D location in the scalp, which
% will be used for visualization purposes.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


% data file with simulated HFOs and artifacts under various noise levels
%%% white Gaussian noise
% load(fullfile(defaultBaseDataPath,'raw','Pt_simulation_10dB_white.mat'));
%%% realistic background noise
load(fullfile(defaultBaseDataPath,'raw','Pt_simulation_10dB_real.mat'));
sigRaw = sigEpochSpike_cat;

% setup channel file
copyfile(fullfile(defaultBaseDataPath,'raw','electrodes.xyz'),...
    pathOutput);
pathChannel = fullfile(pathOutput,'electrodes.xyz');

% save and display dataset information
cd(pathOutput);
save(fullfile(pathOutput,'1_HFOBasic'));
fprintf('Dataset imported.\n');
fprintf('Total %d epochs counted.\n',length(sigEpochInfo));


