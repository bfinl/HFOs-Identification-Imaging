%% HFO detection
% This function prepares data file.
% (1) load dataset with basic pre-processing
% (2) average spikes
% (3) setup electrode file
%
% Note that the dataset contains raw spike epochs exported from Curry 7.
% The epochs were filtered above 1Hz and notch filtered at 60Hz with
% harmonics removed. The spike epochs were averaged and saved for later
% processing and source imaging. The electrode file contains the digitizer
% of the surface channels, which defines the sequence, the name, and the 3D
% location in the scalp, which will be used for visualization purposes.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


% data file are exported from Curry 7 and converted to .mat file
load(fullfile(defaultBaseDataPath,'raw','Pt_sample_raw_data.mat'));
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


