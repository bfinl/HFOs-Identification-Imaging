%% setup path and folder
% The code and data dependency is added for the framework.
% Toolboxes used for pre-processing and topo-map plot are listerd.
%
%%% EEGLab, which can be downloaded and deploy from
%%% https://sccn.ucsd.edu/eeglab/index.php.
%%% The version of EEGLAB tested is version 14.1.1b (recommended) - also
%%% tested on eeglab v13.6.5b and above.
%
%%% jLab, for time-frequency visualization, which can be downloaded and
%%% deployed from https://github.com/jonathanlilly/jLab.
%%% The version of jLab tested is version 1.6.2.
% 
%%% regtools, for imaging regularization, which can be downloaded and
%%% deployed from http://www2.compute.dtu.dk/~pcha/Regutools/.
%%% The version of Regularization Tools tested is version 4.1.
%
%%% Other third-party codes
%%% (1) SSD, Spatio-Spectral Decomposition, is used for extraction and
%%%     denoising of HFO activities, which can be downloaded from
%%%     https://github.com/svendaehne/matlab_SSD/.
%%% (2) Optimal kmeans, kmeansElbow, is used for select optimal kmeans
%%%     clusters for HFOs identification/extraction.
%%%     https://it.mathworks.com/matlabcentral/fileexchange/65823-kmeans_opt
%%% (3) findpeaks, built-in function in Matlab, with minor alterations for
%%%     extra outputs, is used for locating the peaks in a signal.
%%% (4) tightSubplot, is used for setting axes with adjustable margins and
%%%     gaps. This function was implemented by Pekka Kumpulainen and is
%%%     available from the following link:
%%%     https://www.mathworks.com/matlabcentral/fileexchange/27991-
%%%     tight_subplot-nh-nw-gap-marg_h-marg_w.
%
% Please put EEGLab, jLab, and regtools into "toolbox" folder and please
% find the subfolders as "eeglab", "jLab", and "regtools". Please also
% download the other third-party codes with the provided links into the
% "thirdParty" folder. The functions kmeansElbow, findpeaks, and
% tightSubplot with minor modifications has been provided in the
% "thirdParty" folder already.
%
% The user can also put the toolboxes or codes at other places, as long as
% the codes could be searched and called by Matlab.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function defaultBaseDataPath = setupPathAndFolder()
% restoredefaultpath;
% add dependencies to path
mainPath = mfilename('fullpath');
[mainPath,~,~] = fileparts(mainPath);
addpath(genpath(mainPath));

toolboxPath = fullfile(mainPath,'toolbox');
% jLab, for wavelet representation
addpath(fullfile(toolboxPath,'jLab','jwavelet'));
% eeglab, for EEG signal processing
addpath(fullfile(toolboxPath,'eeglab')); eeglab(); close all;
% regtools, for ESI regularization
addpath(fullfile(toolboxPath,'regtools'));
% other third-party codes
addpath(fullfile(toolboxPath,'thirdParty'));

% setup working forders
defaultBaseDataPath = fullfile(mainPath,'data');

