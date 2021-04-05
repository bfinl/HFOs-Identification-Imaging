%% Input/Output function
% This function reads channel list from electrode information file.
% Read channel file and generate channel list as output.
%
% The toolbox EEGLab is needed for reading digitizer, which can be
% downloaded and deployed from https://sccn.ucsd.edu/eeglab/index.php.
% The version of EEGLAB tested is version 14.1.1b (recommended) - also
% tested on eeglab v13.6.5b and above.
% Please also refer to the function "setupPathAndFolder.m".
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [chanLbl,chanLoc] = jc_io_readChannelList(pathChan)

if ~exist('readlocs.m','file')
    % readlocs.m is a function to load digitizer, built-in from eeglab
    eeglab(); % initialize eeglab
end

chan = readlocs(pathChan);
chan = struct2cell(chan);
chanLbl = squeeze(chan(4,1,:));

chanLoc = squeeze(chan(1:3,1,:));
chanLoc = cell2mat(chanLoc)';
chanLoc = [chanLoc(:,2),chanLoc(:,1),chanLoc(:,3)];

