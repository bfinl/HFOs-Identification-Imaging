%% Visualization
% This function plots topo-map for input signal.
% The implementation is based on EEGLab.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [fig,elec_loc] = plotTopoMap(sig,fig,pathChan,badChannel,fNormSig)

% load EEGLab, if needed
if ~exist('topoplot','file') || ~exist('readlocs','file')
    eeglab();
end

% initialize input variables
if nargin<2 || isempty(fig)
    fig = figure;
end
if nargin < 3
    pathChan = [];
end
if nargin < 4
    badChannel = [];
end
if nargin < 5 || isempty(fNormSig)
    fNormSig = true;
end

% load channel locations
elec_loc = readlocs(pathChan);
elec_loc(badChannel)=[];

% setup EEG data
if size(sig,2)>1
    [~,sigMax_idx] = max(rms(sig,1));
    sig = sig(:,sigMax_idx);
end
if fNormSig
    % map to [-1,1]
    sig = 2*(sig - mean(sig))./(max(sig)-min(sig));
end

topoplot(sig,elec_loc,'electrodes','off','maplimits','absmax','style','map');

end

