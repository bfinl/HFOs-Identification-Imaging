%% Time-frequency analysis
% This function uses continuous wavelet transform for time-frequency
% analysis on bio-signals.
% The main implementation is in JLab toolbox.
%
% Inputs:
%   sig: signal in 1D to be decomposed with wavelet decomposition
%   fs: sampling frequency
%   wname: wavelet kernel name
%   flagPlot: logical variable for plotting option
%   filtRange: optional, the frequency range used to extract/reconstruct a
%   specific band signal
%   filtOpt: wavelet density, or overlapping size during WD, default 10
%   sigType: output signal type, 0: real or 1: complex
%   fig: figure handle if any
%   ax: axes handle if any
% Outputs:
%   rSignal: output reconstructed signal, complex or real
%   h: figure handle of WD plot
%   wtPlot: WD matrix in real, used for plotting
%   t: time scales
%   freq: frequency scales
%   wtRaw: WD matrix in complex, for reference
%
% The toolbox jLab is needed for wavlet transformation, which can be
% downloaded and deployed from http://www.jmlilly.net/doc/jLab.html.
% The version of jLab tested is version 1.6.2.
% Please also refer to the function "setupPathAndFolder.m".
%--------------------------------------------------------------------
% 2017.02.11
% Zhengxiang Cai
% Initial implementation.
%
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function  [h,wtPlot,t,freq,wtRaw] = jc_tfa_WaveletSpectrum(...
    sig,fs,wname,flagPlot,filtOpt,fig,ax,cxFactor)

% add jLab path
if ~exist('morsespace','file')
    jlab_addpath;
end
        
% initial input parameters
% wavelet name
if nargin < 3
    wname = 'morse';
end
% sample frequency
if nargin < 2
    fs = 1000;
end
% plot figure option
if nargin < 4
    flagPlot = 1;
end
% filtering Option
if nargin < 5
    filtOpt = 10;
end
% output figure
if nargin < 6 && flagPlot
    fig = figure;
    plot(1);
end
if nargin < 7 && flagPlot
    ax = fig.CurrentAxes;
end
if nargin < 8 && flagPlot
    cxFactor = 0.01;
end

% sample time
t = 0:1/fs:(length(sig)-1)*1/fs;
sig = sig(:);

% Compute wt and plot using jLab
switch wname
    case {'morse'}
        gamma = 5; beta = 2;
    case {'morl'}
        gamma = 3; beta = 2;
end

psi = morsespace(gamma,beta,2*pi,2*pi/1000,filtOpt);
freq = fs/2/2/pi*psi;
wt = wavetrans(sig,{gamma,beta,psi,'bandpass'});
wtRaw = wt';
wtPlot = abs(wt');
if flagPlot
    h = fig; % axes(ax);
    pcolor(ax,t,freq,wtPlot), shading(ax,'interp');
    ylog(fig,ax), set(ax,'tickdir','in');
    set(ax,'ylim',[1 max(freq)]);
    xlabel('Time (sec)');
    ylabel('Pseudo-frequency (Hz)','interpreter','tex');
    title(['Wavelet Spectrum with "',wname,'"']);
    colorbar;
    caxis([quantile(wtPlot(:),cxFactor),quantile(wtPlot(:),1-cxFactor)]);
else
    h = [];
end

end

