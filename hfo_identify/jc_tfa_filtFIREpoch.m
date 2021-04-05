%% Time-frequency analysis
% This function performs filtering on multichannel dataset.
%
% Input:
%   sigOrigEpoch: the unfiltered signal epoch
%   fs: sampling frequency, in Hz
%   freqLb,freqUb: lower and upper bound of filtering band
%   order: order of filter design
%   filtType: bandpass or highpass
%   varargin: specified inputs for filter design
% Output:
%   sigFiltEpoch: the filtered signal epoch
%   bpFilt: designed filter
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [sigFiltEpoch,bpFilt] = jc_tfa_filtFIREpoch(sigOrigEpoch,fs,freqLb,freqUb,...
    order,filtType,varargin)


switch filtType
    case {'bandpass','band'}
        if isempty(varargin)
            bpFilt = designfilt('bandpassfir','FilterOrder',order, ...
                     'CutoffFrequency1',freqLb,'CutoffFrequency2',freqUb, ...
                     'SampleRate',fs);
        else
            bpFilt = designfilt('bandpassfir',varargin{:},...
                     'SampleRate',fs);
        end
    case {'highpass','high'}
        if isempty(varargin)
            bpFilt = designfilt('highpassfir','StopbandFrequency',freqLb, ...
                     'PassbandFrequency',freqUb,'PassbandRipple',0.5, ...
                     'StopbandAttenuation',order,'DesignMethod','kaiserwin',...
                     'SampleRate',fs);
        else
            bpFilt = designfilt('highpassfir',varargin{:},...
                     'SampleRate',fs);
        end
end

sigFiltEpoch = zeros(size(sigOrigEpoch));
[~,~,nEpoch] = size(sigOrigEpoch);
for iEpoch = 1:nEpoch
    sigOrigEpoch_tmp = sigOrigEpoch(:,:,iEpoch);
    sigFiltEpoch(:,:,iEpoch) = filtfilt(bpFilt,sigOrigEpoch_tmp')';
end

end


