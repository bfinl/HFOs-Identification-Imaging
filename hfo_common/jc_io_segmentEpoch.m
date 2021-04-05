%% Input/Output function
% Segment concatenated signal dataset
%
% This function is based on reshape function, and it could deal with
% dataset of permuted dimensions.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [sigSegment,nEpoch] = jc_io_segmentEpoch(sigEpoch,nChannel,nEpochLen)

% check data size
[mSig,nSig] = size(sigEpoch);
if nChannel ~= mSig
    if nSig == nChannel
        sigEpoch = reshape(sigEpoch,mSig,[]);
        sigEpoch = sigEpoch.';
        sigEpoch = reshape(sigEpoch,nChan,mSig,[]);
    else
        error('Signal size is INCORRECT.');
    end
end

% segment signal epoch
sigSegment = reshape(sigEpoch,nChannel,nEpochLen,[]);
nEpoch = size(sigSegment,3);
