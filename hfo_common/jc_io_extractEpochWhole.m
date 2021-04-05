%% Input/Output function
% This function extracts multichannel epochs from data segments according
% to input index list of events.
%
% Example:
% --------|--------*---------|-----|--------*---------|----------
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function sigEpoch = jc_io_extractEpochWhole(sig,sigEpochIdx,nWin)

% check if signal epoch index has the same size as signal channel (row)
[nChan,nLen] = size(sig);
nIdx = length(sigEpochIdx);
if nChan~=nIdx
    error('Signal channel does not match epoch index size.')
end

% extract epochs
sigEpoch = cell(nChan,1);
for iChan = 1:nChan
    % check if epoch index are greater than window size
    % if not, delete the index
    sigEpochIdx_tmp = sigEpochIdx{iChan};
    sigEpochIdx_tmp(sigEpochIdx_tmp<=nWin) = [];
    sigEpochIdx_tmp(sigEpochIdx_tmp>(nLen-nWin)) = [];
    if ~isempty(sigEpochIdx_tmp)
        nEpoch = length(sigEpochIdx_tmp);
        sigEpoch{iChan} = zeros(nChan,2*nWin+1,nEpoch);
        for iEpoch = 1:nEpoch
            idx = (-nWin:nWin)+sigEpochIdx_tmp(iEpoch);
            sigEpoch{iChan}(:,:,iEpoch) = sig(:,idx);
        end
    end
end
