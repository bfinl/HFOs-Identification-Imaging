%% Input/Output function
% Extract epochs from data segments according to input index information.
%
% Input:
%   sig: the signal to be segmented, which can be multi-channeled
%   sigEpochIdx: the index for the center sample of the segmentation
%   nWin: the window size of the epoch to be segmented out
% Output:
%   sigEpoch: the signal cell that has been generated 
%   sigEpochIdx: the updated version of the index input, for those epochs
%   have been extracted, the index would be the same, while for those are
%   not, the index would be null.
% Example:
% --------|--------*---------|-----|--------*---------|----------
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [sigEpoch,sigEpochIdx] = jc_io_extractEpoch(sig,sigEpochIdx,nWin)

% check if signal epoch index has the same size as signal channel (row)
[nChan,nLen] = size(sig);
nIdx = length(sigEpochIdx);
% some legacy definition
if ~iscell(sigEpochIdx)
    sigEpochIdx = sigEpochIdx(:)';
    sigEpochIdx = num2cell(sigEpochIdx);
end
if nChan ~= nIdx
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
        sigEpoch{iChan} = zeros(nEpoch,2*nWin+1);
        for iEpoch = 1:nEpoch
            idx = (-nWin:nWin)+sigEpochIdx_tmp(iEpoch);
            sigEpoch{iChan}(iEpoch,:) = sig(iChan,idx);
        end
    end
    sigEpochIdx{iChan} = sigEpochIdx_tmp;
end
