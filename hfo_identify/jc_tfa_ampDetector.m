%% Time-frequency analysis
% Use defined threshold to find time points that exceed the threshold.
% A defined window size is used to combine events occurring together within
% the samples of this window.
%
%%% License
% We provide our code and data under a CC-BY-NC-SA-4.0 license, "As is" and
% without any guarantee to the scientific community for academic and
% research purposes primarily, not commercial use.
% 
% You should have received a copy of the CC-BY-NC-SA-4.0 license along with
% this program. If not, see https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function sigEpochIdx = jc_tfa_ampDetector(sig,thrd,nWin)

[nChan,nSig] = size(sig);
nThrd = length(thrd);
if nChan ~= nThrd
    error('Channel size is not equal to threshold size.');
end

% process data in row (channel) series
sigEpochIdx = cell(nChan,1);
for iChan = 1:nChan
    % find all possible locations
    sig_tmp = abs(sig(iChan,:)-mean(sig(iChan,:)));
    sigEpochIdx_tmp = find(sig_tmp > thrd(iChan));
    sigEpochIdx_tmp(sigEpochIdx_tmp < nWin) = [];
    iEvent = 1;
    while ~isempty(sigEpochIdx_tmp)
        % use all possible locations as a pool
        % find maximum of the whole, and do segment within a window
        if sigEpochIdx_tmp(1)+nWin <= nSig
            [~,sigEpochIdx_tmp2] = max(abs(sig(iChan,sigEpochIdx_tmp)));
        else
            break
        end
        sigEpochIdx{iChan}(iEvent,1) = sigEpochIdx_tmp(sigEpochIdx_tmp2);
        % delete the other smaller activities within the window
        sigEpochIdx_tmp( (sigEpochIdx_tmp>=(sigEpochIdx{iChan}(iEvent,1)-nWin))...
            & (sigEpochIdx_tmp<=(sigEpochIdx{iChan}(iEvent,1)+nWin)) ) = [];
        % find next event
        iEvent = iEvent+1;
    end
end

