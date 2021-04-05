%% Time-frequency analysis
% Envelope baseline-crossing sieving
% This function defines a baseline for the amplidue of the envelope to
% sieve the dataset using this baseline.
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


function [sig,sig2,sig_idx] = jc_tfa_envelopeSieve(sig,thrd,count,rule,nWin,sig2,sig_idx) %#ok<INUSL>

if ~iscell(sig)
    sig = {sig};
end
if ~exist('sig2','var')
    sig2 = [];
end
if length(nWin) > 1
    nWin_bkg = nWin(1); % background window
    nWin_tgt = nWin(2); % target window
else
    nWin_bkg = nWin;
    nWin_tgt = nWin;
end

nChan = find(cellfun(@isempty,sig)==0)';
for iChan = nChan
    % tmp signal matix
    sig_tmp = sig{iChan};
    nLen = size(sig_tmp,2);

    % calculate envelope
    % analytical signal use hilbert transform , then extract envelope
    [sigEvlpUB,sigEvlpLB] = envelope(sig_tmp');
    sigEvlpUB = sigEvlpUB'; sigEvlpLB = sigEvlpLB';
    % set threshold
    idxThrd = [1:nWin_bkg,(nLen-nWin_bkg+1):nLen];
    thrd_ub = thrd*median(sigEvlpUB(:,idxThrd),2);
    thrd_lb = thrd*median(sigEvlpLB(:,idxThrd),2);
    
    % signal index for zero-crossing examine
    if exist('nWin','var')
        idx = round(nLen/2+(-nWin_tgt/2:nWin_tgt/2));
    else
        idx = 1:nLen;
    end
    % zero-crossing
    %%% positive part
    zcCount_p = jc_tfa_zeroxDetector(sig_tmp(:,idx),thrd_ub);
    %%% negative part
    zcCount_n = jc_tfa_zeroxDetector(-sig_tmp(:,idx),-thrd_lb);
    %%% estimate the overall crossing count
    zcCount = mean([zcCount_p,zcCount_n],2); %#ok<NASGU>
    
    % count zero-crossing
    zcCount_idx = eval(['zcCount',rule,'count']);
    sig{iChan} = sig_tmp(zcCount_idx,:);
    if ~isempty(sig2)
        sig2{iChan} = sig2{iChan}(zcCount_idx,:);
        sig_idx{iChan} = sig_idx{iChan}(zcCount_idx,:);
    end
end

end

