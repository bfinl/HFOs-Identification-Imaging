%% Time-frequency analysis
% Zero-crossing sieving
% This function rejects signal epochs that not satisfy a threshold rule.
% Usually a rule is to exceed or below a defined count of zero-crossing.
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


function [sig,sig2,sig_idx] = jc_tfa_zeroxSieve(sig,thrd,count,rule,nWin,sig2,sig_idx) %#ok<INUSL>

if ~iscell(sig)
    sig = {sig};
end
if ~exist('sig2','var')
    sig2 = [];
end    

nChan = find(cellfun(@isempty,sig)==0)';
for iChan = nChan
    sig_tmp = sig{iChan};
    nLen = size(sig_tmp,2);
    if exist('nWin','var') && ~isempty(nWin)
        idx = round(nLen/2+(-nWin/2:nWin/2));
    else
        idx = 1:nLen;
    end
    % zero-crossing
    zcCount = jc_tfa_zeroxDetector(sig_tmp(:,idx),thrd); %#ok<NASGU>
    % count zero-crossing
    zcCount_idx = eval(['zcCount',rule,'count']);
    sig{iChan} = sig_tmp(zcCount_idx,:);
    if ~isempty(sig2)
        sig2{iChan} = sig2{iChan}(zcCount_idx,:);
        sig_idx{iChan} = sig_idx{iChan}(zcCount_idx,:);
    end
end

end

