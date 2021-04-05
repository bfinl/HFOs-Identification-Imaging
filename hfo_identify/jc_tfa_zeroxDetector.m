%% Time-frequency analysis
% Zero-crossing detector
% This fucntion operates on vector or martix containing several dimensions
% of signals in rows. Then, according to a defined baseline, or zero by
% default, the function finds all crossing points along this baseline.
% The function output is the counts of baseline-crossing.
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


function zcCount = jc_tfa_zeroxDetector(sig,thrd)

if nargin < 2
    thrd = 0;
end

[nChan,~] = size(sig);
% set baseline-threshold to every channel
if isscalar(thrd)
    thrd = repmat(thrd,nChan,1);
end

zcCount = zeros(nChan,1);
for iChan = 1:nChan
    % re-centralize signal to threshold
    sig_tmp = sig(iChan,:)-thrd(iChan);
    % remove zero, just in case there would be fake crossing
    % in this case, the rule would be more strict
    % that is, points just touching zero does not count
    sig_tmp(sig_tmp==0) = [];

    % detect zero-crossing
    sigSgn = sign(sig_tmp);
    sigSgnDiff = diff(sigSgn);
    zcCount(iChan) = length(nonzeros(sigSgnDiff));
end

end

