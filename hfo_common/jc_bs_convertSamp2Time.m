%% Basic functions
% This function converts index (sample) points to time points.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function sigTimeIdx = jc_bs_convertSamp2Time(sigSampIdx,fs)

sigTimeIdx = sigSampIdx./fs.*1000;

