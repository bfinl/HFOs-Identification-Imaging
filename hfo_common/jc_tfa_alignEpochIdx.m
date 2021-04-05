%% Time-frequency analysis
% This function aligns epochs of data according to alignment index.
% Subfuntion for alignEpoch
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function sigAlignedSeg = jc_tfa_alignEpochIdx(sigOrigSeg,sig_peak_idx,...
    nChan,nSamp,nSeg,alignMethod)

if nargin < 6
    alignMethod = 'MinTruncate';
end

alignMethod = lower(alignMethod);
switch alignMethod
    case {'mintruncate','min','minimumtruncate'}
        % find minimum and maximum alignment points
        sigMGFP_peak_idx_min = min(sig_peak_idx);
        sigMGFP_peak_idx_max = max(sig_peak_idx);
        % maximum offset
        sigMGFP_peak_idx_maxOffset = sigMGFP_peak_idx_max - sigMGFP_peak_idx_min;

        % align
        sigAlignedSeg = zeros(nChan,nSamp-sigMGFP_peak_idx_maxOffset,nSeg);
        for iSeg = 1:nSeg
            sigMGFP_offset = sig_peak_idx(iSeg) - sigMGFP_peak_idx_min;
            sigAlignedSeg(:,:,iSeg) = sigOrigSeg(:,sigMGFP_offset+(1:end-sigMGFP_peak_idx_maxOffset),iSeg);
        end
    case {'center','cen','central'}
        % find minimum half length
        [sigMinHalfLength,sigMinHalfLengthIdx] = min([sig_peak_idx(:); nSamp-sig_peak_idx(:)]);
        if sigMinHalfLengthIdx <= nSeg % first half
            sigMinHalfLength = sigMinHalfLength - 1;
        end
        sigMinFullLength = 2*(sigMinHalfLength)+1;
        % align
        sigAlignedSeg = zeros(nChan,sigMinFullLength,nSeg);
        for iSeg = 1:nSeg
            sigMGFP_offset = sig_peak_idx(iSeg) - sigMinHalfLength;
            if nSeg > 1
                sigAlignedSeg(:,:,iSeg) = sigOrigSeg(:,sigMGFP_offset+(0:sigMinFullLength-1),iSeg);
            else
                sigAlignedSeg(:,:) = sigOrigSeg(:,sigMGFP_offset+(0:sigMinFullLength-1));
            end
        end
    case {'zeropad','zp'}
        % aligh to the center
        sigAlignedSeg_tmp = jc_tfa_alignEpochIdx(sigOrigSeg,sig_peak_idx,...
            nChan,nSamp,nSeg,'cen');
        % zero-pad to original length
        nLenZp = nSamp-length(sigAlignedSeg_tmp(1,:,1));
        nLenZp1 = round(nLenZp/2);
        nLenZp2 = floor(nLenZp/2);
        sigAlignedSeg = zeros(nChan,nSamp,nSeg);
        % zero-padding using signal bounds
        for iSeg = 1:nSeg
            sigAlignedSeg(:,:,iSeg) = [repmat(sigAlignedSeg_tmp(:,1,iSeg),1,nLenZp1),...
                sigAlignedSeg_tmp(:,:,iSeg),repmat(sigAlignedSeg_tmp(:,end,iSeg),1,nLenZp2)];
        end
end


