%% Basic function
% This function is to compose event pool for multichannel data.
% It is simply implemented to convert cell epochs to regular data matrix.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function eventPool = jc_bs_composeEventPoolWhole(eventPool_set)

nCell = find(~cellfun(@isempty,eventPool_set));

eventPool = [];
for iCell = nCell'
    eventPool = cat(3,eventPool,eventPool_set{iCell});
end

