%% Basic function
% This function finds the folder and file name according to the input
% string of file path.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [folderFile,fileName] = jc_bs_findFileFolder(filePath)

fileSepIdx = strfind(filePath,'\');
if isempty(fileSepIdx)
    fileSepIdx = strfind(filePath,'/');
    if isempty(fileSepIdx)
        folderFile = './';
        fileName = filePath;
        return
    end
end
fileSepIdx = fileSepIdx(end);

folderFile = filePath(1:fileSepIdx);
fileName = filePath((fileSepIdx+1):end);

end