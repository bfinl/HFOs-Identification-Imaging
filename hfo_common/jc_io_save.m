% Input/Output function
% This function saves dataset with current time tag.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function jc_io_save(varargin)

if isnumeric(varargin{end})
    flagBase = varargin{end};
    varargin(end) = [];
else
    flagBase = 1;
end    

% string expression of saving command
expr = 'save(';
for i = 1:length(varargin)
    expr = [expr,'''',varargin{i},'''']; %#ok<*AGROW>
    if i < length(varargin)
        expr = [expr,','];
    end
end
expr2 = [expr,',''-append'');'];
expr1 = [expr,');'];

% save variables in corresponding workspace 
if flagBase
    % base workspace
    evalin('base','jc_saveTime = datetime;');
    try
        evalin('base',expr2);
    catch
        evalin('base',expr1);
    end
    evalin('base',sprintf('save(''%s'',''jc_saveTime'',''-append'');',varargin{1}));
else
    % calling function workspace
    evalin('caller','jc_saveTime = datetime;');
    try
        evalin('caller',expr2);
    catch
        evalin('caller',expr1);
    end
    evalin('caller',sprintf('save(''%s'',''jc_saveTime'',''-append'');',varargin{1}));
end

