%% Print function
% Print string in a block.
%
%--------------------------------------------------------------------
% Chris Cline
% 2017.06.01
% Initial implementation.
%
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function varargout = jc_print_block(varargin) %#ok<STOUT>

	global sayDateFormat;
    if isempty(sayDateFormat)
		sayDateFormat = 'HH:MM:ss';
    end
	
    % print time
	if verLessThan('matlab','8.4')
        fprintf('---------------------------------------------------------\n');
		fprintf('%s ',datestr(now,13));
    else
        fprintf('---------------------------------------------------------\n');
		fprintf('%s ',datestr(datetime,sayDateFormat));
	end

	% print input string
    for i = 1
        fprintf(' |');
    end
	fprintf('- ');
	fprintf(varargin{:});
	fprintf('\n');
    fprintf('---------------------------------------------------------\n');
    
end

