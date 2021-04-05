% Time-frequency analysis
% This function normalizes the input array to a PDF function, which
% contains elements above zero and summing up to 1.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function x_norm = jc_tfa_normPDF(x,dim)

if nargin < 2
    dim = 2;
end

x_norm = zeros(size(x));
if ~isvector(x)
    [m,n] = size(x);
    switch dim
        case 1
            k = m;
            for i = 1:k
                x_norm(i,:) = normPFD_vector(x(i,:));
            end
        case 2
            k = n;
            for i = 1:k
                x_norm(:,i) = normPFD_vector(x(:,i));
            end
    end
else
    x_norm = normPFD_vector(x);
end

end

%% sub-function
function x_norm = normPFD_vector(x)
    x(isnan(x)) = [];
    % level up x to zero
    x = x - min(x);
    % normalize x to sum up to 1
    x_norm = x./sum(x);
end
