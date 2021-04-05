%% Visualization
% This function visualizes scatter graph in 3D space with group labels.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [fig,colorGrp] = gscatter3(feature1,feature2,feature3,group,featureLabel)

if nargin < 4
    group = ones(length(feature1),1);
end
if nargin < 5
    featureLabel = [];
elseif isnumeric(featureLabel)
    featureLabel = num2str(featureLabel);
end

% setup figure
fig = gcf;
hold all;

% get index for which clusters are filled, not empty
iGrp_set = unique(group);
nGrp = length(iGrp_set);
lgd = cell(nGrp,1);
% setup colors
colorGrp = lines(nGrp);

% loop over all filled clusters
for iGrp = 1:nGrp
    iGrp_idx = find(group==iGrp_set(iGrp));
    fadeFactor = 0.8;
    faceCol = fadeFactor*colorGrp(iGrp,:) + (1-fadeFactor)*[1 1 1];
    edgeCol = colorGrp(iGrp,:);
    mrkSize = 48;
    fig(iGrp) = scatter3(feature1(iGrp_idx),feature2(iGrp_idx),feature3(iGrp_idx),...
        mrkSize,'filled','o',...
        'MarkerFaceColor',faceCol,...
        'MarkerEdgeColor',edgeCol,...
        'LineWidth',1);
    if ~isempty(featureLabel)
        text(feature1(iGrp_idx),...
            feature2(iGrp_idx),...
            feature3(iGrp_idx),...
            featureLabel(iGrp_idx));
    end
    lgd{iGrp} = ['Cls ',num2str(iGrp_set(iGrp))];
end
view(3), box on, grid on;
legend(lgd);

end

