% Curry 7 M-file, Platform: PCWIN, Created on: 11/7/2017 12:47:28 PM
%
% curryloc contains mesh locations
% currytri contains triangle vertex indices
% currylfd(1:3,:) contains leadfield locations
% currylfd(4:6,:) contains leadfield orientations
% currylfd(7:end,:) contains leadfield
%
% load mat file
load ('leadField.mat');
% number of locations
nLoc = length ( curryloc );
% number of values and components per value
nTot = length ( currylfd );
nCom = ceil ( nTot / nLoc );
nVal = nTot / nCom;
% channel (>= 7) to be plotted
nRow = 7;
% basis vector (1..nCom) to be plotted
nBas = 1;
% prepare vector of values to be plotted (zero-padded,transposed)
V = zeros ( nLoc, 1 );
V(1:nVal,1) = currylfd(nRow,nBas:nCom:end)';
% plot using patch command
clf ( 'reset' );
nMin = min(curryloc,[],2)-10;
nMax = max(curryloc,[],2)+10;
axis ( [nMin(1),nMax(1),nMin(2),nMax(2),nMin(3),nMax(3)] );
axis equal;
axis vis3d;
hpatch = patch ( 'vertices',curryloc','faces',currytri','FaceVertexCData',V );
set ( hpatch,'EdgeColor','none','FaceColor','interp','FaceLighting','phong','DiffuseStrength',0.8 );
camlight right;
lighting phong;
colorbar;
