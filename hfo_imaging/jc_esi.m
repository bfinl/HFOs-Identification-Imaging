%% Electrophysiological source imaging class
% This file defines an ESI class with variables and functions for
% performing source imaging analysis on scalp-recorded EEG signals (spikes
% and HFOs).
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
% 2018.01.24
% Initial implementation.
%
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


classdef jc_esi
    
	properties (SetAccess = public)
        % sensor space
        sensorMeasurement
        sensorPhi
        sensorPhiMax
        nChannel
        nEpochLength
        nPhiLength
        signalRange
        channelList
        channelLoc
        channelBadList
        channelBadIndex = []
        pathChannel

        % sensor noise space
        noiseRange
        noiseCovariance

        % source space
        sourceEstimation
        sourceCost
        sourceLocation

        % forward operator
        leadField
        lfdLocation
        lfdOrientation
        nLfd
        meshTri
        meshLoc
        nOrn
        
        % inverse operator
        inverseMethod = 'sLORETA'
        paramMethod
        paramLambda
        nSource
        
        % identification
        patientID = 'P1'
        
        % parameter folders
        fileLeadField
        folderLeadField
        pathLeadField
        pathInput
        pathOutput
        
    end
    
    properties (Constant, Access = protected)
        % constants
        INVERSE_METHOD_LIST = {'sLORETA'};
        PARAM_METHOD_LIST = {'lcurve','dcp'};
    end
    % END of properties
   
    methods
        %% constructor
        % initialize ESI object
        function obj = jc_esi(sensorMeasurement,pathLeadField)
            
            if nargin < 2
                % lead-field data
                [fileLeadField,folderLeadField] = uigetfile('*.mat','Setup Lead-field',pathInputFolder);
                pathLeadField = [folderLeadField,fileLeadField];
            else
                if ~strfind(pathLeadField,'.mat') %#ok<STRIFCND>
                    pathLeadField = [pathLeadField,'.mat'];
                end
                [folderLeadField,fileLeadField] = jc_bs_findFileFolder(pathLeadField);
            end
            obj.pathLeadField   = pathLeadField;
            obj.folderLeadField = folderLeadField;
            obj.fileLeadField   = fileLeadField;
            
            % load leadfield file
            load(obj.pathLeadField);
            % setup leadfield information
            % number of locations
            nMeshLocation = length ( curryloc );
            % number of values and components per value
            nLfd = length ( currylfd ); %#ok<USENS>
            nOrn = ceil ( nLfd / nMeshLocation );
            nSource = nLfd / nOrn;
            % extract lfd matrix
            leadField = currylfd(7:end,:);
            lfdLocation = currylfd(1:3,1:nOrn:end);
            lfdOrientation = currylfd(4:6,1:nOrn:end);

            [nChannelSensor,nEpochLength] = size(sensorMeasurement);
            [nChannelLfd,~] = size(leadField);

            if nChannelLfd==nChannelSensor
               % sensor epoch has its rows equal to lfd rows in dimension
               obj.nChannel = nChannelSensor;
               obj.nEpochLength = nEpochLength;
               obj.nSource = nSource;
               obj.sensorMeasurement = sensorMeasurement;
               obj.leadField = leadField;
               if nChannelLfd==nEpochLength
                   % sensor epoch sizes are equal
                   % by default, rows are channels, columns are time points
                   warning('Sensor epoch has same number of rows and columns.');
               end
            else
               if nChannelLfd~=nEpochLength
                   % sensor epoch has no dimension that equals to lfd rows
                   error('Sensor epoch does not match lead field.');
               else
                   % sensor epoch might be flipped
                   sensorMeasurement = sensorMeasurement.';
                   obj = esi(sensorMeasurement,fileLeadField);
                   return
               end
            end
            
            % setup object variables
            obj.nLfd = nLfd; obj.nOrn = nOrn;
            obj.meshTri = currytri; obj.meshLoc = curryloc;
            obj.lfdLocation = lfdLocation;
            obj.lfdOrientation = lfdOrientation;
            end
        
        %% setup sensor space
        % setup sensor measurement range for signal, noise, and noise covariance
        % note here, noise covariance is in the unit of energy or power
        function obj = setSensorSpace(obj,signalRange,noiseCovariance,pathChan)
            if ~isnumeric(obj.sensorMeasurement)
               error('ESI object needs initialization first.');
            end
            % start setup sensor space signals
            if isvector(signalRange) && numel(signalRange)>1 && numel(signalRange)<4
               obj.signalRange = signalRange;
               if length(signalRange)==2
                   obj.signalRange(3) = floor(abs(signalRange(2)-signalRange(1))/2)+1;
               end
            else
               error('Signal range should be a vector with 2 or 3 elements.');
            end
            [obj.channelList,obj.channelLoc] = jc_io_readChannelList(pathChan);
            obj.noiseCovariance = noiseCovariance;
            obj.sensorPhi       = obj.sensorMeasurement(:,obj.signalRange(1):obj.signalRange(2));
            obj.sensorPhiMax    = obj.sensorPhi(:,obj.signalRange(3));
            obj.nPhiLength      = size(obj.sensorPhi,2);
        end
        
        %% setup inverse model
        % setup inverse method and parameter estimation method
        function obj = setInverseOperator(obj,inverseMethod,paramMethod)
           if nargin == 1
               obj.inverseMethod = 'sLORETA';
               obj.paramMethod = 'lcurve';
           else
               switch inverseMethod
                   case obj.INVERSE_METHOD_LIST
                       obj.inverseMethod = inverseMethod;
                   otherwise
                       obj.inverseMethod = 'sLORETA';
                       warning('Invalid inverse method. Use sLORETA as default.');
               end
               if nargin == 3
                   switch paramMethod
                       case obj.PARAM_METHOD_LIST
                           obj.paramMethod = paramMethod;
                       otherwise
                           obj.paramMethod = 'lcurve';
                           warning('Invalid parameter estimation method. Use lcurve as default.');
                   end
               end
           end
        end
        
        %% inverse solver
        function obj = inverseSolver(obj,varargin)
            % initialize sovler input
            chkSingleNumFcn = @(x) validateattributes(x,{'numeric'},{'integer','size',[1 1]});
            p = inputParser();
			p.addParameter('phi',obj.sensorPhiMax,@ismatrix);
            p.addParameter('fSolveTimeCourse',true,@islogical);
			p.addParameter('Lfd',obj.leadField,@ismatrix);
			p.addParameter('nSource',obj.nSource,chkSingleNumFcn);
			p.addParameter('lambda',[],@iscolumn);
			p.addParameter('C_noise',obj.noiseCovariance,@ismatrix);
			p.addParameter('paramMethod',obj.paramMethod,@ismatrix);
			p.addParameter('sigSpot',obj.signalRange(3),chkSingleNumFcn);
            p.addParameter('pathOutput','../output/',@ischar);
            p.addParameter('meshTri',obj.meshTri,@ismatrix);
            p.addParameter('meshLoc',obj.meshLoc,@ismatrix);
			p.parse(varargin{:});
			s = p.Results;
            
            if s.fSolveTimeCourse
                T     = obj.nPhiLength;
                Phi   = obj.sensorPhi;
                s.phi = std(Phi,[],2); % used to solve for lambda
            else
                T = 1;
                Phi = obj.sensorPhiMax;
            end
            switch obj.inverseMethod
                case {'sLORETA'}
                    % solve for optimal lambda
                    s.lambda = [];
                    if isempty(s.lambda)
                        [~,~,s.lambda,Sj] = jc_sl_sLORETA(s);
                    else
                        Sj = [];
                    end
                    % allocate for solution space
                    J_set = zeros(size(s.Lfd,2),T);
                    cost_set = zeros(T,2);
                    % solver
                    for it = 1:T
                        clc;
                        fprintf(['Solving sLORETA... ',num2str(floor(it/T*100)),'%% done...\n']);
                        [J_set(:,it),cost_set(it,:)] = jc_sl_sLORETA(s,s.lambda,Phi(:,it),Sj);
                    end
            end
            
            % log outputs
            obj.paramLambda      = s.lambda;     % ESI hyperparameter
            obj.sourceEstimation = J_set;        % source imaging estimation
            obj.sourceCost       = cost_set;     % solution cost
            obj.pathOutput       = s.pathOutput; % output path
        end
        
        %% save solutions
        function obj = saveSolution(obj,pathOutput,pathChannel,...
                fSaveFig,fSaveLog,figTitle)
            if nargin < 2 || isempty(pathOutput)
                pathOutput = obj.pathOutput;
            end
            if nargin < 3
                pathChannel = [];
            end
            if nargin < 4 || isempty(fSaveFig)
                fSaveFig = true;
            end
            if nargin < 5 || isempty(fSaveLog)
                fSaveLog = true;
            end
            if nargin < 6
                figTitle = [];
            end
            obj.pathOutput = pathOutput;
            obj.pathChannel = pathChannel;

            if ~isfolder(pathOutput); mkdir(pathOutput); end

            % plot results
            p = struct('j',obj.sourceEstimation,...
                 'pathOutputFile',[obj.pathOutput,filesep,obj.inverseMethod],...
                 'fSaveFig',fSaveFig,...
                 'pathChannel',pathChannel,...
                 'badChannel',obj.channelBadIndex,...
                 'phi',obj.sensorPhiMax./obj.noiseCovariance,...
                 'mapThreshold',0.5,...
                 'fig',figure(),...
                 'figTitle',figTitle,...
                 'fPlotFig',true,...
                 'fPlotTopo',true,...
                 'XYZLim',[1,1;1,1;1,1],...
                 'meshTri',obj.meshTri,...
                 'meshLoc',obj.meshLoc,...
                 'nLfd',obj.nLfd);

            plotSource(obj,p);

            if fSaveLog
                % log file, for test
                evalin('base','save([pathOutput,filesep,esi.inverseMethod,''_log.mat'']);');
                jc_print_block([obj.inverseMethod,' solved / saved.']);
            end
        end
        
        %% sub-function for plotting
        function fig = plotSource(obj,varargin) %#ok<INUSL>
            
            % initialize parameters
            chkSingleNumFcn = @(x) validateattributes(x,{'numeric'},{'real','size',[1 1]});
            p = inputParser();
            p.addParameter('j',[],@(x) isvector(x)||isempty(x));
            p.addParameter('pathOutputFile',[],@ischar);
            p.addParameter('fSaveFig',true,@islogical);
            p.addParameter('pathChannel',[]);
            p.addParameter('badChannel',[]);
            p.addParameter('phi',[],@isvector);
            p.addParameter('fig',[],@ishandle);
            p.addParameter('figTitle',[],@ischar);
            p.addParameter('fPlotFig',true,@islogical);
            p.addParameter('fPlotTopo',true,@islogical);
            p.addParameter('XYZLim',ones(3,2),@ismatrix);
            p.addParameter('fViewPositionFlip',false,@islogical);
            p.addParameter('viewPosition',[],(@(x)(isvector(x)||isempty(x))));
            p.addParameter('alphaFactorCortex',1,chkSingleNumFcn);
            p.addParameter('DiffuseStrength',0.6,chkSingleNumFcn);
            p.addParameter('SpecularStrength',0.2,chkSingleNumFcn);
            p.addParameter('SpecularColorReflectance',0.2,chkSingleNumFcn);
            p.addParameter('mapAmplitude',1,chkSingleNumFcn);
            p.addParameter('mapThreshold',0,chkSingleNumFcn);
            p.addParameter('meshTri',[],@ismatrix);
            p.addParameter('meshLoc',[],@ismatrix);
            p.addParameter('nLfd',[],chkSingleNumFcn);
            p.parse(varargin{:});
            s = p.Results;

            % prepare parameters for plot
            % mesh size
            nLoc = length ( s.meshLoc );
            nTri = length ( s.meshTri );
            % color space
            cmap = cmap_data();

            % define plot section in lim XYZ
            % it is optional to plot only a section of source space
            % say half of the brain, otherwise the whole space, as default
            if ~isempty(s.XYZLim)
                XYZLim_tmp = s.XYZLim.*[min(s.meshLoc,[],2) max(s.meshLoc,[],2)];
                xlim = XYZLim_tmp(1,:);
                ylim = XYZLim_tmp(2,:);
                zlim = XYZLim_tmp(3,:);

                xidx = find(s.meshLoc(1,:)>xlim(1) & s.meshLoc(1,:)<xlim(2));
                yidx = find(s.meshLoc(2,:)>ylim(1) & s.meshLoc(2,:)<ylim(2));
                zidx = find(s.meshLoc(3,:)>zlim(1) & s.meshLoc(3,:)<zlim(2));

                xyzidx = intersect(xidx,yidx);
                xyzidx = intersect(xyzidx,zidx);

                vol_display_tri_idx = cell(length(xyzidx),1);
                for i = 1:length(xyzidx)
                    [~,vol_display_tri_idx{i}] = find(s.meshTri == xyzidx(i));
                end
                vol_display_tri_idx = vertcat(vol_display_tri_idx{1:end});
                vol_display_tri_idx = unique(vol_display_tri_idx);
            else
                vol_display_tri_idx = 1:length(s.meshTri);
            end

            %==========================================================================
            % start plotting
            fig = s.fig;
            if isempty(fig); fig = figure; end
            figure(fig);
            hold all;
            axis tight manual
            ax = gca;
            ax.NextPlot = 'replaceChildren';

            % prepare source space for plot
            % convert amplitude of source space
            % if no estimated sources, then take all-zero source space
            j_amp = s.j;
            if isempty(j_amp)
                j_amp = zeros ( nLoc, 1 );
            else
                if length(j_amp)/nLoc>2.5 % test if it's a 3D case
                    j_amp = bs_dpAmp(j_amp);
                end
                if s.mapAmplitude
                    j_amp = abs(j_amp)/max(abs(j_amp));
                    j_amp(isnan(j_amp)) = 0;
                else
                    j_amp = (j_amp)/max(abs(j_amp));
                    j_amp(isnan(j_amp)) = 0;
                end
            end
            % thresholding source amplitude for plot
            if s.mapThreshold > 0
                j_amp(j_amp < (max(j_amp)*s.mapThreshold)) = 0;
            end

            % prepare axis space for cortex model
            nMin = min(s.meshLoc,[],2)-10;
            nMax = max(s.meshLoc,[],2)+10;
            if sum(sum(s.XYZLim~=1))~=6
                XYZLim_idx = find(s.XYZLim~=1);
                XYZLim_orig = [nMin nMax];
                XYZLim_orig(XYZLim_idx) = XYZLim_tmp(XYZLim_idx);
                nMin = XYZLim_orig(:,1); nMax = XYZLim_orig(:,2); 
            end
            axis ( [nMin(1),nMax(1),nMin(2),nMax(2),nMin(3),nMax(3)] );
            axis equal;
            axis vis3d;

            % prepare source vector to be plotted (zero-padded, transposed)
            % then plot as 3D triangular mesh
            nJLen = length(j_amp);
            if nJLen == nTri
                VSource = j_amp;
                hpatchSource = trisurf(s.meshTri(:,vol_display_tri_idx)',...
                    s.meshLoc(1,:),s.meshLoc(2,:),s.meshLoc(3,:),VSource(vol_display_tri_idx));
                set(hpatchSource,'EdgeColor','None', 'FaceAlpha',s.alphaFactorCortex,...
                    'FaceLighting','gouraud','backfacelighting','reverselit',...
                    'DiffuseStrength',s.DiffuseStrength,...
                    'SpecularStrength',s.SpecularStrength,...
                    'SpecularColorReflectance',s.SpecularColorReflectance);
            else
                VSource = zeros ( nLoc, 1 );
                VSource(1:length(j_amp),1) = j_amp;
                hpatchSource = patch ( 'vertices',s.meshLoc','faces',...
                    s.meshTri(:,vol_display_tri_idx)','FaceVertexCData',VSource );
                set ( hpatchSource,'EdgeColor','none','FaceColor','interp',...
                    'FaceLighting','gouraud','facealpha',s.alphaFactorCortex,...
                    'backfacelighting','reverselit',...
                    'DiffuseStrength',s.DiffuseStrength,...
                    'SpecularStrength',s.SpecularStrength,...
                    'SpecularColorReflectance',s.SpecularColorReflectance);
            end
            axis off;
            colormap(cmap(end/2:end,:));
            caxis([0,0.99]);
            colorbar;
            colorbar('ticks',[-1 1],'ticklabels',{'Min','Max'});

            % set view position and lighting
            j_amp_main_idx = find(abs(j_amp)>0.9*max(j_amp));
            srcLoc_main = s.meshLoc(:,j_amp_main_idx);
            viewPosition = mean(srcLoc_main,2);
            viewAngle = [sign(viewPosition(1))*120,5];
            if s.fViewPositionFlip
                viewPosition(2) = -viewPosition(2);
            end
            light('Position',viewPosition);
            light('Position',[-viewPosition(1),-viewPosition(2),-viewPosition(3)]);
            if isempty(s.viewPosition)
                view(viewAngle);
            else
                view(s.viewPosition);
            end
            title([s.figTitle,' imaging']);
            
            % END of plotting
            %==========================================================================

            % save results
            % save figure
            if ~isempty(s.pathOutputFile) && s.fSaveFig
                savefig(fig,s.pathOutputFile);
                print(fig,[s.pathOutputFile,'_sl'],'-dpng','-r300');
            end
            % scalp map
            if ~isempty(s.phi) && s.fPlotTopo
                figTopo = plotTopoMap(s.phi,[],s.pathChannel,s.badChannel);
                title([s.figTitle,' topo-map']);
                if s.fSaveFig
                    print(figTopo,[s.pathOutputFile,'_topo'],'-dpng','-r300');
                end
            end

            %% sub-function: compute amplitude of sources
            function j_amp = bs_dpAmp(j_opt,nOrn)
                if nargin < 2
                    nOrn = 3;
                end

                if nOrn == 3
                    [nSource,nTime,nEpoch] = size(j_opt); %#ok<*PROPLC>
                    j_opt = reshape(j_opt,nSource,[]);
                    nSource = nSource/3;

                    j_amp = zeros(nSource,nTime*nEpoch);
                    for iTimeEpoch = 1:nTime*nEpoch
                        j_opt_tmp = reshape(j_opt(:,iTimeEpoch),3,[]);
                        j_amp(:,iTimeEpoch) = sqrt(sum(j_opt_tmp.^2,1));
                    end
                    j_amp = reshape(j_amp,nSource,nTime,nEpoch);
                else
                    j_amp = j_opt;
                end
            end % END of source amplitude subfunction
            
            %% sub-function: set color map
            function cmap = cmap_data()
                cmap = [0.0664062500000000,0.132812500000000,0.199218750000000;0.0755657327586207,0.146551724137931,0.206088362068966;0.0938846982758621,0.174030172413793,0.219827586206897;0.112203663793103,0.201508620689655,0.233566810344828;0.130522629310345,0.228987068965517,0.247306034482759;0.148841594827586,0.256465517241379,0.261045258620690;0.167160560344828,0.283943965517241,0.274784482758621;0.185479525862069,0.311422413793103,0.288523706896552;0.203798491379310,0.338900862068966,0.302262931034483;
                    0.222117456896552,0.366379310344828,0.316002155172414;0.231276939655172,0.380118534482759,0.322871767241379;0.249595905172414,0.407596982758621,0.336610991379310;0.267914870689655,0.435075431034483,0.350350215517241;0.286233836206897,0.462553879310345,0.364089439655172;0.304552801724138,0.490032327586207,0.377828663793103;0.322871767241379,0.517510775862069,0.391567887931035;0.332031250000000,0.531250000000000,0.398437500000000;0.375224497126437,0.560704022988506,0.437051005747126;
                    0.418417744252874,0.590158045977012,0.475664511494253;0.440014367816092,0.604885057471264,0.494971264367816;0.483207614942529,0.634339080459770,0.533584770114943;0.526400862068966,0.663793103448276,0.572198275862069;0.569594109195402,0.693247126436782,0.610811781609195;0.612787356321839,0.722701149425287,0.649425287356322;0.655980603448276,0.752155172413793,0.688038793103448;0.699173850574713,0.781609195402299,0.726652298850575;0.742367097701149,0.811063218390805,0.765265804597701;
                    0.785560344827586,0.840517241379310,0.803879310344828;0.807156968390805,0.855244252873563,0.823186063218391;0.850350215517241,0.884698275862069,0.861799568965517;0.893543462643678,0.914152298850575,0.900413074712644;0.936736709770115,0.943606321839080,0.939026580459770;0.958333333333333,0.958333333333333,0.958333333333333;0.951778017241379,0.933459051724138,0.901400862068966;0.945222701149425,0.908584770114943,0.844468390804598;0.938667385057471,0.883710488505747,0.787535919540230;
                    0.932112068965517,0.858836206896552,0.730603448275862;0.928834410919540,0.846399066091954,0.702137212643678;0.922279094827586,0.821524784482759,0.645204741379310;0.915723778735632,0.796650502873563,0.588272270114943;0.909168462643678,0.771776221264368,0.531339798850575;0.902613146551724,0.746901939655172,0.474407327586207;0.896057830459770,0.722027658045977,0.417474856321839;0.889502514367816,0.697153376436782,0.360542385057471;0.882947198275862,0.672279094827586,0.303609913793103;
                    0.876391882183908,0.647404813218391,0.246677442528736;0.873114224137931,0.634967672413793,0.218211206896552;0.866558908045977,0.610093390804598,0.161278735632184;0.863281250000000,0.597656250000000,0.132812500000000;0.849542025862069,0.561018318965517,0.132812500000000;0.835802801724138,0.524380387931035,0.132812500000000;0.822063577586207,0.487742456896552,0.132812500000000;0.808324353448276,0.451104525862069,0.132812500000000;0.794585129310345,0.414466594827586,0.132812500000000;
                    0.780845905172414,0.377828663793103,0.132812500000000;0.773976293103448,0.359509698275862,0.132812500000000;0.760237068965517,0.322871767241379,0.132812500000000;0.746497844827586,0.286233836206897,0.132812500000000;0.732758620689655,0.249595905172414,0.132812500000000;0.719019396551724,0.212957974137931,0.132812500000000;0.705280172413793,0.176320043103448,0.132812500000000;0.691540948275862,0.139682112068966,0.132812500000000;0.677801724137931,0.103044181034483,0.132812500000000;
                    0.664062500000000,0.0664062500000000,0.132812500000000];
            end % END of color map subfunction
            
        end % END of plot subfunction
       
	end
   % END of methods

end
% END of class



