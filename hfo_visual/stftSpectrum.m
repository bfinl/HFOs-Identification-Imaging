%% Short time Fourier transform
% This function is used for time-frequency spectrogram analysis on
% bio-signals. Several time windows are designed.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2017.02.11
% Initial implementation.
%
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function [tFT,lambdaFT,sClnWindFT] = stftSpectrum(sCln,fs,WindName,L,R,M,fPlot,fig,ax)

% Initial input parameters
% sample length
sCln = sCln(:);
N = length(sCln);
% sample time course
dt = 1/fs;
t  = 0:dt:dt*(N-1); t = t';
% window length
if ~exist('L','var') || isempty(L)
    L = floor(max(N/20,fs/10));
end
% figure options
if ~exist('fPlot','var') || isempty(fPlot)
    fPlot = true;
end
% output figure
if (~exist('fig','var') || isempty(fig)) && fPlot
    fig = figure;
    plot(1);
end
if (~exist('ax','var') || isempty(ax)) && fPlot
    ax = fig.CurrentAxes;
end
% DFT length
if ~exist('M','var') || isempty(M)
    M = min(512,20*L);
elseif isempty(M)
    M = min(1024, 2^nextpow2(20*L));
end
% window shift
if ~exist('R','var') || isempty(R)
    R = floor(max(L/50,10));
end

% Window type, WindName
% Rectangular, or Rct
if strcmp(WindName,'Rct') || strcmp(WindName,'Rectangular')
    Wind = rectwin(L)';
end
% Hamming, or Ham
if strcmp(WindName,'Ham') || strcmp(WindName,'Hamming')
    Wind = hamming(L,'periodic')';
end
% Gaussian, or Gauss, or Gabor
if strcmp(WindName,'Gauss') || strcmp(WindName,'Gaussian') || strcmp(WindName,'Gabor')
    Wind = gausswin(L)';
end

% STFT
% adjust length of signals
RAddN = mod(N,R);
if RAddN == 0; RAddN = L-R; else RAddN = L-RAddN; end % add more zeros to finish last few points
sCln = [sCln; zeros(RAddN,1)];

% initial data matrix
sClnTrnc = zeros(floor(N/R),L);
sClnWind = zeros(floor(N/R),L);
sClnWindZP = zeros(floor(N/R),M);
sClnWindFT = zeros(floor(N/R),M);
tFT = zeros(floor(N/R),1);
for i = 1:floor(N/R)
% truncate
sClnTrnc(i,:) = sCln((i-1)*R+1:(i-1)*R+L)';
sClnWind(i,:) = sClnTrnc(i,:).*Wind;
sClnWindZP(i,:) = [sClnWind(i,:) zeros(1,M-length(sClnWind(i,:)))];

sClnWindFT(i,:)  = fft(sClnWindZP(i,:),size(sClnWindZP,2));
tFT(i) = t((i-1)*R+1);
end

lambdaFT = linspace(0,1,size(sClnWindFT,2)/2)*(fs/2);
% sClnWindFT = 20*log(abs(sClnWindFT(:,1:end/2)))';
sClnWindFT = (abs(sClnWindFT(:,1:end/2)))';
if fPlot
    h = fig;
    pcolor(ax,tFT,lambdaFT,sClnWindFT), shading interp;
    set(gca,'ydir','normal');
    set(gca,'ylim',[min(lambdaFT) max(lambdaFT)]);
%     xlabel('Time (sec)');
%     ylabel('Frequency (Hz)','interpreter','tex');
    ylog, set(gca,'tickdir','in');
    colorbar;
    caxis([quantile(sClnWindFT(:),0.01),quantile(sClnWindFT(:),0.99)]);
end

end
