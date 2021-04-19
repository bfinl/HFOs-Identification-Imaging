%% High-frequency Oscillations Imaging
% main function
% This script images identified and extracted HFOs in scalp M/EEG
% recordings to the cortical space.
%
% The implementation first reads the modeled individual head model (BEM)
% and the extracted HFOs/spikes signals. Then, the cortical source
% distribution is estimated for the input signals using sLORETA, 
% Standardized low resolution brain electromagnetic tomography. The
% estimated results are shown in the cortical space.
% 
% The input leadfiled from individual BEM is provided in leadField.mat,
% which was generated in Curry 7 (Compumedics, NC, USA). The input signals
% (scalp EEG measurements) are extracted and prepared by the
% main_hfoDetection.m. As an alternative, the user can simply skip the
% detection process and run this script by using the prepared HFOs/spikes
% data, provided in ".\data\process", which has been identified and
% extracted by running the main_hfoDetection.m.
%
%%% REFERENCE
% The main imaging method was implemented with the sLORETA algorithm.
% R. D. Pascual-Marqui, "Standardized low-resolution brain electromagnetic
% tomography (sLORETA): technical details." Methods Find Exp Clin Pharmacol
% 24 Suppl D: 5-12, 2002.
%
%%% CITATION
% Please cite the following paper in your publications or presentations if
% this project or part of the codes or the data, provided here, helps your
% research.
% Z. Cai, A. Sohrabpour, H. Jiang, S. Ye, B. Joseph, B. H. Brinkmann, 
% G. A. Worrell, and B. He, "Noninvasive High-frequency Oscillations 
% Riding Spikes Delineates Epileptogenic Sources." PNAS April 27, 2021 
% 118 (17) e2011130118; https://doi.org/10.1073/pnas.2011130118.
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
% 2020.08.21
% Document and commit for repository.


%% clean workspace
close all; clear; clc;

%% setup path and folder
% code and data path
defaultBaseDataPath = setupPathAndFolder();
% raw data folder
pathRaw = fullfile(defaultBaseDataPath,'raw');
% input folder
pathInput = fullfile(defaultBaseDataPath,'process');
% output folder
pathOutput = fullfile(defaultBaseDataPath,'imaging');
if ~isfolder(pathOutput); mkdir(pathOutput); end
cd(pathInput);

%% setup imaging data
setupImagingData;

%% ESi inverse solver - HFOs
% initialize esi, source imaging object
esi = jc_esi(sig_hfo,fileLeadFieldLocation);
% setup sensor space
esi = setSensorSpace(esi,sigSpot_hfo,noiseCov_hfo,pathChannel);
% setup solver
esi = setInverseOperator(esi,'sLORETA',paramMethod);
% solve for inverse 
esi = inverseSolver(esi,'pathOutput',pathOutput,'fSolveTimeCourse',false);
% plot and save solution
%%% to save figure and imaging results
% saveSolution(esi,pathOutput,pathChannel,true,true,'HFOs');
%%% for visualization only
saveSolution(esi,pathOutput,pathChannel,false,false,'HFOs');

% output log files
jc_print_block('ESI finished for HFOs.');

%% ESi inverse solver - spikes
% initialize esi, source imaging object
esi = jc_esi(sig_spk,fileLeadFieldLocation);
% setup sensor space
esi = setSensorSpace(esi,sigSpot_spk,noiseCov_spk,pathChannel);
% setup solver
esi = setInverseOperator(esi,'sLORETA',paramMethod);
% solve for inverse 
esi = inverseSolver(esi,'pathOutput',pathOutput,'fSolveTimeCourse',false);
% plot and save solution
%%% to save figure and imaging results
% saveSolution(esi,pathOutput,pathChannel,true,true,'Spike');
%%% for visualization only
saveSolution(esi,pathOutput,pathChannel,false,false,'Spike');

% output log files
jc_print_block('ESI finished for spikes.');

