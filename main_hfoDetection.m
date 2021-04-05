%% High-frequency Oscillations Detection
% main function
% This script detects HFOs in scalp M/EEG recordings based on a
% detection-discrimination framework.
% This framework applys to both simulation and clinical analysis, and the
% sample data can be selected accordingly.
%
%%% BRIEF ABOUT THE FRAMEWORK
% In the first stage, a high sensitive detector is used to create a big HFO
% pool. Then, the candidates are sieved to exclude those noisy activities.
% Next, distinguishing features are extracted from the time-frequency 
% ddomain of the raw and high-passed EEG data, with which a Gaussian
% mixture model operates to isolate putative HFO events from other spurious
% activities.
%
% The main idea here is to use data-specific features to help isolate HFO
% events and avoid confounding HF activities from noise, spikes, or
% filtering process.
%
% This script mainly operates in the sensor domain, and a separate script,
% main_hfo_Imaging, deals with source imaging, which projects the sensor
% level activities to the source domain for estimation of the brain sources.
%
% For simulation dataset, there are two sets of data included. The user can
% switch the dataset in loadData_simulation.m.
%
%%% REFERENCE
% The main identification framework was inspired by several studies.
% S. Liu et al., "Exploring the time-frequency content of high frequency
% oscillations for automated identification of seizure onset zone in
% epilepsy." J Neural Eng 13, 026026 (2016).
%
% J. A. Blanco et al., "Unsupervised classification of high-frequency
% oscillations in human neocortical epilepsy and control patients." J
% Neurophysiol 104, 2900-2912 (2010).
%
% V. V. Nikulin et al., "A novel method for reliable and fast extraction of
% neuronal EEG/MEG oscillations on the basis of spatio-spectral
% decomposition." Neuroimage 55, 1528-1535 (2011).
%
%%% CITATION
% Please cite the following paper in your publications or presentations if
% this project or part of the codes or the data, provided here, helps your
% research.
% Z. Cai et al., "Noninvasive High-frequency Oscillations Riding Spikes
% Delineates Epileptogenic Sources." Proc. Natl. Acad. Sci., 2021 (DOI).
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
% setup output subfolder
pathOutput = fullfile(defaultBaseDataPath,'process');

%% load data file
%%% clinical dataset
loadData_clinical; fClinical = true;
%%% simulation dataset
% loadData_simulation; fClinical = false;

%% initial detection
main_hfoDetection_detection;

%% feature extraction
main_hfoDetection_feature;

%% clustering
main_hfoDetection_clustering;

%% save hfo results
if fClinical
    %%% for clinical data analysis
    main_hfoDetection_saveHFO;
else
    %%% for simulation evaluation
    main_hfoDetection_simEval;
end

%%% process finished
jc_print_block('HFO identification finished.');








