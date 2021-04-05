To run the codes, Matlab (The MathWorks, Inc., MA, USA) has to be installed on the computer. The current codes have been developed with Matlab 2018a and tested also with Matlab 2017b, 2019b, and 2020b.

In addition, some of the codes require functions from third-party toolboxes for signal processing and visualization purposes, as listed in the following.

1. EEGLab is used for visualization of the topo plots, which can be downloaded and deployed from https://sccn.ucsd.edu/eeglab/index.php. The version tested is v13.6.5b and v14.1.1b (recommended) and probably compatible with versions above.
2. jLab is used for time-frequency visualization, which can be downloaded and deployed from https://github.com/jonathanlilly/jLab. The version of jLab tested is v1.6.2.
3. regtools is used for regularization of the source imaging process, which can be downloaded and deployed from http://www2.compute.dtu.dk/~pcha/Regutools/. The version of Regularization Tools tested is v4.1.
4. SSD, Spatio-Spectral Decomposition, is used for extraction and denoising of HFO activities, which can be downloaded from https://github.com/svendaehne/matlab_SSD/.
5. Optimal kmeans, kmeansElbow, is used for selecting optimal kmeans clusters for HFOs identification/extraction. This function was implemented by Sebastien De Landtsheer, https://it.mathworks.com/matlabcentral/fileexchange/65823-kmeans_opt.
6. findpeaks, built-in function in Matlab, with minor alterations for extra outputs, is used for locating the peaks in a signal.
7. tightSubplot, is used for setting axes with adjustable margins and gaps. This function was implemented by Pekka Kumpulainen and is available from the following link: https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w.

For the codes to correctly locate the third-party toolboxes and functions, the corresponding toolbox codes could be deployed into the "toolbox" folder and separated into the subfolders, namely, "eeglab", "jLab", "regtools", and “thirdParty”.
The functions kmeansElbow, findpeaks, and tightSubplot with minor modifications has been provided in the "thirdParty" folder already. These toolboxes will be loaded automatically when running the main workflow codes.

As an alternative, the user can also put the toolboxes or codes at other places, as long as the codes could be searched and called by Matlab. It should also be noted that if there already exist other versions of the same toolboxes in the Matlab search path, it is recommended to remove the duplicated ones before running the codes to avoid potential conflicts due to overwritten directories and compatibility.
