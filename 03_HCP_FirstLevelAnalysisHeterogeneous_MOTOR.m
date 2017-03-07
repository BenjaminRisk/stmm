%------------------------------------------------------------------
% Benjamin Risk
% Estimate temporal variance parameters for STMM1.
% Includes cross-validation to smooth temporal parameters
% 
% Questions? bbr28@cornell.edu
% Please report bugs
%------------------------------------------------------------------

%% Begin user input ----------->
cd <CD TO stmm DIRECTORY>
addpath './mfunctions'

taskfolder = 'tfMRI_MOTOR';
task = 'MOTOR';
nCov = 40; %number of covariates: 3*(number tasks)+2*6+2*5 (6 motion parameters
           %for each session, 5 terms in spline basis for each session
nSession = 284; %number of timepoints in session
nVertex = 29716; % equals 29716 for right cortex and 29696 for left cortex
mySubjectDir = ['./Data/' task '/Subjects'];
myResultsDir = ['./ResultsAndFigures/' task];
unix(['mkdir -p ' myResultsDir]);


% Create a contrast matrix of dimension (n contrast) x (number of
% covariates including nuisance)
contrastMat = [-1/5 0 0 -1/5 0 0 1 0 0 -1/5 0 0 -1/5 0 0 -1/5];
contrastMat = [contrastMat zeros(1,nCov-length(contrastMat))];

load(['./Data/' task '/subjects.mat'])
load ./supportingdatafiles/latlong.mat
load ./supportingdatafiles/matrix_GordonRev2_CII_CORTEX_RIGHT.mat
%% <------------ end user input





%% Main estimation loop:
tic;
fitfirstlevel(subjects,mySubjectDir,myResultsDir);
toc
% 75 minutes

%% Conduct MUMM analysis using GLS estimators
load([myResultsDir '/UnsmoothedARparameters.mat']);
secondlevelMUMM(contrastMat,dhatsGLSAllSubjects,myResultsDir)

% Create gifti files: 
load([myResultsDir '/TMUMMAllSubjects'])

% NOTE: One will need to specify "CortexRight" when these GII files are
% loaded into wb_view
createFuncCIItoGIIgeneric(squeeze(subjResultsContrast(:,1,:))','CORTEX_RIGHT',[myResultsDir '/' task '_TMUMM_subjectContrasts']);
createFuncCIItoGIIgeneric(resultsBetasTMUMM,'CORTEX_RIGHT',[myResultsDir '/' task '_TMUMM_Betas'])
createFuncCIItoGIIgeneric(resultsTstatBetasTMUMM,'CORTEX_RIGHT',[myResultsDir '/' task '_TMUMM_TstatBetas'])
createFuncCIItoGIIgeneric(resultsContrastsTMUMM,'CORTEX_RIGHT',[myResultsDir '/' task '_TMUMM_Contrasts'])
createFuncCIItoGIIgeneric(resultsTstatContrastsTMUMM,'CORTEX_RIGHT',[myResultsDir '/' task '_TMUMM_TstatContrasts'])

createFuncCIItoGIIgeneric(squeeze(phiAllSubjects(:,1,:))','CORTEX_RIGHT',[myResultsDir '/' task '_AllSubjects_Unsmoothed_AR3_phi1'])
createFuncCIItoGIIgeneric(squeeze(phiAllSubjects(:,2,:))','CORTEX_RIGHT',[myResultsDir '/' task '_AllSubjects_Unsmoothed_AR3_phi2'])
createFuncCIItoGIIgeneric(squeeze(phiAllSubjects(:,3,:))','CORTEX_RIGHT',[myResultsDir '/' task '_AllSubjects_Unsmoothed_AR3_phi3'])
createFuncCIItoGIIgeneric(varInnovAllSubjects','CORTEX_RIGHT',[myResultsDir '/' task '_AllSubjects_Unsmoothed_InnovVariance'])


%% Determine optimal bw for smoothing the first-level variance parameters:


hvec = [1.15 1.25 1.5 2 2.25 2.5 3 4];
% NOTE: 
%   1.15: median 2 neighbors, many with zero;  1.25: median 4 neighbors;
%   1.5:  median 6 neighbors; 2: median 8 neighbors; 2.25: median 12;
%   2.5: 16; 3: 19; 4: 38;

tic;
[mseCortex,hvecCortex,h_ob] = gcv_ARparametersAllSubjects(hvec,latR,longR,varResAllSubjects,phiAllSubjects);
toc

clf;
a=figure;
plot(hvecCortex',mseCortex');
line([h_ob h_ob],[0 5],'Color','r')
title('5-fold cross-validation across right cortex for 30 subjects')
xlabel('bandwidth (mm)')
ylabel('MSE')
ylim([0,2])


%% Smooth AR parameters with optimal bandwidth and calculate K_n' Psi_n K_n:
tic;
smoothARparameters(h_ob,myResultsDir,mySubjectDir,subjects,contrastMat,latR,longR,'CORTEX_RIGHT');
toc

% Create gifti files:
load([myResultsDir '/SmoothedARparameters.mat'])

createFuncCIItoGIIgeneric(squeeze(smPhiAllSubjects(:,1,:))','CORTEX_RIGHT',[myResultsDir '/' task '_AllSubjects_Smoothed_AR3_phi1'])
createFuncCIItoGIIgeneric(squeeze(smPhiAllSubjects(:,2,:))','CORTEX_RIGHT',[myResultsDir '/' task '_AllSubjects_Smoothed_AR3_phi2'])
createFuncCIItoGIIgeneric(squeeze(smPhiAllSubjects(:,3,:))','CORTEX_RIGHT',[myResultsDir '/' task '_AllSubjects_Smoothed_AR3_phi3'])
createFuncCIItoGIIgeneric(smVarInnovAllSubjects','CORTEX_RIGHT',[myResultsDir '/' task '_AllSubjects_Smoothed_InnovVariance'])

     


