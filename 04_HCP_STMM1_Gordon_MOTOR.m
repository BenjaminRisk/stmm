%--------------------------------------------------
% Benjamin Risk
% Fit the STMM1
% Questions? bbr28@cornell.edu
%--------------------------------------------------

%% Begin user input ----------->
cd ~/Dropbox/SpatioTemporalFMRI_FixedVertex/FullAnalysis/stmm
addpath './mfunctions'

load ./supportingdatafiles/matrix_GordonRev2_CII_CORTEX_RIGHT.mat
load ./supportingdatafiles/latlong.mat

task='MOTOR';
side='RIGHT';
contrasts =  [-1/5 -1/5 1 -1/5 -1/5 -1/5]; %Dimension: (number of contrasts) x (number of main effects in contrast);
                     % Should sum to 1.
load(['./Data/' task '/subjects.mat'])
myResultsDir = ['./ResultsAndFigures/' task];
unix(['mkdir ' myResultsDir '/Covariograms']);

%% <------------ end user input

tic;
secondlevelSTMM(contrasts,myResultsDir,latR,longR,'CORTEX_RIGHT');
toc
% 40 minutes

load([myResultsDir '/STMMAllSubjects'])

% randomly chosen subject:
rng(453)
sindex = randsample(length(subjects),1);

sindex = 5;
subjects(sindex)

createFuncCIItoGIIgeneric(resultsBetas,'CORTEX_RIGHT',[myResultsDir,'/',...
    task,'_STMM1_Betas']);
createFuncCIItoGIIgeneric(resultsTstatBetas,'CORTEX_RIGHT',[myResultsDir,'/',...
    task,'_STMM1_TstatBetas'])
createFuncCIItoGIIgeneric(resultsContrasts,'CORTEX_RIGHT',[myResultsDir,'/',...
    task,'_STMM1_Contrasts'])
createFuncCIItoGIIgeneric(resultsTstatContrasts,'CORTEX_RIGHT',[myResultsDir,'/',...
    task,'_STMM1_TstatContrasts'])


% create file for subject random effects for lh task:
[nVertex,nParcelm] = size(matrix_GordonRev2_CII_CORTEX_RIGHT);
nSubject = length(subjects);

subjectRandomEffects = zeros(nVertex,nSubject);
for r = 2:nParcelm
    nVertexParcelm = sum(matrix_GordonRev2_CII_CORTEX_RIGHT(:,r));
    for n=1:nSubject
        subjectRandomEffects(matrix_GordonRev2_CII_CORTEX_RIGHT(:,r),n) = eBLUPsSubject(r,n,3);
    end
end



createFuncCIItoGIIgeneric(squeeze(subjResultsContrast)','CORTEX_RIGHT',[myResultsDir,'/',...
    task,'_STMM1_subjectContrasts'])

createFuncCIItoGIIgeneric(subjectRandomEffects,'CORTEX_RIGHT',[myResultsDir,'/MOTOR_STMM1_LH_eBLUPsSubject'])

createFuncCIItoGIIgeneric(eBLUPsSubjectVertex(:,:,3),'CORTEX_RIGHT',[myResultsDir,'/MOTOR_STMM1_LH_eBLUPsSubjectVertex'])





