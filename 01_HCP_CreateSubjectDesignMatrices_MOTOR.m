%---------------------------------------------------------------------
% Benjamin Risk
% Description: create design matrices.
% Here, design matrices are created for the HCP motor task. 
% Questions: bbr28@cornell.edu
%---------------------------------------------------------------------


%% Begin user input ---->

% In the example that follows, a design matrix is constructed for two runs
% from a paradigm, which corresponds to the HCP data.

% Provide path in which the files in SpatiotemporalFMRI_FixedVertex have
% been saved:
cd ~/Dropbox/SpatioTemporalFMRI_FixedVertex/FullAnalysis/stmm
addpath './mfunctions'

%% 1.
% The subject design matrices need to be saved to different folders. The code below creates a list
% of these folders for each subject.

% Specify the number of subjects:
myNSubject=30;

% Specify the path to the HCP .func.gii data:
myHCP_path = '/media/ben/UR100-20140805'; %ASSUMED TO BE ORGANIZED IN THE SAME DIRECTORY STRUCTURE AS AMAZON S3>

% Specify task folder directory:
task = 'MOTOR';

% Extract list of subjects:
load('~/Dropbox/SpatioTemporalFMRI_FixedVertex/FullAnalysis/stmm/Data/MOTOR/subjectsMOTOR.mat')


%The code below can be adapted to automate the creation of folders for each subject:
mySubjectSaveDir = cell(1,myNSubject);
for n=1:myNSubject
    mySubjectSaveDir(n) = {['./Data/' task '/Subjects/',char(subjects(n))]};
end

for n=1:myNSubject
    unix([ 'mkdir -p ' char(mySubjectSaveDir(n))]);
end

%% 2.
% Provide location of data for the first run; for the HCP data,
% use the function <make_HCPlist>
myParadigmRun1 = 'tfMRI_MOTOR_RL';
myTaskDirListInRun1 = make_HCPtaskDirList(myHCP_path,myParadigmRun1,subjects);

%% 3. 
% Provide location of data for the second run; for the HCP data,
% use the function <make_HCPlist>
myParadigmRun2 = 'tfMRI_MOTOR_LR';
myTaskDirListInRun2 = make_HCPtaskDirList(myHCP_path,myParadigmRun2,subjects);

%% 4.
% Provide a list of the paths to the movement regressors for run 1:
myMotionDirListInRun1 = make_HCPmotionDirList(myHCP_path,myParadigmRun1,subjects);

%% 5.
% Provide a list of the paths to the movement regressors for run2:
myMotionDirListInRun2 = make_HCPmotionDirList(myHCP_path,myParadigmRun2,subjects);

%% 6.
% Provide the file name of the movement regressors without the file extension.
% The name is assumed to be the same for all subjects and all runs. The
% file extension is assumed to be .txt.
myMovementFile = 'Movement_Regressors';

%% 7.
% 
% For the task onset and duration times, this code assumes there exists a text file 
% in the fsl "three column" format. The names of the files without file extension
% need to be provided. The file extension is assumed to be .txt.
% For HCP data, the task text files are found in the ./EVs folder. 
myTasks = {'cue','lf','lh','rf','rh','t'}; %for Motor
%myTasks = {'mental','rnd'}; %for ToM

%% 8.
%repetition time in seconds
myTR = 0.72; 

%% 9. 
% number of scans in each run;
myNTR = 284; %frames for Motor 
%myNTR = 274;  %frames for ToM

%% 10.
% number of timepoints used per TR in the discretization of the integral 
% used in convolving HRF with tasks. The actual weight (mesh size) used in
% the approximate integration is TR/nBins;      
myNBin = 16; 

%% 11. 
% rank of piece-wise linear spline basis to capture drift; note a separate basis 
% is used for each run
myNKnot = 5;    
%% <------------ end user input

%% 12.
% Parameters of the hrf. See documentation in <hrf.m> Commented out
% because the function defaults to the spm defaults.
% p = [6,16,1,1,0,6,32];

for n=1:myNSubject
    createDesignMat(mySubjectSaveDir(n),myTaskDirListInRun1(n),myTaskDirListInRun2(n),myMotionDirListInRun1(n),...
        myMotionDirListInRun2(n),myMovementFile,myTasks,myTR,myNTR,myNBin,myNKnot);
end

