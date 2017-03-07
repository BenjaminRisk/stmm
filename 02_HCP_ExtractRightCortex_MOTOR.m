%-----------------------------------------------------------
% Benjamin Risk
% Description: This code extracts surface data from the CIFTI
% files. It concatenates the two phase encoding runs (with RL first and LR second).
% The code extracting the cortex from the CIFTI format requires
% wb_command. 

% NOTE: Depending on how you installed wb_command, you may need to modify
% <readincortex.m> so that it finds the command, see <readincortex.m>
% if it produces an error 
%
% Questions? bbr28@cornell.edu
%-----------------------------------------------------------


%% Begin user input------>
cd ~/Dropbox/SpatioTemporalFMRI_FixedVertex/FullAnalysis/stmm
addpath './mfunctions'

% Provide path in which the files in SpatiotemporalFMRI_FixedVertex have
% been saved:
myHCP_path = '/media/ben/UR100-20140805';

% Directory containing surface fMRI data:
foldername = 'tfMRI_MOTOR';

% Task:
task = 'MOTOR';

% Load list of subjects saved in 01_CreateSubjectDesignMatrices_MOTOR.m:
load(['./Data/' task '/subjectsMOTOR.mat'])

% Base directory in which to save cortical data in .mat format:
myHCPforMatFiles = ['./Data/' task '/Subjects'];
%% <-------- End user input



% The following code takes 8.5 minutes for thirty subjects on 3.66 GHz x 8:
myNSubject=length(subjects);
tic;
for j = 1:myNSubject
     subject = char(subjects(j));
     file1 = [myHCP_path, '/', subject, '/MNINonLinear/Results/',foldername,'_RL/' foldername '_RL_Atlas.dtseries.nii'];
     file2 = [myHCP_path, '/',subject, '/MNINonLinear/Results/',foldername,'_LR','/' foldername '_LR_Atlas.dtseries.nii'];
     fData = readincortex(file1,file2,'CORTEX_RIGHT');
     save([myHCPforMatFiles,'/',subject,'/matData_RIGHT_CORTEX.mat'],'fData');    
end
toc
