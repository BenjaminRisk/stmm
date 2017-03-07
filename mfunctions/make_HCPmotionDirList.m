function pathsList = make_HCPmotionDirList(HCP_path,paradigm_phaseencoding,subjects)
%This function creates a list of the directories containing the files with
%the 6 parameters from the affine registration to be used as movement
%regressors. The subjects are chosen in the order in which they appear when using the
%<ls> unix command.
%   This function aids in the creation of lists that can then be input into
%   functions that build matrices of motion parameters for the special case of
%   the HCP data.
% INPUT:
%   HCP_path:  string containing the path to the HCP sampler data
%
%   paradigm_phaseencoding:  string containing the name of the folder
%   containing the results of the desired paradigm along with 
%   the phase-encoding direction (either RL or LR), e.g., "tfMRI_MOTOR_RL"
%
%   subjects: cell containing the subject IDs (strings)
%
% OUTPUT:
%   pathsList: a list of length nsubjects of the directories containing the
%   HCP motion data for each subject. (The subjects chosen are the first nsubjects
%   when listed using ls, which can be viewed by printing pathsList.)
%

nSubject = length(subjects);
pathsList = cell(1,nSubject);
for n=1:nSubject
    subject = char(subjects(n));
    pathsList(n) = {[HCP_path,'/',subject,'/MNINonLinear/Results/',paradigm_phaseencoding]};
end
end
