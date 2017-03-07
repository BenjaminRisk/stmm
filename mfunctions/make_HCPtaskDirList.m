function pathsList = make_HCPtaskDirList(HCP_path,paradigm_phaseencoding,subjects)
%This function creates a list of the directories containing the text files
%specifying the onset and duration times (in sec) for a user-specified paradigm 
%for some number of subjects. The subjects are chosen in the order in which they appear 
%when using the <ls> unix command.
%   This function aids in building design matrices for the special case of
%   the HCP data.
% INPUT:
%   HCP_path:  string containing the path to the HCP sampler data
%
%   paradigm_phaseencoding:  string containing the name of the folder
%   containing the results of the desired paradigm along with 
%   the phase-encoding direction (either RL or LR), e.g., "tfMRI_MOTOR_RL"
%
%   subjects: cell containing subject IDs
%
% OUTPUT:
%   pathsList: a list of length nsubjects of the directories containing the
%   HCP data for each subject. (The subjects chosen are the first nsubjects
%   when listed using ls, which can be viewed by printing pathsList.
%

%create list of all subjects:
nSubject = length(subjects);
pathsList = cell(1,nSubject);
for n=1:nSubject
    subject = char(subjects(n));
    pathsList(n) = {[HCP_path,'/',subject,'/MNINonLinear/Results/',paradigm_phaseencoding,...
        '/EVs']};
end
end
