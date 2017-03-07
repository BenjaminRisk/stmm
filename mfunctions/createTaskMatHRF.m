function [designMat,varNames] = createTaskMatHRF(pathTasks,tasks,TR,nTR,nBin,p)
%createTaskMatHRF This function convolves the HRF specified by the
%parameters in p with the task onsets and durations. 
%INPUT:
% pathTasks: a string pointing to the directory containing
% text files that are in the "three column" format for fsl with the first time denoting
% the start times (with the first scan at time zero, in seconds) and the 
% second column denotes the task duration (sec). Note the third column denotes
% the value but is ignored -- we assume an indicator variable such
% that the value is equal to one during the task.
%
% tasks: a list of the task (i.e., cue) names corresponding to the names
% of the text files that appear in the directory pathTasks. For HCP, these
% files appear in the "./EVs" folder.
%
% TR: time between scans
%   
% nTR: number of scans in a run
%
% nBin - This parameter specifies the number of points in each TR that are
% used in the discrete approximation to the integration in the convolution of the HRF 
% with stimulus onsets and durations. TR/nBin = mesh size; a larger nBin leads to 
% a more accurate approximation, although the computation time increases
% somewhat. Defaults to 16. 
%
% p    - parameters of the response function (two Gamma functions).
% Defaults are taken from spm.
%The following is copied exactly from <spm_hrf.m> --->
%        p(1) - delay of response (relative to onset)          6
%        p(2) - delay of undershoot (relative to onset)       16
%        p(3) - dispersion of response                         1
%        p(4) - dispersion of undershoot                       1
%        p(5) - ratio of response to undershoot                6
%        p(6) - onset (seconds)                                0
%        p(7) - length of kernel (seconds)                    32
% <---- end copy
% Brisk NOTE: p(7) is NOT a parameter of the HRF but rather
% the time at which the hrf function is no longer calculated.
% The support of the function is [0,Inf), so this actually corresponds to
% the points at which we truncate the HRF.
%
%OUTPUT:
% designMat: number of scans (nTR) by 3Q (includes onset and dispersion
% derivatives for each task)
%
% varNames: column headings for designMat

Q = length(tasks);
runTime = TR*(nTR-1); % duration (sec) of the run; the first scan occurs at time 0.
meshSize = TR/nBin;
time = 0:meshSize:runTime;      
nTime = length(time); %equal to (nTR-1)*nBin+1

%default parameters of the HRF:
if nargin<6
    p = [6 16 1 1 6 0 32];
end

%create HRF function including onset and dispersion derivatives
hrf_mat= hrf(TR,nBin,p);

% convolve HRF and derivatives with each task:
designMatLong = zeros(nTime,Q*3);
for q=1:Q
        indicatorLong = zeros(nTime,1);
        taskData = dlmread([char(pathTasks),'/',char(tasks(q)),'.txt']);
        startTask = taskData(:,1);
        endTask = startTask+taskData(:,2);        
        %create time points corresponding to the resolution of the discretized integral:         
        %flag timepoints that fall within a task:
        for t=1:nTime
            indicatorLong(t)=logical(sum(time(t)>=startTask & time(t)<endTask));
        end
        for j=1:3
            temp = meshSize*conv(indicatorLong,hrf_mat(:,j));
            designMatLong(:,(3*(q-1)+j)) = temp(1:nTime);
        end
end
    designMat = designMatLong(1:nBin:nTime,:);
    varNames = cell(1,3*Q);
    for q=1:Q
        index = 3*(q-1)+1;
        varNames(index) = tasks(q);
        varNames(index+1) = {['d' char(tasks(q)) '_dt']};
        varNames(index+2) = {['d' char(tasks(q)) '_dd']};
    end
end


