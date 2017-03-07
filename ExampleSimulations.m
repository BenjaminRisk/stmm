%-------------------------------------
% Benjamin Risk
% Simulate data and estimate STMM.
%-------------------------------------

% load the subject design matrices and some of the parameters
% used in simulating data:
cd <CD TO stmm DIRECTORY>
addpath ./mfunctions
load ./supportingdatafiles/filesforsim.mat

[nSubject,nTime,nCov] = size(allSubjDesignMat);
myBetaFixed = [31,0]; 

% NOTE: myInvMAllSubjects contains matrices used to calculate Worsley's
% reduced bias AR parameters; see equation S.4 and S.5; it is calculated outside of
% the estimation loop to improve speed.

%parameters correspond to lo-hi-hi
[fData,sData] = simSTMM1(allSubjDesignMat,myBetaFixed,[423 423],[2346 2346],[0.23,0.23],...
        smPhiAllSubjectsParcel,smVarInnovAllSubjectsParcel,parcelDist,dhatParcel);
% fData: 30 x 215 x 548 are the simulated time series for each location and
% each subject
% sData: 30 x 215 x 2 are the true subject effects for "mental" and
% "random"

mkdir ./Simulation/Covariograms
tic;
[estSTMM,estMUMM] = fitSTMM1(fData,parcelDist,[1.5,2,2.5,3],allSubjDesignMat,[1 0 0 1 zeros(1,24)],...
    [1 -1],myInvMAllSubjects,'./Simulation','./Simulation/Covariograms/simulationCovario');
toc
% Takes 2.5 minutes on an i5 1.8x4 GHz
% Takes 23 seconds on i7 3.6x8 GHz

% Compare MSEs
temp = permute(sData,[1,3,2]);
mean((temp(:) - estSTMM.eBLUPsBetaHat(:)).^2)
mean((temp(:) - estMUMM.dhats(:)).^2)


% plot an example subject:
clf;
subplot(1,3,1)
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), sData(1,:,1)','EdgeColor','None');
axis equal
colorrange = caxis();
title('Truth for a subject')

subplot(1,3,2)
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), squeeze(estSTMM.eBLUPsBetaHat(1,1,:))','EdgeColor','None');
axis equal
title('STMM')
caxis(colorrange);

subplot(1,3,3)
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),squeeze(estMUMM.dhats(1,1,:))','EdgeColor','None');
axis equal
title('MUMM')
caxis(colorrange);

