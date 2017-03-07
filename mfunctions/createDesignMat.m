function [subjDesignMat,varNamesDesignMat] = createDesignMat(subjectSaveDir,taskDir1,taskDir2,...
    motionDir1,motionDir2,movementFileName,taskList,TR,nTR,nBin,nKnot,p)
%Create design matrix including tasks convolved with HRF, their
%derivatives, motion variables, and basis to capture drift. Assumes
%two runs as in the HCP data.
%Input:
%   subjectSaveDir
%   taskDir1
%   taskDir2
%   motionDir1
%   motionDir2
%   movementFileName
%   taskList
%   TR
%   nTR
%   nBin
%   nKnot
%   p

if nargin<12
    p = [6 16 1 1 6 0 32];
end

 % create tasks convolved with HRF along with their derivatives with
    % respect to onset and dispersion
    [designMatHRFrun1,varNames] = createTaskMatHRF(taskDir1,taskList,TR,nTR,nBin,p);
    designMatHRFrun2 = createTaskMatHRF(taskDir2,taskList,TR,nTR,nBin,p);
        
    %create motion covariates
    nuisanceRun1 = dlmread([char(motionDir1),'/',movementFileName,'.txt']);
    nuisanceRun2 = dlmread([char(motionDir2),'/',movementFileName,'.txt']);
    
    %HCP data contain additional columns that we are not interested in;
    %subset to the six parameters from affine registration.
    %We will allow the motion parameter coefficients to vary by run, i.e.,
    %an interaction with run.
    nuisMat = [[nuisanceRun1(:,1:6);zeros(nTR,6)],[zeros(nTR,6);nuisanceRun2(:,1:6)]];
    nuisVarNames = {'zTransX_run1','zTransY_run1','zTransZ_run1','zRotX_run1',...
        'zRotY_run1','zRotZ_run1','zTransX_run2','zTransY_run2','zTransZ_run2',...
            'zRotX_run2','zRotY_run2','zRotZ_run2'};
        
    % create piecewise linear spline to capture scanner drift. Create
    % separate splines for each run.
    driftBasis = createLinSpline(TR,nTR,nKnot);
    driftNames = {'zBasis1_sess1','zBasis2_sess1', 'zBasis3_sess1', 'zBasis4_sess1', ...
        'zBasis5_sess1','zBasis1_sess2','zBasis2_sess2', 'zBasis3_sess2', ...
        'zBasis4_sess2', 'zBasis5_sess2'};

    % full design matrix:
    subjDesignMat = [[designMatHRFrun1;designMatHRFrun2],nuisMat,...
        [driftBasis;zeros(nTR,nKnot)],[zeros(nTR,nKnot);driftBasis]];
    varNamesDesignMat = [varNames,nuisVarNames,driftNames];

    save([char(subjectSaveDir),'/subjDesignMat'],'subjDesignMat','varNamesDesignMat');
end

