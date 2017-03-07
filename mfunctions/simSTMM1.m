function [fSimdata,subjSimData,sampleVar] = simSTMM1(allSubjDesignMat,betaFixed,Ssigma,Bsigma,Btheta,...
    phiArray,varInnovMatrix,distMat,dhatArray)
%simSTMM1 simulate the spatial temporal mixed model with fixed population vertex effects
% nuisance parameters are set equal to first level estimates from ToM analysis of HCP dataset
%INPUT: 
%   allSubjDesignMat: N x T x (Q+R) array of design matrices for all
%   subjects
%   betaFixed: Vx2 vector of the fixed effect coefficients for the overall effect of the
%       tasks at each vertex
%       -if size 1x2, then the same coefficients are used for all locations
%   Ssigma: 1x2 vector of variances of subject random effects
%   Bsigma: 1x2 vector of variances for the vertex-subject random effects
%   Btheta: 1x2 vector of dependence parameters for the exponential
%       variogram for the vertex random effects
%   phiArray: nSubject x 3 x nParcelm array 
%				- this can also be input as a 1 x 3 vector, in which case
%				  the same AR parameters are used for all vertices and all subjects
%   varInnovMatrix: nSubject x nParcelm matrix of innovation variance of
%                   ARp process 
%					- this can also be input as a scalar, in which case the
%					  same variance is used for all subjects and all vertices
%   distMat: nParcelm x nParcelm array
%   dhatArray: nSubject x (Q+R) x nParcelm -- determines the coefficients for the nuisance parameters
% OUTPUT:
%   fSimdata: N x V x T dataset, where N is the number of subjects, V is the number of vertices
%       (determined by mParcel), and T=548 is the number of timepoints.
%   subjSimData: N x V x 2: total vertex activation for each subject, sums
%   fixed and random effects
%   sampleVar: sample variance of true random effects; "oracle" estimator
%   of variance components
nParcelm = size(distMat,1);
nSubject = size(allSubjDesignMat,1);

if numel(varInnovMatrix) ~= nSubject*nParcelm
	varInnovMatrix = varInnovMatrix*ones(nSubject,nParcelm);
end


if numel(phiArray) ~= nSubject*3*nParcelm
	phiArray2 = kron(ones(nSubject*nParcelm,1),phiArray);
    phiArray2 = reshape(phiArray2,[nSubject,nParcelm,3]);
	phiArray = permute(phiArray2,[1,3,2]);
end

if numel(betaFixed)==2
    betaFixed = ones(nParcelm,1)*betaFixed;
end

%% set-up parameters of subject-vertex random effects:
funk = @(h,thetas)thetas(1)*exp(-thetas(2)*h)+thetas(3);

bThetaMental = [Bsigma(1),Btheta(1),0];
OmegaMental = funk(distMat,bThetaMental);

bThetaRandom = [Bsigma(2),Btheta(2),0];
OmegaRandom = funk(distMat,bThetaRandom);

%% Generate BOLD signal:
fSimdata = zeros(nSubject,nParcelm,548);
subjSimData = zeros(nSubject,nParcelm,2);
%Generate subject random effects:
    sMental = sqrt(Ssigma(1))*randn(nSubject,1);
    sRandom = sqrt(Ssigma(2))*randn(nSubject,1);
    sMat = [sMental,sRandom];

    %Generate subject-vertex interaction random effects:
    bMental = mvrandn(nSubject,OmegaMental);
    bRandom = mvrandn(nSubject,OmegaRandom);
          
    sampleVar.Ssigma =  var(sMat);
    sampleVar.Bsigma = [var(bMental(:)),var(bRandom(:))];
    for n=1:nSubject    
        subjDesignMat = squeeze(allSubjDesignMat(n,:,:));
        xMat = subjDesignMat(:,[1 4]);
        zMat = subjDesignMat(:,[2 3 5:end]);
        Btemp = [bMental(:,n),bRandom(:,n)];
        gammasn = squeeze(dhatArray(n,[2 3 5:end],:));
        tempVar = varInnovMatrix(n,:);
        % generate AR errors:
        parfor v=1:nParcelm
            phivector = phiArray(n,:,v);
            a_nv = simAR3_2sess(548,phivector,tempVar(v));
            subjdata = betaFixed(v,:)'+Btemp(v,:)'+sMat(n,:)';
            subjSimData(n,v,:) = subjdata;
            fSimdata(n,v,:) = xMat*subjdata+zMat*gammasn(:,v)+a_nv;
        end
    end
end

