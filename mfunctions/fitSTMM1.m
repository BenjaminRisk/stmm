function [outSTMM,outMUMM] = fitSTMM1(fdata,distmat,bandwidth,allsubjdesignmat,xindices,contrastmat,invMAllSubjects,resultsDir,savecovario)
%fitSTMM1 
% 	Benjamin Risk
%	November 10, 2015
%   Fit first and second level models of STMM1. Also produces estimates
% 	 from MUMM. The input is typically data for a parcel rather than the 
%    whole cortex, such that V = number of vertices in parcel   
%INPUT: 
% fdata - N x V x T data matrix
% distmat - V x V dense distance matrix
% bandwidth - either scalar or vector. If scalar, no cross-validation is
%             performed and the data are smoothed at bw. If a vector,
%             then gcv is used to select the optimal value from those in
%             the vector bandwidth.
% xindices - indices in the design matrices corresponding to
%   covariates-of-interest
% contrastmat nContrasts x nCOI where nCOI is number of covariates of
%   interest. Ordering corresponds to the non-zero entries of xindices
% invMAllSubjects
% resultsDir
% savecovario
nCov=size(allsubjdesignmat,3);
nQ = sum(xindices~=0);
nSession=2;
[nSubject,nVertex,nTR] = size(fdata);
nTR = nTR/nSession;
if numel(xindices)==nCov 
    xindices = logical(xindices);
end
dhatsGLSAllSubjects = zeros(nSubject,nCov,nVertex);
dhatsOLSAllSubjects = zeros(nSubject,nQ,nVertex);
varResAllSubjects = zeros(nSubject,nVertex);
varInnovAllSubjects = zeros(nSubject,nVertex);
phiAllSubjects = zeros(nSubject,3,nVertex);
rhoAllSubjects = zeros(nSubject,3,nVertex);

for n = 1:nSubject
    subjDesignMat = squeeze(allsubjdesignmat(n,:,:));
    invM = squeeze(invMAllSubjects(n,:,:));
    parfor v=1:nVertex    
        [betasTemp,phiTemp,sigmasqTemp,covTemp,rhoTemp] = fitarpOLSYW(squeeze(fdata(n,v,:)),subjDesignMat,3,invM);
        dhatsOLSAllSubjects(n,:,v)=betasTemp(xindices);
        invTemp = covTemp\subjDesignMat;
        intTemp = subjDesignMat'*invTemp;
        invInvTemp = intTemp\invTemp';
        dhatsGLSAllSubjects(n,:,v) = invInvTemp*squeeze(fdata(n,v,:));
        varResAllSubjects(n,v) = covTemp(1,1);
        varInnovAllSubjects(n,v) = sigmasqTemp;
        phiAllSubjects(n,:,v) = phiTemp;
        rhoAllSubjects(n,:,v) = rhoTemp;
    end
end

nbw = length(bandwidth);
mseCortex = zeros(nSubject,nbw);
if nbw>1
    for n=1:nSubject
     [tempMSE1] = gcv_temporalparameters(varResAllSubjects(n,:), squeeze(phiAllSubjects(n,:,:)),bandwidth,distmat);
     mseCortex(n,:) = tempMSE1;
    end
    mseSum = sum(mseCortex);
    h_ob = bandwidth(mseSum == min(mseSum));
else
    h_ob = bandwidth;
end

%% smooth temporal variance parameters:
smVarQAllSubjects = zeros(nSubject,nCov,nCov,nVertex);
smVarResAllSubjects = zeros(nSubject,nVertex);
smPhiAllSubjects = zeros(nSubject,3,nVertex);
smVarInnovAllSubjects = zeros(nSubject,nVertex);
    
parfor v=1:nVertex    
    dist = distmat(:,v);
    index = (dist<=h_ob);
    wts = biweight(dist(index)/h_ob)/h_ob;
    wts = wts/sum(wts);

    phiNeighbors = phiAllSubjects(:,:,index);
    rhoNeighbors = rhoAllSubjects(:,:,index);
    sigmasqTall = varResAllSubjects(:,index)*wts;
    
    phiTall1= squeeze(phiNeighbors(:,1,:))*wts;
    phiTall2= squeeze(phiNeighbors(:,2,:))*wts;
    phiTall3= squeeze(phiNeighbors(:,3,:))*wts;
    phiTall = [phiTall1, phiTall2, phiTall3];
    rhoTall1 = squeeze(rhoNeighbors(:,1,:))*wts;
    rhoTall2 = squeeze(rhoNeighbors(:,2,:))*wts;
    rhoTall3 = squeeze(rhoNeighbors(:,3,:))*wts;
    rhoTall = [rhoTall1, rhoTall2, rhoTall3];
    
    smVarResAllSubjects(:,v) = sigmasqTall;
    smPhiAllSubjects(:,:,v) = phiTall;
    smVarInnovAllSubjects(:,v) = varInnovAllSubjects(:,index)*wts;
    p=3;
    for n = 1:nSubject
        subjDesignMat = squeeze(allsubjdesignmat(n,:,:));
        cap = (subjDesignMat'*subjDesignMat)\subjDesignMat';   
        corrTemp = zeros(nTR,1);
        corrTemp(1) = 1;
        corrTemp(2:(p+1)) = rhoTall(n,:);  
        for t=(p+2):nTR
             corrTemp(t) = phiTall(n,:)*corrTemp((t-1):-1:(t-p));
        end
        covTemp = sigmasqTall(n)*corrTemp;
        covmat = toeplitz(covTemp);
        covmat = kron(eye(2),covmat);
        varQ = cap*covmat*cap';
        % this results in a matrix that is not quite symmetric; symmetrize:
        varQ = (varQ + varQ')/2; 
        smVarQAllSubjects(n,:,:,v) = varQ;
    end
end

dTemp = dhatsGLSAllSubjects(:,xindices,:);
vTemp = smVarQAllSubjects(:,xindices,xindices,:);
clear('dhatsAllSubjects','varQAllSubjects','smVarQAllSubjects');

outMUMM.betaHatVertex = squeeze(mean(dTemp,1))';
outMUMM.varPop = squeeze(var(dTemp))'./nSubject;
outMUMM.t_betaHatVertex = outMUMM.betaHatVertex./sqrt(outMUMM.varPop);
outMUMM.dhats = dTemp;
outMUMM.contrast_betaHatVertex = outMUMM.betaHatVertex*contrastmat';

nContrast = size(contrastmat,1);
subjResultsContrast = zeros(nSubject,nVertex,nContrast);
for n=1:nSubject
    subjResultsContrast(n,:,:) = contrastmat*squeeze(dTemp(n,:,:));
end

varContrast = squeeze(var(subjResultsContrast)./nSubject);
outMUMM.t_contrast_betaHatVertex = outMUMM.contrast_betaHatVertex./sqrt(varContrast)';

outSTMM = estSTMM1(dhatsOLSAllSubjects,vTemp,distmat,1,contrastmat,resultsDir,savecovario);
outSTMM.varResAllSubjects = smVarResAllSubjects;
outSTMM.varInnovAllSubjects = smVarInnovAllSubjects;
outSTMM.phiAllSubjects = smPhiAllSubjects;
outSTMM.bw_mse = mseCortex;
end
