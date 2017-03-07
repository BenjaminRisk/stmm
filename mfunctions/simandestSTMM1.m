function [out] = simandestSTMM1(nsim,allSubjDesignMat,betaFixed,mySsigma,...
    myBsigma,myBtheta,myPhiarray,mySigmasq,hvec,parcelDist,dhatParcel,myInvMAllSubjects)
%simandestSTMM1 simulate and estimate STMM1
%   Detailed explanation goes here
meanSTMM = zeros(18,nsim);
mseSTMM = zeros(15,nsim);
biasDhatSTMM = zeros(3,nsim);
nParcelm = size(dhatParcel,3);
meanMUMM = zeros(9,nsim);
mseMUMM = zeros(6,nsim);
biasDhatMUMM = zeros(3,nsim);

if numel(betaFixed)==2
    betaFixed = ones(nParcelm,1)*betaFixed;
end

for r=1:nsim
    [fData,sData] = simSTMM1(allSubjDesignMat,betaFixed,mySsigma,myBsigma,...
        myBtheta,myPhiarray,mySigmasq,parcelDist,dhatParcel);
    
    [stmmresults,mummresults] = fitSTMM1(fData,parcelDist,hvec,allSubjDesignMat,...
        [1 0 0 1 zeros(1,24)],[1 -1],myInvMAllSubjects,['~/Dropbox/SpatioTemporalFMRI_FixedVertex/TempResults/simcovario' num2str(r)]);

    % compile statistics for MUMM:
    meanMUMM(1:2,r) = mean(mummresults.betaHatVertex);
    meanMUMM(3,r) = mean(mummresults.contrast_betaHatVertex);

    mseMUMM(1,r) = mean((mummresults.betaHatVertex(:,1)-betaFixed(:,1)).^2);
    mseMUMM(2,r) = mean((mummresults.betaHatVertex(:,2)-betaFixed(:,2)).^2);
    meanMUMM(3,r) = mean(mummresults.contrast_betaHatVertex);
    tempbias = mummresults.contrast_betaHatVertex - betaFixed(:,1) + betaFixed(:,2);
    mseMUMM(3,r) = mean(tempbias(:).^2);
    
    meanMUMM(7:8,r) = mean(abs(mummresults.t_betaHatVertex)>2.04523);
    meanMUMM(9,r) = mean(abs(mummresults.t_contrast_betaHatVertex)>2.04523);
    
    difftemp = permute(mummresults.dhats,[1,3,2])-sData;
    d1 = squeeze(difftemp(:,:,1));
    d2 = squeeze(difftemp(:,:,2));
    biasDhatMUMM(1,r) = mean(d1(:));
    biasDhatMUMM(2,r) = mean(d2(:));
    biasDhatMUMM(3,r) = biasDhatMUMM(1,r) - biasDhatMUMM(2,r);
    
    meanMUMM(4,r) = mean(d1(:));
    meanMUMM(5,r) = mean(d2(:));
    meanMUMM(6,r) = mean(d1(:)-d2(:));
    mseMUMM(4,r) = mean(d1(:).^2);
    mseMUMM(5,r) = mean(d2(:).^2);
    mseMUMM(6,r) = mean((d1(:)-d2(:)).^2);
    
    % compile statistics for STMM:
    meanSTMM(1:2,r) = mean(stmmresults.betaHatVertex);
    meanSTMM(3,r) = mean(stmmresults.contrast_betaHatVertex);
    
    mseSTMM(1,r) = mean((stmmresults.betaHatVertex(:,1) - betaFixed(:,1)).^2);
    mseSTMM(2,r) = mean((stmmresults.betaHatVertex(:,2) - betaFixed(:,2)).^2);
    mseSTMM(3,r) = mean((stmmresults.contrast_betaHatVertex - betaFixed(:,1) + betaFixed(:,2)).^2);
    
    meanSTMM(4:5,r) = stmmresults.Shat;
    meanSTMM(6:7,r) = stmmresults.Bhat;
    meanSTMM(8:9,r) = [stmmresults.bTheta(1,2),stmmresults.bTheta(2,2)];

    meanSTMM(10,r) = mean(stmmresults.varInnovAllSubjects(:));
    mseSTMM(10,r) = mean((stmmresults.varInnovAllSubjects(:) - mySigmasq(:)).^2);
       
    temp1 = stmmresults.phiAllSubjects(:,1,:);
    if numel(myPhiarray)>3
        phi1 = myPhiarray(:,1,:);
        phi2 = myPhiarray(:,2,:);
        phi3 = myPhiarray(:,3,:);
    else
        phi1 = myPhiarray(1);
        phi2 = myPhiarray(2);
        phi3 = myPhiarray(3);
    end
    meanSTMM(11,r) = mean(temp1(:));
    mseSTMM(11,r) = mean((temp1(:)-phi1(:)).^2);
    
    temp2 = stmmresults.phiAllSubjects(:,2,:);
    meanSTMM(12,r) = mean(temp2(:));
    mseSTMM(12,r) = mean((temp2(:) - phi2(:)).^2);
    
    temp3 = stmmresults.phiAllSubjects(:,3,:);
    meanSTMM(13,r) = mean(temp3(:));
    mseSTMM(13,r) = mean((temp3(:)-phi3(:)).^2);
    
    tempSTMM = permute(stmmresults.eBLUPsBetaHat,[1,3,2]) - sData;
    temp1 = squeeze(tempSTMM(:,1,:));
    temp2 = squeeze(tempSTMM(:,2,:));
    temp3 = temp1-temp2;
    biasDhatSTMM(1,r) = mean(temp1(:));
    biasDhatSTMM(2,r) = mean(temp2(:));
    biasDhatSTMM(3,r) = mean(temp3(:));
    
    temp4 = stmmresults.eBLUPsBetaHat(:,1,:);
    temp5 = stmmresults.eBLUPsBetaHat(:,2,:);
    temp6 = temp4 - temp5;
    meanSTMM(14,r) = mean(temp4(:));
    meanSTMM(15,r) = mean(temp5(:));
    meanSTMM(16,r) = mean(temp6(:));
    
    mseSTMM(14,r) = mean(temp1(:).^2);
    mseSTMM(15,r) = mean(temp2(:).^2); 
    mseSTMM(16,r) = mean(temp3(:).^2);
    
    meanSTMM(17:18,r) = mean(abs(stmmresults.t_betaHatVertex)>1.96);
    meanSTMM(19,r) = mean(abs(stmmresults.t_contrast_betaHatVertex>1.96));
end

meanEstimatesMUMM = mean(meanMUMM,2);
BiassqEstimatesMUMM = zeros(6,1);
BiassqEstimatesMUMM(1:2) = (meanEstimatesMUMM(1:2)' - mean(betaFixed)).^2;
BiassqEstimatesMUMM(3) = (meanEstimatesMUMM(3) - mean(betaFixed*[1;-1]))^2;
BiassqEstimatesMUMM(4:6) = mean(biasDhatMUMM,2).^2;
MSEestimatesMUMM = mean(mseMUMM,2);

BiassqEstimatesSTMM = zeros(16,1);    
meanEstimatesSTMM = mean(meanSTMM,2);
BiassqEstimatesSTMM(1:2) = (meanEstimatesSTMM(1:2)' - mean(betaFixed)).^2;
BiassqEstimatesSTMM(3) = (meanEstimatesSTMM(3) - mean(betaFixed*[1;-1]))^2;

MSEestimatesSTMM = mean(mseSTMM,2);
BiassqEstimatesSTMM(14:16) = mean(biasDhatSTMM,2).^2;

trueValues = [mean(betaFixed),mean(betaFixed*[1;-1]),mySsigma,myBsigma,myBtheta,mean(mySigmasq(:)),mean(phi1(:)),mean(phi2(:)),mean(phi3(:)),...
    mean(betaFixed),mean(betaFixed*[1;-1])]';

for j=4:9
    BiassqEstimatesSTMM(j) = (meanEstimatesSTMM(j) - trueValues(j))^2;
    MSEestimatesSTMM(j) = mean((meanSTMM(j,:) - trueValues(j)).^2);
end

for j = 10:13
    BiassqEstimatesSTMM(j) = (meanEstimatesSTMM(j) - trueValues(j))^2;
end

out.meanEstimatesSTMM = meanEstimatesSTMM(1:16);
out.rejectRatesSTMM = meanEstimatesSTMM(17:19);
out.VarEstimatesSTMM = MSEestimatesSTMM(1:16) - BiassqEstimatesSTMM;
out.BiassqEstimatesSTMM = BiassqEstimatesSTMM;
out.MSEestimatesSTMM = MSEestimatesSTMM;
out.seMSEestimatesSTMM = sqrt(var(mseSTMM,[],2)/nsim);
out.seMeanEstimatesSTMM = sqrt(var(meanSTMM,[],2)/nsim);

out.meanEstimatesMUMM = meanEstimatesMUMM(1:6); 
out.rejectRatesMUMM = meanEstimatesMUMM(7:9);
out.VarEstimatesMUMM = MSEestimatesMUMM - BiassqEstimatesMUMM;
out.BiassqEstimatesMUMM = BiassqEstimatesMUMM;
out.MSEestimatesMUMM = MSEestimatesMUMM;
out.seMSEestimatesMUMM = sqrt(var(mseMUMM,[],2)/nsim);
out.seMeanEstimatesMUMM = sqrt(var(meanMUMM,[],2)/nsim);
out.trueValues = trueValues;
end