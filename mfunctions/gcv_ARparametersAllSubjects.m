function [mseMat,hvecMat,optbw] = gcv_ARparametersAllSubjects(hvec,lat,long,varResAll,phiAll)
[nSubject,nVertex] = size(varResAll);
hvecMat = zeros(nSubject,length(hvec));
mseMat = zeros(nSubject,length(hvec));
distMat = zeros(nVertex,nVertex);
parfor v=1:nVertex
    distMat(v,:) = distance(lat(v),long(v),lat,long);
    %distMatR(v,:) = mygreatcirc(latR(v),longR(v),latR,longR);
end

for n=1:nSubject
    [tempMSE1] = gcv_temporalparameters(varResAll(n,:), squeeze(phiAll(n,:,:)),hvec, distMat);
    hvecMat(n,:) = hvec;
    mseMat(n,:) = tempMSE1;
end

[~,mS] = min(mseMat,[],2);
optbw = zeros(nSubject,1);
for n=1:nSubject
    optbw(n) = hvecMat(n,mS(n));
end
optbw = mean(optbw);
