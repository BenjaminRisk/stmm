function [hvec,mse] = cv_5fold_temporalparameters(varVec, phiArray, hvec, lat, long)
%Use cross-validation to select a single bandwidth for all vertices
%Input:
% varArray:    1 x V vector of residuals (V locations)
% phiArray:    3 x V matrix of AR3 parameters
% hvec: vector of bandwidths to evaluate
% lat: V x 1 vector of latitudes
% long: V x 1 vector of longitudes

%Output:
% hvec: bandwidths at which MSE is calculated
% mse: vector of errors of length hvec

nbw = length(hvec);
% 
% if nargin==5
%     ROI = 1:length(varVec);
%     nVertexROI = length(varVec);
% end

arArray = [(varVec-mean(varVec))./std(varVec); (phiArray(1,:)-mean(phiArray(1,:)))./std(phiArray(1,:)); ...
    (phiArray(2,:)-mean(phiArray(2,:)))./std(phiArray(2,:));(phiArray(3,:)-mean(phiArray(3,:)))./std(phiArray(3,:))];
K = 5;

% use cross validation with FAMILY as the unit of observation
nVertex = length(varVec);
cvpart = cvpartition(nVertex,'Kfold',K);
mse = zeros(nbw,K,max(cvpart.TestSize));

for k=1:K
    tIndex = find(cvpart.test(k));
    for v=1:cvpart.TestSize(k)
        vTemp = tIndex(v);
        dist = double(mygreatcirc(lat(vTemp),long(vTemp),lat,long));    
        for j=1:nbw
            h = hvec(j);
            index = (dist<=h & cvpart.training(k));
            subAR = arArray(:,index);
            wts = biweight(dist(index)/h)/h;
            wts = wts/sum(wts);
            tempError=0;
            for n=1:4
                tempAR = subAR(n,:)*wts;
                tempError = tempError+(arArray(n,vTemp) - tempAR)^2;
            end
            mse(j,k,v) =  tempError/cvpart.TestSize(k);
        end
    end
end
mse = sum(mse,3);
mse = sum(mse,2);
end

