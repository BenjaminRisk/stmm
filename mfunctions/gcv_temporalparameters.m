function [mse] = gcv_temporalparameters(varVec,phiArray,hvec,distmat)
% Use gcv approximation to leave-one-out cross-validation to estimate 
% spatial smoothing bandwidth for AR parameters
%INPUT:
% varVec:    1 x V vector of marginal variance of the AR process
% phiArray:    3 x V matrix of AR3 parameters
% hvec: vector of bandwidths to evaluate
% distmat: VxV dense matrix of distances
%
%OUTPUT:
% mse: 1 x length(hvec)
nVertex = size(distmat,1);
fortrace = zeros(nVertex,1);

% standardize four parameters so they are weighted evenly:
arArray = [(varVec-mean(varVec))./std(varVec); (phiArray(1,:)-mean(phiArray(1,:)))./std(phiArray(1,:)); ...
    (phiArray(2,:)-mean(phiArray(2,:)))./std(phiArray(2,:));(phiArray(3,:)-mean(phiArray(3,:)))./std(phiArray(3,:))];

mse = zeros(length(hvec),1);
mseTemp = zeros(nVertex,1);
if hvec(1)==0
    hvec(1) = 1e-4;
end

for j = 1:length(hvec)
    for v=1:nVertex
           h = hvec(j);
           index = (distmat(:,v)<=h);
           subAR = arArray(:,index);
           wts = biweight(distmat(index,v)/h)/h;
           wts = wts/sum(wts);
           fortrace(v) = max(wts); %obtain diagonal element of "S" matrix,
                                    % eg 7.5.2, p.244, ESLII
           tempError=0;
           for n=1:4
                tempAR = subAR(n,:)*wts;
                tempError = tempError+(arArray(n,v) - tempAR)^2;
            end
           mseTemp(v) =  tempError;
    end
    mse(j) = sum(mseTemp)/(1-sum(fortrace)/nVertex)^2;
 end
mse = mse./nVertex;
end
