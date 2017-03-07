function [ecovariogram,stderr,count,avedist,cloud,distCloud] = ecovariogram(data,distMat,hbin,retCloud)

% Benjamin Risk
% 25 March 2015
% Create empirical COVARIOGRAM for two dimensional data. Uses d_{nvq}*d_{nv'q}
% Input:
%   data        vectorized two-dimensional data in column-major order
%   distMat     distance matrix with dimensions nrow*ncol x nrow*ncol in
%               column major order. Can only be a distance matrix where
%               zeros correspond to "too far" as well as 0 distance, which
%               corresponds to the my old function evariogramDistThresh
%   hbin        mesh for binning observations (equal to the number of bins
%               plus one such that the last element represents the upper 
%               bound at which distances will be calculated)
%   retCloud    return each product that is used in the covariogram
% Output:
%   ecovariogram  the empirical covariogram
%   stderr      standard error for each value of h, where h is the distance
%   count       the number of data points used for each of the distances in
%               hbin
%   avedist     the average distance between points for all pairs that fall
%               within a bin.
%   cloud       vector of all elements included in the covariogram
%   distCloud   vector all distances corresponding to cloud
if hbin(1)<0
    warning('My Warning: hbin(1) should be greater than or equal to zero');
end
hsize = length(hbin);
ecovariogram = zeros(hsize-1,1);
stderr = ecovariogram;
nRow = length(data);
count = zeros(hsize-1,1);
avedist = count;
%break ties by randomly adding tiny bit of noise:
%distMat(distMat~=0) = randn(nnz,1)/(tol*1000);
cloud=0; 
distCloud=0;
for j=1:(hsize-1)
   [iIndex,jIndex] = find(distMat>hbin(j) & distMat <=hbin(j+1));
   allIndex = nRow*(jIndex-1)+iIndex;
   distances = distMat(allIndex);
   avedist(j) = mean(distances);
   % note: mean(data(iIndex))==mean(data(jIndex))
   % lmean = mean(data(iIndex));
   % covtemp = (data(iIndex)-lmean).*(data(jIndex)-lmean);
   covtemp = data(iIndex).*data(jIndex);
   count(j) = length(covtemp);
   if retCloud
       cloud = [cloud;covtemp];
       lIndex = nRow*(jIndex-1)+iIndex;
       distCloud = [distCloud;distMat(lIndex)];
   end
   ecovariogram(j) = mean(covtemp);
   stderr(j) = std(covtemp)/(sqrt(length(covtemp)/2)); %each obs appears twice
 end
end
