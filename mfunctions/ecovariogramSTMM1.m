function [ecovariogram,stddev,count,avedist] = ecovariogramSTMM1(subjdata,meandata,distMat,hbin)

% Benjamin Risk
% 8 July 2015
% Create empirical COVARIOGRAM for two dimensional data. Uses
% d_{nvq}*d_{nv'q} - \bar{d}_{\cdot v q} \bar{d}_{\cdot v' q}
% Input:
%   subjdata    vectorized spatial data for a single subject in 
%               column-major order
%   meandata    vectorized spatial data of mean across subjects in
%               column-major order
%   distMat     distance matrix with dimensions nrow*ncol x nrow*ncol in
%               column major order. 
%   hbin        mesh for binning observations (equal to the number of bins
%               plus one such that the last element represents the upper 
%               bound at which distances will be calculated)
% Output:
%   ecovariogram  the empirical covariogram
%   stddev      standard deviation for each value of h, where h is the distance
%				NOTE: These standard deviations are biased due to dependence between
%				pairs of points and are not used.
%   count       the number of data points used for each of the distances in
%               hbin
%   avedist     the average distance between points for all pairs that fall
%               within a bin.
if hbin(1)<0
    warning('My Warning: hbin(1) should be greater than or equal to zero');
end
hsize = length(hbin);
ecovariogram = zeros(hsize-1,1);
stddev = ecovariogram;
nRow = length(subjdata);
count = zeros(hsize-1,1);
avedist = count;
for j=1:(hsize-1)
   [iIndex,jIndex] = find(distMat>hbin(j) & distMat <=hbin(j+1));
   allIndex = nRow*(jIndex-1)+iIndex;
   distances = distMat(allIndex);
   avedist(j) = mean(distances);
   covtemp = subjdata(iIndex).*subjdata(jIndex)-meandata(iIndex).*meandata(jIndex);
   count(j) = length(covtemp);
   ecovariogram(j) = mean(covtemp);
   stddev(j) = std(covtemp,1); %the flag "1" uses n as the divisor; each observation
    % appears twice but this standard deviation is equal to the std calculated with
    % each observation appearing once.
 end
end
