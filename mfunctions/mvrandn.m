function [z] = mvrandn(n,covMat)
%MVRANDN generate n iid realizations of a multivariate normal random vector 
%with covariance matrix covMat
[eigvec,eigval] = eig(covMat);
if sum(diag(eigval)<=0)
    error('Covariance matrix not positive definite');
end
V = size(covMat,1);
sqrtCov = eigvec*sqrt(eigval)*eigvec';
z = sqrtCov*randn(V,n);
end

