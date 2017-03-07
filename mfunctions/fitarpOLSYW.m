function [betahat,phihat,sigmasqhat,covmat,rhohat] = fitarpOLSYW(y,xData,p,invM)
% Benjamin Risk 
% March 2015
%This function estimates betas using OLS and estimates the 
% AR parameters using the Yule-Walker equations on the OLS residuals.
% Estimation of the AR parameters accounts for two independent session
% NOTE: covmat assumes two independent sessions
% NOTE: Assumes there are two sessions with an equal number of timepoints
%
% Input:    y - timeseries
%           xData - covariates
%           p - AR order to model
%           invM - inverse of M matrix used in calculating reduced biased
%               sample autocorrelations. The calculation of M is slow, but only
%               needs to be done once per subject, whereas this function is called for
%               every vertex. M can be calculated using <createM_reducedbias.m> 

% Output:
%           betahat - covariate coefficients
%           phihat - AR coefficients, i.e., PACF at lags 1:L
%           sigmasqhat - estimate of the variance of the INNOVATION process
%           covmat - T x T unconditional covariance matrix
%           rhohat - unconditional correlations, i.e., ACF at lags 1:L
% Edits: 19 May 2015
% added code to calculate reduced-bias autocorrelations in the residuals.
% uses the function "acf_reducedbias.m"
betahat = xData\y;
res = y - xData*betahat;
sessLength = length(y)/2;

% calculate autocorrelation for two sessions:
L = size(invM,1)-1;
[rhohat,sigmaMar] = acf_reducedbias_2sess(res,false,L,invM);
rhohat = rhohat(2:(p+1));
Rhat = eye(p);
    for i=1:p
        for  j=1:p
         if i ~= j
            Rhat(i,j) = rhohat(abs(i-j)); 
         end
        end
    end    
phihat = Rhat\rhohat;

if(any(abs(roots([-phihat(p:-1:1);1]))<1))
    warning('There exists a root inside the unit circle')
end
 
sigmasqhat = sigmaMar*(1-rhohat'*phihat);

corrTemp = zeros(sessLength,1);
corrTemp(1) = 1;
corrTemp(2:(p+1)) = rhohat;
%     t = p+1;
%     check = 0.5;
%     while abs(check)>1e-06 && t<sessLength
%        t=t+1;
%        check = rho'*corrTemp((t-1):-1:(t-p));
%        corrTemp(t) = check;
%     end
     for t=(p+2):sessLength
         corrTemp(t) = phihat'*corrTemp((t-1):-1:(t-p));
     end
    covTemp = sigmaMar*corrTemp;
    covmat = toeplitz(covTemp);
    covmat = kron(eye(2),covmat);
end
