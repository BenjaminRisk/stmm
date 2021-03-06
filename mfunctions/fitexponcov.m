function [thetahat,f,fval,eflag] = fitexponcov(h,ecovario,lb,ub,thetainit)
%Fit exponential covariogram to the empirical covariogram
%
% INPUT:
%   h: distances at which ecovario was cal
%   ecovario: empirical covariogram
%
% OUTPUT:
%   thetahat(1): variance parameter of exponential covariogram
%   thetahat(2): dependence parameter
%   thetahat(3): baseline dependence -- constant for all voxels
%   f:           covariogram function handle: f(h,thetahat)
if size(h,1)<size(h,2)
    h = h';
end
if size(ecovario,1)<size(ecovario,2)
    ecovario = ecovario';
end

f1 = @(x,ecovario,h) sum((ecovario - x(1)*exp(-x(2)*h)-x(3)).^2);
f2 = @(x)f1(x,ecovario,h);

% default number of MaxFunEvals = 100*numberofvariables
options = optimset('MaxIter',1000,'MaxFunEvals',20000,'Display','off','TolFun',1e-5,'LargeScale','off');
%options = optimset('Algorithm','interior-point','MaxFunEvals',10000);
[thetahat,fval,eflag] = fmincon(f2,thetainit,[],[],[],[],lb,ub,[],options);

f = @(h) thetahat(1)*exp(-thetahat(2)*h)+thetahat(3);
%fittedh = f(h,thetahat);
end
