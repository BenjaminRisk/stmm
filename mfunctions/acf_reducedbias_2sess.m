function [acf,var0] = acf_reducedbias_2sess(resids,Hatmat,L,invM)
%acf_reducedbias Estimate reduced bias sample autocorrelations for p lags 
%from residuals from OLS estimation. This code modifies acf_reducedbias.m
%to accomodate two temporally independent sessions.
% INPUT:
%   resids - Tx1 residuals from OLS estimation.
%   Hatmat - TxT hat matrix: X(X'X)^{-1} X'; not actually used if M is
%   supplied.
%   p      - number of ACFs to calculate. 
%   invM   - optional argument equal to the inverse of the matrix with
%            entries equal to M_{lj} =  tr((I-H)D_l(I-H)(D_j + D_j'))
%
% OUTPUT:
%   acf:   - reduced bias autocorrelation function evaluated at lags 0 to L
%   var0:  - reduced bias estimate of the variance
T = size(resids,1);
if T==1
    error('resids should be a column vector');
end
a = zeros(L+1,1);
resids1 = resids(1:(T/2));
resids2 = resids(((T/2)+1):end);
for l=0:L
    a(l+1) = resids1((l+1):end)'*resids1(1:(end-l)) + resids2((l+1):end)'*resids2(1:(end-l));
end

if nargin<4
   M = createM_reducedbias_2sess(Hatmat,L);
   invM = inv(M);
end
intacf = invM*a;
acf = intacf./intacf(1);
var0 = intacf(1);
end
        