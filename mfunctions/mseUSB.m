function [mse] = mseUSB(ddata)
% Benjamin Risk
% 21 April 2015
% Calculate sum of squares decomposition for the second-level of the STMM
% 
% Input:
%   ddata: N x V matrix
%
% Output:
%   An object "mse" with three fields:
%   mse.U: mean square from the location effects
%   mse.S: mean square from the subject effects
%   mse.B: mean square from the location-subject interaction.
nParcelm = size(ddata,2);
nSubject = size(ddata,1);
%if nSubject > nParcelm
%    warning('NOTE: There are more subjects than vertices. Make sure that subjects comprise rows and vertices comprise columns');
%end
ddataMean = mean(ddata(:));
ddataU = mean(ddata,1) - ddataMean;
ddataS = mean(ddata,2) - ddataMean;
%here we subtract the overall mean because we are using the centered effects
ddataB = ddata - ones(nSubject,1)*ddataU - ddataS*ones(1,nParcelm) - ddataMean; 
mse.U = nSubject*sum(ddataU(:).^2)/(nParcelm-1);
mse.S = nParcelm*sum(ddataS(:).^2)/(nSubject-1);
mse.B = sum(ddataB(:).^2)/((nParcelm-1)*(nSubject-1));
end
