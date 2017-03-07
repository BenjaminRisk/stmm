function [bf,time] = hrf(TR,nBin,p)
%Generate the hrf as the difference between two gamma pdfs. This
%corresponds to the canonical HRF in SPM. 
%DETAILS:
%   Returns a matrix with the first column corresponding to the HRF, the
%   second column corresponding to the derivative of the HRF with respect
%   to the "onset" parameter (shift in when the HRF begins), and the third 
%   column corresponding to the derivative with respect to the dispersion parameter. 
%   The derivatives are approximated using delta=0.01.
%----------------------------------------------------------------------
% Input:
% 
% TR   - time between scans
%
% nBin - number of points at which to evaluate hrf for each volume. That
% is, the hrf function will be evaluated at the points 0:(TR/nBin):(length
% of HRF). This parameter (together with the TR) will determine the accuracy of
% the discrete approximation to the integration in the convolution of the HRF 
% with stimulus onsets and durations.
%
% p    - parameters of the response function (two Gamma functions).
%The following is copied exactly from <spm_hrf.m> --->
%        p(1) - delay of response (relative to onset)          6
%        p(2) - delay of undershoot (relative to onset)       16
%        p(3) - dispersion of response                         1
%        p(4) - dispersion of undershoot                       1
%        p(5) - ratio of response to undershoot                6
%        p(6) - onset (seconds)                                0
%        p(7) - length of kernel (seconds)                    32
% <---- end copy
% Brisk NOTE: p(7) is NOT a parameter of the HRF but rather
% the time at which the hrf function is no longer calculated.
% The support of the function is [0,Inf), so this actually corresponds to
% the points at which we truncate the HRF.
%
%---------------------------------------------------------------------
% Output:
% bf   - A matrix with three columns corresponding to the hrf, its
% derivative with respect to the onset parameter (p(6)), and its
% derivative with respect to the dispersion of response (p(3))
if nargin <2
    nBin = 16;
end
if nargin<3
    p   = [6 16 1 1 6 0 32];
end
mesh = TR/nBin;
time = 0:mesh:p(7);
bf = subhrf(time,p);

% Approximate derivative with respect to onset parameter.
% spm12:spm_hrf uses p(6)/mesh for the onset parameter; I simplify and
% let p(6) be the time in seconds of the hrf delay.
% spm reverses the signs on the dispersal derivative, i.e.,
% (y(x) - y(x+delta))/delta; I correct this and use (y(x+delta) - y(x))/delta. 
delta=0.01;
p(6) = p(6)+delta;
dbf_dp6 = (subhrf(time,p) - bf)./delta;
p(6) = p(6) - delta;

% Approximate derivative wrt to dispersion parameter:
p(3) = p(3) + delta;
dbf_dp3 = (subhrf(time,p) - bf)./delta;

bf = [bf' dbf_dp6' dbf_dp3'];
end

%---------------------------------------
% subfunction to calculate the HRF only:
function bf = subhrf(x,p)
% Note that the function spm12:spm_Gpdf.m parameterizes the gamma pdf with the
% shape and rate parameter, whereas stats:gampdf parameterizes with the shape and
% scale parameter. This function tries to use the stats package but if not
% found uses spm.

% spm uses p(6) in the "temporal derivative" calculation. Here, it 
% DELAYS the start of the HRF p(6). spm subtracts mesh/p(6),
% but here we simplify and simply add p(6).
x = x + p(6); 

try
    bf = gampdf(x,p(1)/p(3),p(3)) - gampdf(x,p(2)/p(4),p(4))/p(5);
catch
    try
    bf = spm_Gpdf(x,p(1)/p(3),1/p(3)) - spm_Gpdf(x,p(2)/p(4),1/p(4))/p(5);        
    catch
        error('stmm:Undefined function -- addpath to either stats:gampdf or spm12:spm_Gpdf to calculate gamma pdf');
     end
end
end