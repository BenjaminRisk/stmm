function fbasis = createLinSpline(TR,nTR,nKnot)
%create a basis for piecewise linear spline with equally spaced knots
% the first column is all ones
% the second column is the time
% columns 3 to nKnots are shifted time
%INPUT:
% TR: time between scans
% nTR: number of scans
% nKnots: number of knots in the spline; equaly to the rank of the basis.
%
%OUTPUT:
% Rank-nKnot piecewise-linear basis
 time = 0:TR:((nTR-1)*TR);
 fbasis = [ones(nTR,1),time',zeros(nTR,(nKnot-2))];
 breakpts = linspace(0,time(end),nKnot);
 for r=3:nKnot
    fbasis(:,r) = (time-breakpts(r-1)).*(time>breakpts(r-1));
 end
end

