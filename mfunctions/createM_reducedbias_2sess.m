function [M] = createM_reducedbias_2sess(Hatmat,L)
% createM_reducedbias_2sess
% Create the M matrix used in calculating the reduced bias sample autocorrelations.
% This function modifies createM_reducedbias to accomodate two temporally
% independent sessions with an equal number of timepoints. 
%   INPUT:
%       Hatmat - the hat matrix used in the OLS 
%       p      - the number of sample autocorrelations to calculate
%   OUTPUT:
%       the "M" matrix in Appendix A of Worsley 2002

% NOTE: This code could be made more efficient
nTime = size(Hatmat,1)/2;
R = eye(nTime*2) - Hatmat;
M = zeros(L+1,L+1);
    for l = 0:L
        if l==0
            Dl = eye(2*nTime);
        else
            temp = zeros(nTime,1);
            temp(l+1)=1;
            Dl = toeplitz(zeros(nTime,1),temp);
            Dl = kron(eye(2),Dl);
        end
        for j=0:L
            if j==0
                M(l+1,j+1) = trace(R*Dl);
            else
                temp = zeros(nTime,1);
                temp(j+1) = 1;
                DjtDj = kron(eye(2),toeplitz(temp));
                M(l+1,j+1) = trace(R*Dl*R*DjtDj);
            end
        end
    end
end

