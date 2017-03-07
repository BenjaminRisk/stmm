function [rt] = simAR3_2sess(n, phivec, sigmasq)
%SIMAR3 simulate AR3 model with two independent n/2 time series
% 18 May 2015: Added burnin = 100
        burnin=100;
        innov1 = sqrt(sigmasq)*randn(n/2+burnin,1);
        innov2 = sqrt(sigmasq)*randn(n/2+burnin,1);
        rt1 = zeros(n/2+burnin,1);
        rt2 = rt1;
        rt1(1) = innov1(1);
        rt1(2) = phivec(1)*rt1(1)+innov1(2);
        rt1(3) = phivec(1)*rt1(2)+phivec(2)*rt1(1)+innov1(3);

        rt2(1) = innov2(1);
        rt2(2) = phivec(1)*rt2(1)+innov2(2);
        rt2(3) = phivec(1)*rt2(2)+phivec(2)*rt2(1)+innov2(3);

        for i = 4:(n/2+burnin)
            rt1(i) = phivec*rt1((i-1):-1:(i-3))+innov1(i);
            rt2(i) = phivec*rt2((i-1):-1:(i-3))+innov2(i);
        end
        rt = [rt1((burnin+1):(n/2+burnin));rt2((burnin+1):(n/2+burnin))];
end
