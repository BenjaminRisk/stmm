function [out] = combineoutsim(out1,nsim1,out2,nsim2)
nsim12 = nsim1+nsim2;
out.meanEstimatesSTMM = (nsim1*out1.meanEstimatesSTMM+nsim2*out2.meanEstimatesSTMM)./nsim12;
out.rejectRatesSTMM = (nsim1*out1.rejectRatesSTMM+nsim2*out2.rejectRatesSTMM)./nsim12;
out.MSEestimatesSTMM = (nsim1*out1.MSEestimatesSTMM+nsim2*out2.MSEestimatesSTMM)./nsim12;

out.rejectRatesMUMM = (nsim1*out1.rejectRatesMUMM+nsim2*out2.rejectRatesMUMM)./nsim12;
out.MSEestimatesMUMM = (nsim1*out1.MSEestimatesMUMM+nsim2*out2.MSEestimatesMUMM)./nsim12;

end

