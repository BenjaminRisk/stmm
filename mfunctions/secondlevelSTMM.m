function [] = secondlevelSTMM(contrastmat,resultsDir,lat,long,hemisphere)

if strcmp(hemisphere,'CORTEX_RIGHT')
    load('./supportingdatafiles/matrix_GordonRev2_CII_CORTEX_RIGHT.mat')
    matrix_Gordon = matrix_GordonRev2_CII_CORTEX_RIGHT;
elseif strcmp(hemisphere,'CORTEX_LEFT')
    load('./supportingdatafiles/matrix_GordonRev2_CII_CORTEX_LEFT.mat')
    matrix_Gordon = matrix_GordonRev2_CII_CORTEX_LEFT;
else 
    error('Specify hemisphere equals CORTEX_RIGHT or CORTEX_LEFT')
end

[nVertex,nGordon] = size(matrix_Gordon);
load([resultsDir '/dhatsAllSubjectsR_2']);
nSubject = size(dhatsAllSubjectsParcel,1);


ntask=size(contrastmat,2);
ncontrast = size(contrastmat,1);
resultsBetas = zeros(nVertex,ntask);
resultsContrasts = zeros(nVertex,ncontrast);
resultsShats = zeros(nGordon,ntask);
resultsBhats = zeros(nGordon,ntask);
resultsThetas = zeros(nGordon,ntask);
resultsTstatBetas = resultsBetas;
resultsTstatContrasts = resultsContrasts;
subjResultsContrast = zeros(nSubject,size(contrastmat,1),nVertex);
resultsAllThetas = zeros(nGordon,ntask,3);

ntask=size(contrastmat,2);
ncontrast = size(contrastmat,1);
resultsBetas = zeros(nVertex,ntask);
resultsContrasts = zeros(nVertex,ncontrast);
resultsShats = zeros(nGordon,ntask);
resultsBhats = zeros(nGordon,ntask);
resultsThetas = zeros(nGordon,ntask);
resultsTstatBetas = resultsBetas;
resultsTstatContrasts = resultsContrasts;
subjResultsContrast = zeros(nSubject,size(contrastmat,1),nVertex);
eBLUPsSubject = zeros(nGordon,nSubject,ntask);
eBLUPsSubjectVertex = zeros(nVertex,nSubject,ntask);

for iParcel=2:nGordon    
   indices = find(matrix_Gordon(:,iParcel));
   latParcel = lat(indices);
   longParcel = long(indices);
   nVp = length(indices);
   parcelDist = zeros(nVp,nVp);
   for v=1:nVp
        %parcelDist(v,:) = mygreatcirc(latParcel(v),longParcel(v),latParcel,longParcel);
        parcelDist(v,:) = distance(latParcel(v),longParcel(v),latParcel,longParcel);
   end
    load([resultsDir '/dhatsAllSubjectsR_',num2str(iParcel)]);
    
    saveTemp = [resultsDir '/Covariograms/CovariogramsGordonNetwork',num2str(iParcel)];
    
    est = estSTMM1(dhatsAllSubjectsParcel,varQAllSubjectsParcel,parcelDist,iParcel,contrastmat,resultsDir,saveTemp);
    
    disp(['Gordon Parcel: ', num2str(iParcel)]);
   
    % write to master results matrix:
    resultsBetas(indices,:) = est.betaHatVertex;
    resultsContrasts(indices,:) = est.contrast_betaHatVertex;
    resultsTstatBetas(indices,:) = est.t_betaHatVertex;
    resultsTstatContrasts(indices,:) = est.t_contrast_betaHatVertex;    
    subjResultsContrast(:,:,indices) = est.contrast_eBLUPsBetaHat;
    resultsThetas(iParcel,:) = est.bTheta(:,2);
    resultsShats(iParcel,:) = est.Shat;
    resultsBhats(iParcel,:) = est.Bhat;
    resultsAllThetas(iParcel,:,:) = est.bTheta;
    eBLUPsSubjectVertex(indices,:,:)=permute(est.eBLUPsSubjectVertex,[3,1,2]);
    eBLUPsSubject(iParcel,:,:) = est.eBLUPsSubject;    
end
 
save([resultsDir '/STMMAllSubjects'],...
    'resultsBetas','resultsContrasts','resultsTstatBetas','resultsTstatContrasts','subjResultsContrast',...
    'resultsThetas','resultsShats','resultsBhats','resultsAllThetas','eBLUPsSubjectVertex','eBLUPsSubject');
end
