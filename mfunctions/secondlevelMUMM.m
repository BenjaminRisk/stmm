function []=secondlevelMUMM(contrastMat,dhatsAllSubjects,resultsDir)
    nContrast = size(contrastMat,1);
    [nSubject,~,nVertex] = size(dhatsAllSubjects);
    
    subjResultsContrast = zeros(nSubject,nContrast,nVertex);

    for n=1:nSubject
        subjResultsContrast(n,:,:) = contrastMat*squeeze(dhatsAllSubjects(n,:,:));
    end

    resultsBetasTMUMM = squeeze(mean(dhatsAllSubjects,1))';
    varPop = squeeze(var(dhatsAllSubjects))'./nSubject;
    resultsTstatBetasTMUMM = resultsBetasTMUMM./sqrt(varPop);

    resultsContrastsTMUMM = squeeze(mean(subjResultsContrast,1));
    varTemp = squeeze(var(subjResultsContrast)./nSubject);
    resultsTstatContrastsTMUMM = resultsContrastsTMUMM./sqrt(varTemp);

    save([resultsDir '/TMUMMAllSubjects'],...
        'subjResultsContrast','resultsBetasTMUMM','resultsTstatBetasTMUMM','resultsContrastsTMUMM','resultsTstatContrastsTMUMM');
end