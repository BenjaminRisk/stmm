function [] = smoothARparameters(optbw,resultsDir,subjectDir,subjects,contrastmat,lat,long,hemisphere)

    load([resultsDir,'/UnsmoothedARparameters.mat']);
    [nSubject,nVertex] = size(varResAllSubjects);
    if strcmp(hemisphere,'CORTEX_RIGHT')
            load('./supportingdatafiles/matrix_GordonRev2_CII_CORTEX_RIGHT');
            matrix_Gordon = matrix_GordonRev2_CII_CORTEX_RIGHT;
    elseif strcmp(hemisphere,'CORTEX_LEFT')
            load('./supportingdatafiles/matrix_GordonRev2_CII_CORTEX_LEFT');
            matrix_Gordon = matrix_GordonRev2_CII_CORTEX_LEFT;
    else error('hemisphere = CORTEX_RIGHT or CORTEX_LEFT');
    end
        
    % create array of all design matrices so that parfor can be used:
    load([subjectDir,'/',char(subjects(1)),'/subjDesignMat.mat']);
    [nTime,nCov] = size(subjDesignMat);
    nSession = nTime/2;
    allDesignMat = zeros(nSubject,nTime,nCov);
    for n=1:nSubject
        load([subjectDir,'/',char(subjects(n)),'/subjDesignMat.mat']);
        allDesignMat(n,:,:) = subjDesignMat;
    end
    %NOTE: the implementation using parfor requires a lot of memory so if this is being run with many
    %subjects, it is best to use 'for'

    distMat = zeros(nVertex,nVertex);
    parfor v=1:nVertex
        distMat(v,:) = distance(lat(v),long(v),lat,long);
        %distMatR(v,:) = mygreatcirc(latR(v),longR(v),latR,longR);
    end

    wts = biweight(distMat./optbw)/optbw;
    clear distMat;

    smwts = sum(wts)';
    wts2 = zeros(nVertex,nVertex);
    for v=1:nVertex
        wts2(:,v) = wts(:,v)./smwts;
    end
    clear wts;

    smPhiAllSubjects = zeros(nSubject,3,nVertex);
    smRhoAllSubjects = zeros(nSubject,3,nVertex);
 

    smPhiAllSubjects(:,1,:) = squeeze(phiAllSubjects(:,1,:))*wts2;
    smPhiAllSubjects(:,2,:) = squeeze(phiAllSubjects(:,2,:))*wts2;       
    smPhiAllSubjects(:,3,:) = squeeze(phiAllSubjects(:,3,:))*wts2;
        
    smRhoAllSubjects(:,1,:) = squeeze(rhoAllSubjects(:,1,:))*wts2;
    smRhoAllSubjects(:,2,:) = squeeze(rhoAllSubjects(:,2,:))*wts2;
    smRhoAllSubjects(:,3,:) = squeeze(rhoAllSubjects(:,3,:))*wts2;
    
    smVarResAllSubjects = varResAllSubjects*wts2;
    smVarInnovAllSubjects = varInnovAllSubjects*wts2;
       
    clear wts2;
    smVarQAllSubjects = zeros(nSubject,nCov,nCov,nVertex);

        p=3;
        for v=1:nVertex
            for n = 1:nSubject
                subjDesignMat = squeeze(allDesignMat(n,:,:));
                cap = (subjDesignMat'*subjDesignMat)\subjDesignMat';   
                corrTemp = zeros(nSession,1);
                corrTemp(1) = 1;
                corrTemp(2:(p+1)) = smRhoAllSubjects(n,:,v);  
                for t=(p+2):nSession
                     corrTemp(t) = smPhiAllSubjects(n,:,v)*corrTemp((t-1):-1:(t-p));
                end
                covTemp = smVarResAllSubjects(n,v)*corrTemp;
                covmat = toeplitz(covTemp);
                covmat = kron(eye(2),covmat);
                varQ = cap*covmat*cap';
                % this results in a matrix that is not quite symmetric; symmetrize:
                varQ = (varQ + varQ')/2; 
                smVarQAllSubjects(n,:,:,v) = varQ;
            end
       end


    save([resultsDir '/SmoothedARparameters.mat'],'dhatsOLSAllSubjects','smVarInnovAllSubjects',...
        'smPhiAllSubjects','smVarResAllSubjects');

    % save smoothed varQ with dhats for each parcel
    % include covariates that will be included in the second level analysis for each revised2 Gordon Parcel
    nGordon = size(matrix_Gordon,2);
    for m=2:nGordon
        indices = find(matrix_Gordon(:,m));
        dhatsAllSubjectsParcel = dhatsOLSAllSubjects(:,logical(contrastmat),indices);
        varQAllSubjectsParcel = smVarQAllSubjects(:,logical(contrastmat),logical(contrastmat),indices);
        smVarResAllSubjectsParcel = smVarResAllSubjects(:,indices);
        smVarInnovAllSubjectsParcel = smVarInnovAllSubjects(:,indices);
        smPhiAllSubjectsParcel = smPhiAllSubjects(:,:,indices);
        save([resultsDir '/dhatsAllSubjectsR_',num2str(m)],...
            'dhatsAllSubjectsParcel','varQAllSubjectsParcel','smVarResAllSubjectsParcel',...
            'smVarInnovAllSubjectsParcel','smPhiAllSubjectsParcel');
    end
end