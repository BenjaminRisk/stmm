function [] = smoothARparameters(myResultsDir,mySubjectDir,subjects,hemisphere)

    load([myResultsDir,'/UnsmoothedARparameters.mat']);
    
    if strcmp(hemisphere,'CORTEX_RIGHT')
            load('./supportingdatafiles/matrix_GordonRev2_CII_CORTEX_RIGHT');
            matrix_Gordon = matrix_GordonRev2_CII_CORTEX_RIGHT;
    elseif strcmp(hemisphere,'CORTEX_LEFT')
            load('./supportingdatafiles/matrix_GordonRev2_CII_CORTEX_LEFT');
            matrix_Gordon = matrix_GordonRev2_CII_CORTEX_LEFT;
    else error('hemisphere = CORTEX_RIGHT or CORTEX_LEFT');
    end
        
    smVarResAllSubjects = zeros(nSubject,nVertex);
    smPhiAllSubjects = zeros(nSubject,3,nVertex);
    smVarQAllSubjects = zeros(nSubject,nCov,nCov,nVertex);
    smVarInnovAllSubjects = zeros(nSubject,nVertex);

    % create array of all design matrices so that parfor can be used:
    load([mySubjectDir,'/',char(subjects(1)),'/subjDesignMat.mat']);
    allDesignMat = zeros(nSubject,size(subjDesignMat,1),size(subjDesignMat,2));
    for n=1:nSubject
        load([mySubjectDir,'/',char(subjects(n)),'/subjDesignMat.mat']);
        allDesignMat(n,:,:) = subjDesignMat;
    end
    %NOTE: the implementation using parfor requires a lot of memory so if this is being run with many
    %subjects, it is best to use 'for'


    distMat = zeros(nVertex,nVertex);
    parfor v=1:nVertex
        distMat(v,:) = distance(lat(v),long(v),lat,long);
        %distMatR(v,:) = mygreatcirc(latR(v),longR(v),latR,longR);
    end

    %parfor v=1:nVertex 
    for v=1:nVertex
        dist = distMat(:,v);
        index = (dist<=h_ob);
        wts = biweight(dist(index)/h_ob)/h_ob;
        wts = wts/sum(wts);

        phiNeighbors = phiAllSubjects(:,:,index);
        rhoNeighbors = rhoAllSubjects(:,:,index);
        sigmasqTall = varResAllSubjects(:,index)*wts;

        phiTall1= squeeze(phiNeighbors(:,1,:))*wts;
        phiTall2= squeeze(phiNeighbors(:,2,:))*wts;
        phiTall3= squeeze(phiNeighbors(:,3,:))*wts;
        phiTall = [phiTall1, phiTall2, phiTall3];
        rhoTall1 = squeeze(rhoNeighbors(:,1,:))*wts;
        rhoTall2 = squeeze(rhoNeighbors(:,2,:))*wts;
        rhoTall3 = squeeze(rhoNeighbors(:,3,:))*wts;
        rhoTall = [rhoTall1, rhoTall2, rhoTall3];

        smVarResAllSubjects(:,v) = sigmasqTall;
        smPhiAllSubjects(:,:,v) = phiTall;
        smVarInnovAllSubjects(:,v) = varInnovAllSubjects(:,index)*wts;
        
        
        p=3;
        for n = 1:nSubject
            subjDesignMat = squeeze(allDesignMat(n,:,:));
            cap = (subjDesignMat'*subjDesignMat)\subjDesignMat';   
            corrTemp = zeros(nSession,1);
            corrTemp(1) = 1;
            corrTemp(2:(p+1)) = rhoTall(n,:);  
            for t=(p+2):nSession
                 corrTemp(t) = phiTall(n,:)*corrTemp((t-1):-1:(t-p));
            end
            covTemp = sigmasqTall(n)*corrTemp;
            covmat = toeplitz(covTemp);
            covmat = kron(eye(2),covmat);
            varQ = cap*covmat*cap';
            % this results in a matrix that is not quite symmetric; symmetrize:
            varQ = (varQ + varQ')/2; 
            smVarQAllSubjects(n,:,:,v) = varQ;
        end
    end
    % 42 minutes without parfor for motor

    save([myResultsDir '/SmoothedARparameters.mat'],'dhatsOLSAllSubjects','smVarInnovAllSubjects',...
        'smPhiAllSubjects','smVarResAllSubjects');

    % save smoothed varQ with dhats for each parcel
    % include covariates that will be included in the second level analysis for each revised2 Gordon Parcel
    nGordon = size(matrix_Gordon,2);
    for m=2:nGordon
        indices = find(matrix_Gordon(:,m));
        dhatsAllSubjectsParcel = dhatsOLSAllSubjects(:,logical(contrastMat),indices);
        varQAllSubjectsParcel = smVarQAllSubjects(:,logical(contrastMat),logical(contrastMat),indices);
        smVarResAllSubjectsParcel = smVarResAllSubjects(:,indices);
        smVarInnovAllSubjectsParcel = smVarInnovAllSubjects(:,indices);
        smPhiAllSubjectsParcel = smPhiAllSubjects(:,:,indices);
        save([myResultsDir '/dhatsAllSubjectsR_',num2str(m)],...
            'dhatsAllSubjectsParcel','varQAllSubjectsParcel','smVarResAllSubjectsParcel',...
            'smVarInnovAllSubjectsParcel','smPhiAllSubjectsParcel');
    end
end