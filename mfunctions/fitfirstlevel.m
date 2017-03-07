function [] = fitfirstlevel(subjects,subjectDir,resultsDir)

    load([subjectDir,'/',char(subjects(1)),'/matData_RIGHT_CORTEX.mat']);
    load([subjectDir,'/',char(subjects(1)),'/subjDesignMat.mat']);

    nSubject = length(subjects);
    nCov = size(subjDesignMat,2);
    nVertex = size(fData,1);
    dhatsOLSAllSubjects = zeros(nSubject,nCov,nVertex);
    dhatsGLSAllSubjects = dhatsOLSAllSubjects;
    varResAllSubjects = zeros(nSubject,nVertex);
    varInnovAllSubjects = zeros(nSubject,nVertex);
    phiAllSubjects = zeros(nSubject,3,nVertex);
    rhoAllSubjects = phiAllSubjects;

    for n = 1:nSubject
        subject = char(subjects(n)); 
        if n>1 
            load([subjectDir,'/',subject,'/matData_RIGHT_CORTEX.mat']);
            load([subjectDir,'/',subject,'/subjDesignMat.mat']);
        end
        cap = (subjDesignMat'*subjDesignMat)\subjDesignMat';
        hatter = subjDesignMat*cap;
        myInvMsubj = inv(createM_reducedbias_2sess(hatter,20));
        for v = 1:nVertex
            fdatatemp = fData(v,:)';
            [dhatsOLSAllSubjects(n,:,v),phiTemp,sigmasqTemp,covTemp,rhoTemp] = fitarpOLSYW(fdatatemp,subjDesignMat,3,myInvMsubj);
            invTemp = covTemp\subjDesignMat;
            intTemp = subjDesignMat'*invTemp;
            invInvTemp = intTemp\invTemp';
            dhatsGLSAllSubjects(n,:,v) = invInvTemp*squeeze(fdatatemp);
            varResAllSubjects(n,v) = covTemp(1,1);
            varInnovAllSubjects(n,v) = sigmasqTemp;
            phiAllSubjects(n,:,v) = phiTemp;
            rhoAllSubjects(n,:,v) = rhoTemp;
        end
    end
    save([resultsDir,'/UnsmoothedARparameters.mat'],...
        'dhatsOLSAllSubjects','dhatsGLSAllSubjects','varInnovAllSubjects','phiAllSubjects','rhoAllSubjects','varResAllSubjects');
end
