function [invMAllSubjects] = createInvM(subjectDesignMatBaseDir,subjects)
    nSubject = length(subjects);
    invMAllSubjects = zeros(nSubject,21,21);
    tic;
    for n=1:nSubject
        load([subjectDesignMatBaseDir,'/',char(subjects(n)),'/subjDesignMat.mat']);
        hatter = subjDesignMat/(subjDesignMat'*subjDesignMat)*subjDesignMat';
        invMAllSubjects(n,:,:) = inv(createM_reducedbias_2sess(hatter,20));
    end
    toc
    %save([invMsaveDir,'/myInvMAllSubjects_L20.mat'],'myInvMAllSubjects');   
end
