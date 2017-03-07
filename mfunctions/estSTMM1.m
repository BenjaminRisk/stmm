function [out] = estSTMM1(dhat,vardhat,distmat,iParcel,contrastmat,resultsDir,savecovario)
% estSTMM1: Estimate the second-level of the spatiotemporal model with
%   fixed vertex effects.
% Input:
%   dhat: N x Q x V
%   vardhat: N x Q x Q x V covariance matrix of K_n Psi_{nv} K_n'
%   distmat: V x V
%   iParcel: parcel Index
%   contrastmat: Number of contrasts X Q vector of linear combinations used in
%   defining contrasts
%   savecovario: path to save figure of covariograms
% Output:
%   out.Shat
%   out.Bhat
%   out.Btheta
%   out.eBLUPsSubject
%   out.eBLUPsSubjectVertex - BLUPs for the subject-vertex interaction RE
%   out.eBLUPsBetaHat - total predicted activation for a subject
%   out.betaHatVertex
%   out.t_betaHatVertex - vertex-wise t-statistics
%   out.contrast_betaHatVertex 
%   out.t_contrast_betaHatVertex
%   out.contrast_eBLUPsBetaHat        
        
    [nSubject,Q,nParcelm] = size(dhat);
    mindist = min(distmat(distmat~=0));
    hbin = linspace((mindist-eps),max(distmat(:))/2,21);
        
    % Refine hbin to only include bins with sufficient data: use data from arbitrary subject and
    % arbitrary covariate -- the number of observations per bin will be the same
    [~,~,cnt,avedist] = ecovariogram(squeeze(dhat(1,1,:))',distmat,hbin,false);   
    check = cnt<10;
    % this approach eliminates zero bins or bins with fewer than 10 cnts, but keeps same distances
    if sum(check)
        hbin = hbin(logical([~check;1]));
        [~,~,~,avedist] = ecovariogram(squeeze(dhat(1,1,:))',distmat,hbin,false);   
    end
    csvwrite([resultsDir '/Covariograms/avedist_Parcel' num2str(iParcel) '.csv'],avedist);

    
    meanKnPsiKn = mean(squeeze(mean(vardhat,1)),3);    
    bThetaMat = zeros(Q,3);
    Omegahat = zeros(Q,nParcelm,nParcelm);
    Shat = zeros(Q,1);
    Bhat = zeros(Q,1);
    %% BEGIN LOOP OVER q
   % Estimate Omega:
 for q=1:Q
    dhatTemp = squeeze(dhat(:,q,:))';
    dbarv = mean(dhatTemp,2);
    ecovarioAll = zeros(length(avedist),nSubject);
    mseTask = mseUSB(dhatTemp');
    parfor n=1:nSubject
        dhatsubj = dhatTemp(:,n);        
        ecovarioAll(:,n) = (nSubject/(nSubject-1))*ecovariogramSTMM1(dhatsubj,dbarv,distmat,hbin);
    end
    % Estimate the covariogram for xMental:
    csvwrite([resultsDir '/Covariograms/ecovarioAll_task',num2str(q),'_Parcel',num2str(iParcel),'.csv'],ecovarioAll);

    meanEcovario = mean(ecovarioAll,2);
    
    % Calculate reasonable optimization space for fmincon (internal of fitexponcov):
    crudeS = mseTask.S/nParcelm;
    crudeB = mseTask.B;
    lbb = [0,1e-03,-crudeS];
    ubb = [2*crudeB,10,2*crudeS];
    
    % use two initializations, one from nearly zero dependence and the other from high dependence,
    % then choose the best:
    initb = [crudeB,9,crudeS];
    [bTheta,bFunc,fval1] = fitexponcov(avedist,meanEcovario,lbb,ubb,initb);
    initb = [crudeB,0.05,crudeS];
    [bTheta2,bFunc2,fval2] = fitexponcov(avedist,meanEcovario,lbb,ubb,initb);
    if fval2<fval1
        bTheta = bTheta2;
        bFunc = bFunc2;
    end
    
   
    % create object of thetas for output:
    bThetaMat(q,:) = bTheta;
    
    %% Estimate variance B:
    Omegahat(q,:,:) = exp(-bTheta(2)*distmat);
    omegasum = sum(sum(squeeze(Omegahat(q,:,:))));
    Bhat(q) = (mseTask.B - meanKnPsiKn(q,q))/(nParcelm/(nParcelm-1) - omegasum/(nParcelm^2 - nParcelm));
    
    if sum(Bhat(q)<0)
        disp(['MOM estimator of Bhat(' num2str(q) ') for parcel ' num2str(iParcel) ' is negative -- replacing with 1e-06']);
        Bhat(Bhat<0) = 1e-6;
    end

    %% Estimate variance S:
    Shat(q) = mseTask.S/nParcelm - omegasum*Bhat(1)/nParcelm^2 - meanKnPsiKn(q,q)/nParcelm;
    if Shat(q)<0
        disp(['MOM estimator of Shat(' num2str(q) ') for parcel ' num2str(iParcel) ' is negative -- replacing with 1e-06']);
        Shat(q) = 1e-6;
    end
    
    %% Create figures to assess covariogram fit
        avedistplus = kron(ones(nSubject,1),avedist);
        a=subplot(1,2,1);
        plot(avedist,mean(ecovarioAll,2),'k','Linewidth',3)
        xlim([0,max(avedist)*1.1]);
        hold on;
        plot(avedist,bFunc(avedist),'r','Linewidth',3)
        text(0.5,min(ecovarioAll(:)),['\theta = ' num2str(bTheta(2))])
        title(['Gordon Network ',num2str(iParcel),' q=',num2str(q),...
            ' Subject*Vertex'],'FontSize',12)
        gscatter(avedistplus,ecovarioAll(:),kron([1:nSubject]',ones(length(avedist),1)))
        legend('Mean','Fitted covariogram','Emp covario per subject');
        
        if Bhat(q)==1e-06
           text(0.5,(min(ecovarioAll(:))+max(ecovarioAll(:)))/2,'NOTE: \sigma_{bq}^2 = 0');
        end
        ylabel('\delta(h)')
        xlabel('Distance (h)')
        hold off;
    
        subplot(1,2,2);
        plot(avedist,meanEcovario,'k','Linewidth',3)
        xlim([0,max(avedist)*1.1]);
        hold on;
        plot(avedist,bFunc(avedist),'r','Linewidth',3)
        text(0.5,min(meanEcovario)+0.05*abs(min(meanEcovario)),['\theta = ' num2str(bTheta(2))])
        ylabel('\delta(h)')
        xlabel('Distance (h)')
        hold off;
    
    if nargin==7
            saveas(a,[savecovario '_q=',num2str(q),'.pdf']);
    end
 end  
    %% END LOOP OVER q
    
    
    % memory clean-up:
    clear('distmat');
     
    %% Create covariance and precision matrices:
    % NOTE: we are reordering the indices such that the dimensions 
    % are N x Q x V:
    %subject-specific components:
    subCov = zeros(nSubject,nParcelm*Q,nParcelm*Q);

    for n=1:nSubject
        for v=1:nParcelm
           for q=1:Q
                iIndex = nParcelm*(q-1)+v;
                for qprime=1:Q
                    jIndex = nParcelm*(qprime-1)+v;
                    subCov(n,iIndex,jIndex) = vardhat(n,q,qprime,v);
                end
           end
        end
         for q=1:Q
               starting = nParcelm*(q-1)+1;
               ending  = starting + nParcelm-1;
               subCov(n,starting:ending,starting:ending) = subCov(n,starting:ending,starting:ending) + ...
                   Omegahat(q,:,:).*Bhat(q)+Shat(q);             
         end
    end
    clear('vardhat');
    
    %create inverse of subject-independent effects:
    % with 182 vertices, this takes 2.54 seconds with sparse matrices versus 0.5 seconds with full. 
    subPrec = zeros(nSubject,nParcelm*Q,nParcelm*Q);
    
    for n=1:nSubject
        subPrec(n,:,:) = inv(squeeze(subCov(n,:,:)));
    end
    clear('subCov');
    
    subPrecWide = zeros(nParcelm*Q,nParcelm*Q*nSubject);
    for n=1:nSubject
        starting = (n-1)*(nParcelm*Q)+1;
        subPrecWide(:,starting:(nParcelm*Q*n))=squeeze(subPrec(n,:,:));
    end

    % create covariance matrix for the fixed vertex effects 
    % in the no-intercept parameterization:
    precBetaHatVertex = squeeze(sum(subPrec,1));
    covBetaHatVertex = inv(precBetaHatVertex);
        
  
    dhat = permute(dhat, [3 2 1]);
    dhat = dhat(:);
        
    betaHatVertex = covBetaHatVertex*subPrecWide*dhat;
    %betaHatVertex = precBetaHatVertex\subPrecWide*dhat;
    resids = dhat - kron(ones(nSubject,1),betaHatVertex);
    
    
    %% Calculate random effects:
    eBLUPsSubjectVertex = zeros(nSubject,Q,nParcelm);
    eBLUPsSubject = zeros(nSubject,Q);
    eBLUPsBetaHat = eBLUPsSubjectVertex;
    for n=1:nSubject
        starting = 1+(n-1)*Q*nParcelm;
        ending = n*Q*nParcelm;
        blupInt = squeeze(subPrec(n,:,:))*resids(starting:ending);
        for q=1:Q
            starting=1+(q-1)*nParcelm;
            ending = q*nParcelm;
            tempBLUP = blupInt(starting:ending);
            eBLUPsSubjectVertex(n,q,:) = Bhat(q)*squeeze(Omegahat(q,:,:))*tempBLUP;
            eBLUPsSubject(n,q) = Shat(q)*sum(tempBLUP);
            eBLUPsBetaHat(n,q,:) = squeeze(eBLUPsSubjectVertex(n,q,:))+...
                betaHatVertex(starting:ending)+eBLUPsSubject(n,q);
        end
    end
    
    
        tStatVertex = betaHatVertex./sqrt(diag(covBetaHatVertex));        
        out.betaHatVertex=reshape(betaHatVertex,[nParcelm,Q]);
        out.t_betaHatVertex = reshape(tStatVertex,[nParcelm,Q]);
        out.contrast_betaHatVertex = out.betaHatVertex*contrastmat';
        
        nContrast = size(contrastmat,1);
        out.t_contrast_betaHatVertex = zeros(nParcelm,nContrast);
        for v=1:nParcelm
            indices = v:nParcelm:(nParcelm*Q);
            for r=1:nContrast
                varTemp = contrastmat(r,:)*covBetaHatVertex(indices,indices)*contrastmat(r,:)';
                out.t_contrast_betaHatVertex(v,r) = out.contrast_betaHatVertex(v,r)/sqrt(varTemp);
            end
        end
        
        
        out.eBLUPsSubjectVertex=eBLUPsSubjectVertex;
        out.eBLUPsSubject = eBLUPsSubject;
        out.eBLUPsBetaHat = eBLUPsBetaHat;
        
        out.contrast_eBLUPsBetaHat = zeros(nSubject,nContrast,nParcelm);
        for n=1:nSubject
            out.contrast_eBLUPsBetaHat(n,:,:) = contrastmat*squeeze(eBLUPsBetaHat(n,:,:));
        end
        
        out.Shat = Shat;
        out.Bhat = Bhat;
        out.bTheta = bThetaMat;        
end
