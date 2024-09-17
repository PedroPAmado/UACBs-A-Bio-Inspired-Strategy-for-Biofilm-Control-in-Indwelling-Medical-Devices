function [] = PTVvsHaematocrit(inputVideoPath, inputVideoNames, UserSettingsGUI)

UserSettingsGUIvalue = UserSettingsGUI{2};
domainROI = UserSettingsGUIvalue.domainROI;

for k=1:length(inputVideoNames) %Loop over all analyzed videos
    disp(['PTV vs Haematocrit of: ' inputVideoNames{k}])
    %Load HtcResults:
    try
        load([inputVideoPath 'Htc_' inputVideoNames{k}(1:(end-4)) '.mat'])
    catch
        disp(['Error: Could not load "' inputVideoPath 'Htc_' inputVideoNames{k}(1:(end-4)) '.mat"'])
        return
    end
    Ndomains = size(HtcResults.Area,1);
    %Load particleInfo:
    try
        load([inputVideoPath 'PTV_' inputVideoNames{k}(1:(end-4)) '.mat'])
    catch
        disp(['Error: Could not load "' inputVideoPath 'PTV_' inputVideoNames{k}(1:(end-4)) '.mat"'])
        return
    end
    
    PTVsettings=PTVsettings;
    
    %Allocate memory for output of different Ndomains Outputs
    error_mat = zeros(Ndomains,size(HtcResults.Ncells,2));
    CRID_vec = zeros(Ndomains,1);
    for j=1:Ndomains
        disp(['     Domain: ' num2str(j)])
        nFrames = length(HtcResults.Ncells(j,:)); %Number of frames
        %Allocate vector containing n. of particles per frame detected by PTV
        PTVparticles = zeros(1,nFrames);
        
        %Retain in particleInfoCut only partilces that are in ROI
        in = zeros(length(particleInfo),1);
        %Convert domainROI to polygon edges
        XY  = bbox2points(domainROI(j,:));
        for i1=1:length(particleInfo)
            in(i1) = inpolygon(particleInfo(i1,2),particleInfo(i1,3),XY(:,1),XY(:,2));
        end
        in = find(in==1); %det indexes of the row of particles that are inside the ROI
        particleInfoCut = particleInfo(in,1:7);
        
        for i=1:nFrames
            PTVrows = find(particleInfoCut(:,1)==i);
            if(isempty(PTVrows)) %No particles found in this frame
                PTVparticles(i) = 0;
            else %Particles found
                PTVparticles(i) = size(PTVrows,1);
            end
        end
        
        %Particle tracking tracks only until second last frame.
        %Here we assume, that the last frame contains as many particles as
        %the second last one to minimize error to the Haematocrit script.
        PTVparticles(end) = PTVparticles(end-1);
        
        %Detected particles PTV vs Haematocrit [%]
        %This calculation gives wrong results if there are no particles
        detectedParticles = PTVparticles./round(HtcResults.Ncells(j,:))*100;
        %coefficient of variance: CV = std/mean
        
        %Coefficient of relative integral difference
        framestep = 1;
        diffIntegral = (PTVparticles-round(HtcResults.Ncells(j,:)))*framestep;
        absDiffIntegral = sum(abs(diffIntegral));
        intergralHaematocrit = sum(round(HtcResults.Ncells(j,:))*framestep);
        CRID = absDiffIntegral/intergralHaematocrit;
        disp(['       CRID = ' num2str(CRID)])
        %Root mean squared difference
        diffIntegral = sum(sqrt((PTVparticles-round(HtcResults.Ncells(j,:))).^2));
        
        %Mean Absolute Error MAE
        error = PTVparticles-round(HtcResults.Ncells(j,:));
        error_per=abs(error./round(HtcResults.Ncells(j,:)));
        disp('      error statistics:')
        fprintf('       median error:  %d \n',median(error))
        fprintf('       mean error:  %3.2f \n',median(error))
        fprintf('       std of error:  %3.2f \n',std(error))
        fprintf('       max error:     %d \n',max(error))
        fprintf('       min error:     %d \n',min(error))
        MAE = sum(abs(PTVparticles-round(HtcResults.Ncells(j,:))))/nFrames;
              
        %Plots
        figure
        hold on
        plot(1:nFrames,round(HtcResults.Ncells(j,:)),'b')
        plot(1:nFrames,mean(round(HtcResults.Ncells(j,:)))*ones(1,nFrames),'b--')
        plot(1:nFrames,PTVparticles,'r')
        plot(1:nFrames,mean(PTVparticles)*ones(1,nFrames),'r--')
        legend('Htc Cal.','mean Htc Cal.','PTV','mean PTV')
        title({['Video: ' inputVideoNames{k} '; Domain ' num2str(j)], 'Number of cells per frame'})
        xlabel('frame')
        ylabel('cells')
        savefig([inputVideoPath 'PTVvsHtc_D' num2str(j) '_' inputVideoNames{k}(1:(end-4)) '.fig'])
        disp(['     saved Htc figure to: ' 'PTVvsHtc_D' num2str(j) '_' inputVideoNames{k}(1:(end-4)) '.fig'])
        
        %Save Output variables
        error_mat(j,:) = error;
        CRID_vec(j,1) = CRID;
        
    end
    %Add error and CRID measure to the PTV results file:
    Ncells=HtcResults.Ncells;    
    error = error_mat;
    CRID = CRID_vec;
    
    save([inputVideoPath 'PTV_' inputVideoNames{k}(1:(end-4)) '.mat'],'particleInfo','Ncells','PTVsettings','error','error_per','CRID')
    disp(['Added error and CRID to PTV results in: ' 'PTV_' inputVideoNames{k}(1:(end-4)) '.mat'])
    disp('-------------------------------------------------------')
       
end
end