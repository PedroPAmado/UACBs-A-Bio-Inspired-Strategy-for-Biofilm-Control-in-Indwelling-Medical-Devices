%% image_analysis.m

%Authors: ARTORG-Center, CVE-Group;
% - Damiano Schirinzi
% - Aurelia Bucciarelli

%% tabula rasa

clearvars -except inputVideoPath inputVideoNames

%Load Inputs from GUI into workspace and delete file containing them 
load('.\Temporary\UserSettingsGUItemp.mat')
UserSettingsGUI = {UserSettingsGUIstring,UserSettingsGUIvalue};

%% Haematocrit calculation
if(UserSettingsGUIvalue.haematocritCheckBox)
    if(UserSettingsGUIvalue.nDomainsEdit == 0)
        disp('Skipped haematocrit calculation; nr. of domains = 0')
    else
        haematocrit(inputVideoPath,inputVideoNames,UserSettingsGUI);
    end
end

%% PTV Calculation
if(UserSettingsGUIvalue.PTVCheckBox)
    %Retrieve PTV settings
    
    PTVsettings = struct('ParticleDetectionAlgorithm',UserSettingsGUIvalue.ParticleDetectionAlgorithmDropDown,...
        'corrthreval',UserSettingsGUIvalue.corrThresholdEdit,...
        'sigmasize',UserSettingsGUIvalue.sigmaEdit,...
        'intthreval',UserSettingsGUIvalue.intensityThresholdEdit,...
        'PTValgorithm',UserSettingsGUIvalue.PTValgorithmDropDown,...
        'Crosscorrlparameters', UserSettingsGUIvalue.MethodDropDown,...
        'minParticles',UserSettingsGUIvalue.nrParticlesEdit,...
        'interrArea',UserSettingsGUIvalue.interrogationAreaEdit,...
        'minCorr',UserSettingsGUIvalue.minCorrelationEdit,...
        'simNeigh',UserSettingsGUIvalue.similarityNeighEdit,...
        'visualPTVControl',UserSettingsGUIvalue.visualPTVControlCheckBox);
    
    %Call PTV Algorithm
    for i=1:length(inputVideoNames)
        %Load Frames
        Frames = load([inputVideoPath inputVideoNames{i}(1:(end-4)) '.mat' ]);
        %Load Haematocrit analysis results with particle diameter and pass it
        %to the PTV settings if stated so in the GUI
        if(UserSettingsGUIvalue.setSigmaAutoCheckBox)
            load([inputVideoPath inputVideoNames{i}(1:(end-4)) '.mat' ]);
            UserSettingsGUIvalue.sigmaEdit = num2str(HtcResults.Diameter);
            PTVsettings.sigmasize = round(median(HtcResults.Diameter));
        end
        disp(['PTV processing: ' inputVideoNames{i}(1:(end-4))])
        disp(['sigmasize: ' num2str(PTVsettings.sigmasize)])
        %Call PTV Algorithm
        particleInfo = PTVlab_artorg(Frames.FrameMat, Frames.UserSettings, Frames.ProcessingSteps, PTVsettings);
        %Add 7th column to particleInfo with the magnitude of velocity vector:
        %v = sqrt(v_x^2+v_y^2), leave 6th column of particleInfo empty
        particleInfo(:,7) = sqrt(particleInfo(:,4).^2+particleInfo(:,5).^2);
        %Save particleInfo and PTVsettings
        save([ inputVideoPath 'PTV_' inputVideoNames{i}(1:(end-4)) '.mat'],'particleInfo','PTVsettings')
        disp(['saved PTV results to: PTV_' inputVideoNames{i}(1:(end-4)) '.mat'])
        disp('-------------------------------------------------------')
    end
end

%% PTV vs Haematocrit
if(UserSettingsGUIvalue.haematocritVsPTVCheckBox)
    PTVvsHaematocrit(inputVideoPath, inputVideoNames, UserSettingsGUI)
end
%% average RBC velocity / RBC line density / RBC flux
if(UserSettingsGUIvalue.generalAnalysisCheckBox)
   generalAnalysis(inputVideoPath, inputVideoNames, UserSettingsGUI)
end
%CountD


