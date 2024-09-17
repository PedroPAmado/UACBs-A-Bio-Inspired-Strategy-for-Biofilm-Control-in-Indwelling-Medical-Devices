%% tabula rasa
clc
clear all
close all

%Clear temporary directory
addpath('H:\Cornel_ First Experiments\data_analysis\scripts_24.05.23\PTV\Image_analysis_20221211_PEdro\Temporary')
delete('Temporary\*')

%Add Library with all subfunctions to Matlab working path
allpaths = genpath('./Library'); %get path to dir Library and subfolders
addpath(allpaths)

%% GUI
%Dialogue for choosing a directory:
inputVideoPath = uigetdir('','Select the folder containing the videos to be processed:');
if(inputVideoPath(end)~='\')
    inputVideoPath = [inputVideoPath '\'];
end
%Get the names of all video files in the chosen folder:
inputVideoNames = dir([inputVideoPath '*.tif']); 
inputVideoNames = {inputVideoNames.name};
%Start user input gui
warning('off','all')
guiHandle = CVE_GUI;
waitfor(guiHandle)
%%
try
    load('.\Temporary\UserSettingsGUItemp.mat')   
    if(UserSettingsGUIvalue.selectNewROIs)
        %ROI slection GUI
        [angle_rotation,ROI,rect,domainROI,domainMask,GAdomains,...
            GAvertices] = ROIselector(UserSettingsGUIvalue, UserSettingsGUIstring);
        
        UserSettingsGUIvalue.domainMask = domainMask;
        UserSettingsGUIvalue.GAdomains = GAdomains;
        UserSettingsGUIvalue.domainROI = domainROI;
        UserSettingsGUIvalue.GAvertices = GAvertices;
        UserSettingsGUIvalue.ROI = ROI;
        UserSettingsGUIvalue.LIF_ROI = rect;
        %UserSettingsGUIvalue.rotation = rotation;
        UserSettingsGUIvalue.angle_rotation = angle_rotation;
        %Save the ROIs to the runtime UserSettings
        save('.\Temporary\UserSettingsGUItemp.mat','UserSettingsGUIvalue','UserSettingsGUIstring')
        %Save the altered UserSettings to the UserSetings file in the video folder
        save([inputVideoPath 'UserSettingsGUI.mat'],'UserSettingsGUIvalue','UserSettingsGUIstring')
    end
catch
    disp('GUI closed before it could save the user settings.')  
    return
end

%% Image processing
if(UserSettingsGUIvalue.imageProcessingCheckBox)
    disp('******************************************************************')
    disp('IMAGE PROCESSING**************************************************')
    image_processing
    
end


%% Image Analysis
if(UserSettingsGUIvalue.imageAnalysisCheckBox)
    disp('******************************************************************')
    disp('IMAGE ANALYSIS****************************************************')
    image_analysis
end

%% END
%delete('.\Temporary\*')
disp('******************************************************************')
disp('script terminated successfully')
