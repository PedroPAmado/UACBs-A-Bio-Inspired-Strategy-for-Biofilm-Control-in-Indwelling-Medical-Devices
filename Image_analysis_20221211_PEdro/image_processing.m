%% image_processing.m
%@brief: This script does the image processing of video material gathered
%from the microscope. In the USER INPUT section you can set all the 
%processing steps that should be performed on the frames.
%This script allows you to select a specified area of the video to be
%processed. You can also rotate the image if needed.
%It allows for inversion of grayscales
%It allows to perform background detection and removal with several 
%techniques. 
%It allows for light intensity and noise filtering.
%It will output the processed frames in matrix form .mat and or will create
%a processed video if desired.

%Authors: ARTORG-Center, CVE-Group;
% - Francesco Clavica
% - Alberto Mantegazza
% - Damiano Schirinzi
% - Aurelia Bucciarelli


%% Load GUI Data

%Load Inputs from GUI into workspace and delete file containing them 
load('.\Temporary\UserSettingsGUItemp.mat')


%% Insert GUI INPUT into script variables
automated = true;
video=([inputVideoPath inputVideoNames{1}]);

% Frame selection
startframe = UserSettingsGUIvalue.startframeEditField;
frame_step = UserSettingsGUIvalue.averageEdit;   
endframe = UserSettingsGUIvalue.endframeEditField;
framerange = endframe-startframe+1; 
frameRate=UserSettingsGUIvalue.framerateEdit;

% Realignment
realignmentdef=UserSettingsGUIvalue.Realignment;

%Synchronization
synchronizationdef=UserSettingsGUIvalue.synchronization;


%Grayscale Inversion:
%(for the processing need bright RBCs on dark background)
grayScaleInversion = UserSettingsGUIvalue.grayscaleInversion;
removeBackground = UserSettingsGUIvalue.backgroundRemoval;
%Three back ground removal methods:
%   - first     :subtract the first frame as bg from all subsequent ones
%   - min       :subtract minimum background calculate over all frames
%   - avg       :subtract the average of all frames as backgorung
%   - run_avg   :subtract the running average
%   - gausMix   :use gaussian mixture models for foreground detection

    if strcmp(UserSettingsGUIvalue.backgroundPopUp,'first frame') %'first frame'
        method_bg_rem = 'first';
    elseif strcmp(UserSettingsGUIvalue.backgroundPopUp,'minimum') %'minimum'
        method_bg_rem = 'min';
    elseif strcmp(UserSettingsGUIvalue.backgroundPopUp,'average') %'average'
        method_bg_rem = 'avg';
    elseif strcmp(UserSettingsGUIvalue.backgroundPopUp,'running average') %'running average'
        method_bg_rem = 'run_avg';
    elseif strcmp(UserSettingsGUIvalue.backgroundPopUp,'gaussian mixture model') %'gaussian mixture model'
        method_bg_rem = 'gausMix';
    else
       disp('Error: Wrong background removal method detected')     
    end

plotting_background = false;  
    %Settings for run_avg:
    runAvgFrameRange = UserSettingsGUIvalue.runAvgEdit;
    
%Light Intensity Adjustment
lightIntensityFilter = UserSettingsGUIvalue.LIF;
manualROI_Selection = UserSettingsGUIvalue.roiSelection;

%Threshold for binary image conversion to black & white
binarization= UserSettingsGUIvalue.binarization; 
binarizationPopUp = UserSettingsGUIvalue.binarizationPopUp;    
if strcmp(binarizationPopUp,'manual')
    autoBinarization = false; 
else
    autoBinarization = true;
end
threshold_BW = UserSettingsGUIvalue.binarizationTHedit;        
    %Threshold for manual binarization
    %manual binarization does not work with gausMix, as the image
    %intensities are scaled down to range [0 1], therefore need to
    %normalize also the manual threshold:
    if(strcmp(method_bg_rem,'gausMix'))
        threshold_BW = threshold_BW/255;
    end

%Noise Filtering with Mathematical Morphology Filter
noiseFilter = UserSettingsGUIvalue.NF;
reconstruction = UserSettingsGUIvalue.reconstruction;      %Opening by reconstruction?

    if strcmp(UserSettingsGUIvalue.filteringElementPopUp,'disk')
        elementForm = 'disk';       %Form of filtering Element
    else
        disp('Error: chosen wrong Filtering Element!')
    end
    
elementRadius = UserSettingsGUIvalue.NFelementRadiusEdit; % in [px]
elementorientation = 90;    %just needed for nonsymmetric shapes.
NFelement = strel(elementForm,elementRadius);
NFproperties = struct('elementForm',elementForm, 'elementRadius', elementRadius);

imagefill = UserSettingsGUIvalue.IMFILL;


% Automated
if(automated)
    matFileName = {};
    for i=1:length(inputVideoNames)
        matFileName{i} = [inputVideoPath inputVideoNames{i}(1:end-4) '.mat']; %
        FlashFileName{i} = [inputVideoPath inputVideoNames{i}(1:end-4) '_Flash.mat']; %
    end
end

%Want to convert created images to avi? 
%   Yes:Set pearameters below
%   No: Set conversion2avi to false
conversion2avi = UserSettingsGUIvalue.outputVideo;
    PlayVideo = false;       %Play video after saving it?
    addFrameNumber = false;  %true if want frame number to be added to images
    % Automated
    if(automated)
        outputVideoName = {};
        for i=1:length(inputVideoNames)
            outputVideoName{i} = [inputVideoPath inputVideoNames{i}(1:end-4) '_processed.avi'];
        end
    end
    
%Visual Control by plotting after runtime?
visualControl = UserSettingsGUIvalue.visualControl;
%%%%%%%%%%%%%%%%%%%%END USER INPUT
%% ROI selection

rotation = UserSettingsGUIvalue.Rotation;
angle_rotation = UserSettingsGUIvalue.angle_rotation;
ROI = UserSettingsGUIvalue.ROI;
rect = UserSettingsGUIvalue.LIF_ROI;
domainROI = UserSettingsGUIvalue.domainROI;
domainMask = UserSettingsGUIvalue.domainMask; 
FrameRate=UserSettingsGUIvalue.framerateEdit;


for k=1:length(inputVideoNames) %Loop over all videos contained in inputVideoNames
    %Display Log outputs on command window:
    disp(['Processing of:' inputVideoNames{k}])
    disp('-------------------------------------')
    %threshold_BW = threshold_BW_vec(k);  %% REMOVE (used for thresholdstudy)
    
    %Load video
    video=([inputVideoPath inputVideoNames{k}]);
    %Run time tracker
    tstart = tic; 
    %Create a struct with all USER settings used to append to FrameMat:
    
    %User Settings as chosen above and useful inputs for subsequent analysis: 
    UserSettings = struct('inputVideoPath',inputVideoPath,...
        'inputVideoName',inputVideoNames{k},'frame1',imread(video,startframe),...
        'startframe',startframe,...
        'frame_step',frame_step,'framerange',framerange,'endframe',endframe,...
        'manualROI_Selection',manualROI_Selection,'ROI',ROI,'LIF_ROI',rect,...
        'Ndomains',UserSettingsGUIvalue.nDomainsEdit,...
        'domainROI',domainROI,'domainMask',domainMask,...
        'method_bg_rem',method_bg_rem,...
        'runAvgFrameRange',runAvgFrameRange,'lightIntensityFilter',lightIntensityFilter,'synchronization',synchronizationdef,...
        'binarization',binarization,'autoBinarization',autoBinarization,'threshold_BW',threshold_BW,...
        'noiseFilter',noiseFilter,'reconstruction',reconstruction,...
        'NFproperties',NFproperties,'NFelement',NFelement,...
        'RBCfill',imagefill,...
        'conversion2avi',conversion2avi,'outputVideoName',outputVideoName{k},...
        'frameRate',frameRate,...
        'rotation',rotation,'rotationAngle',angle_rotation,...
        'grayScaleInversion',grayScaleInversion,'realignment',realignmentdef);
    
    %Tracker of Processing Steps done to the images; must not be equal to
    %the UserSettings if this script is run manually jumping over some
    %functions. But in general this is just an ulterior check. But does
    %also contain the metrics of the autobinarization if needed.
    ProcessingSteps = struct('realignment',false,'synchronization',false,'grayScaleInversion',false,'rotation',false,'LIF',false,...
        'backgroundRemoval',false,'binarization',false,'autothreshold',false,'emetric',false,'NF',false);

    %% Load video to matrix FrameMat
    
    % FrameMat(:,:,i) contains the ith frame 
    disp('Reading in frames from video...')
    h = waitbar(0,'Loading frames:');
    

        %Allocation of memory
        %image = imcrop((imread(video,1)),ROI); %Load and crop image to ROI
        %FrameMat = uint8(zeros([size(image) framerange]));
        for i=1:framerange %Loop through all frames
            j=startframe-1+i;
            image =im2uint8(imread(video,j));
            %Rotation of image
            if rotation
                angle_rotation = UserSettings.rotationAngle;
                image_rot =imrotate(image,angle_rotation,'bilinear','crop');
            else
                image_rot = image;
            end
            %Crop to ROI
            FrameMat(:,:,i) = imcrop(image_rot,ROI);
            waitbar(i/framerange,h,sprintf('Loading frames: %2.0f %%', i/framerange*100));
        end
        %Saving processing step
        ProcessingSteps.rotation = true;

    delete(h)
    FrameMatGrey = FrameMat;
    disp('-------------------------------------')
    
%% Calculation of synchronization
    
    if (synchronizationdef)
        disp('Synchronization')
        [Flash,Flash_info,ProcessingSteps]=synchronization(FrameMat,startframe,endframe,FrameRate, inputVideoPath, inputVideoNames{k},UserSettings,ProcessingSteps);
        disp('-------------------------------------')
    end
    
    %% Realignment
    if (realignmentdef)
        disp('Realignment')
        [FrameMat,ProcessingSteps]=realignment(FrameMat,startframe,endframe, inputVideoPath, inputVideoNames{k},ProcessingSteps);
        disp('-------------------------------------')
    end

    %% Inverion of grayscale and rotation of image
    if(grayScaleInversion)
        disp('Grayscale Inversion ...')
        FrameMat= imcomplement(FrameMat);
        %Saving processing step
        ProcessingSteps.grayScaleInversion = true;
        disp('-------------------------------------')
    end
    
    FrameMatInv = FrameMat;
    
    %% Background Calculation
    if(removeBackground)
        disp('Background detection...')
        [background] = backgroundDetection(FrameMat, UserSettings);
        disp('-------------------------------------')
    end
    %% LIF - Light Intensity Filter
    %Must come after background detection and before the background removal
    %Needs background as reference frame to do the light intensity
    %adjustments
    if(lightIntensityFilter && removeBackground)
        disp('LIF...')
        [FrameMat,ProcessingSteps] = lightIntensityFiltering(FrameMat, UserSettings, ProcessingSteps, background);
        disp('-------------------------------------')
    end
    
    %% Background Removal
    if(removeBackground)
        disp('Background removal...')
        [FrameMat,ProcessingSteps] = backgroundRemoval(FrameMat, UserSettings, ProcessingSteps, background);
        disp('-------------------------------------')
    end
    
    %% Binarization
    if(binarization)
        disp('Binarization...')
        [FrameMat, ProcessingSteps] = imageBinarization(FrameMat, UserSettings, ProcessingSteps);
        disp('------------------------------------')
    end
    
    %% NF - Noise Filtering
    if(noiseFilter)
        disp('Noise Filtering...')
        [FrameMat, ProcessingSteps] = noiseFiltering(FrameMat, UserSettings, ProcessingSteps);
        disp('-------------------------------------')
    end
    
    %% RBC filling
    if(imagefill)
        disp('RBC filling...')
        for i=1:size(FrameMat,3)
            FrameMat(:,:,i) = imfill(FrameMat(:,:,i));
        end
        disp('-------------------------------------')
    end
    
    %% Saving Data
    
            disp(['Saved frames to: ' matFileName{k}]) %save mat
            save(matFileName{k},'FrameMat','FrameMatGrey','UserSettings','ProcessingSteps','-v7.3');
             %save(FlashFileName{k},'Flash','Flash_info','-v7.3');

        if(conversion2avi) %Create Movie
            disp(['Creating movie:' outputVideoName{k}])
            mat2avi(FrameMat,frameRate,outputVideoName{k},PlayVideo)
        end
        disp('-------------------------------------')
 
    
    %% Finish for automated loop
    %Elapsed time
    toc(tstart)
    disp('-------------------------------------')
    disp('-------------------------------------')
end
%% Visual Control
if(visualControl)
    fig = figure(20);
    i= 1;
    stoppedByUser = false;
    title('Visual Control');
    stopButton = uicontrol(fig, 'Style', 'togglebutton', 'String', 'Continue'); %toggle button to stop
    while ((i<=size(FrameMat,3)) && stoppedByUser~=true)
        if(addFrameNumber) %Add box with frame number in upper left corner 
            %of frame before displaying it.
            b = insertText(FrameMat(:,:,i), [10 10], sprintf('frame: %6.0f',i), 'BoxOpacity', 1,'FontSize', 20);
            subplot(1,2,1); imshow(b,[]);
        else
            subplot(1,2,1); imshow(FrameMat(:,:,i),[])
        end
        %Show original frame in right subplot
        subplot(1,2,2); imshow(FrameMatInv(:,:,i),[])
        drawnow
        if (get(stopButton, 'Value')==1)
            break; 
        else
            i = i+1;
        end
    end
    close(gcf)
end

%% END of script
disp('image processing terminated successfully!')







