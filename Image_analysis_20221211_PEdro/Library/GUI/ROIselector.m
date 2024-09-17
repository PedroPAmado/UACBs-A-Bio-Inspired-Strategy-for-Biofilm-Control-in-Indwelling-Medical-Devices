 function [angle_rotation,ROI,rect,domainROI,domainMask,...
     GAdomains,GAvertices] = ROIselector( ...
     UserSettingsGUIvalue, UserSettingsGUIstring)
%ROIselector This function contains all the GUIs to select the image
%processing ROI, the LIF ROI and the image analysis ROIs
%   Note: This function works as subfunction of the main.m script and is
%   not inteded as a run alone function. 

%% Variables for ROI selection
manualROI_Selection = UserSettingsGUIvalue.roiSelection;
synchronizationdef = UserSettingsGUIvalue.synchronization;
lightIntensityFilter = UserSettingsGUIvalue.LIF;
Rotation = UserSettingsGUIvalue.Rotation;
Ndomains = UserSettingsGUIvalue.nDomainsEdit;
NGAdomains = UserSettingsGUIvalue.GAdomainsEdit;
inputVideoPath = UserSettingsGUIvalue.selectedVideodir;
inputVideoName = UserSettingsGUIstring.SelectedVideoList{1};
startframe = UserSettingsGUIvalue.startframeEditField;
video = [inputVideoPath inputVideoName];
image = imread(video,startframe);

if(size(image,3)>1)
    image = rgb2gray(image);
end

%% ROI and rotation and LIF selection
if(manualROI_Selection)
    figure
    imshow(image,[]) %Load the start frame
    
    %%Rotation Dialogue
    if Rotation
    angle_rotation=Rotate(image);
    else 
        angle_rotation=0;
    end
    %rotate image
    image_rot = imrotate(image,angle_rotation,'bilinear','crop');
    imshow(image_rot,[])

    
    %%Manual ROI and LIF_ROI selection   
    uiwait(msgbox('Select ROI for image processing.'));
    [image_ROI,ROI] = imcrop; %Save manually selected ROI in variable ROI
    imshow(image_ROI,[])
    uiwait(msgbox('Select the reference region for Linght Intensity Filtering or Syncronization'));
    [~,rect] = imcrop; %Output rect are the cooardinates for LIF ROI    
    image = image_ROI; %Needed for domainSelection
   
else
    %Check if there are previous ROI Settings
    if(isempty(UserSettingsGUIvalue.ROI)) %No previous ROI settings
        ROI = []; % the all picture
        rect = [];
        
    else %Use the old ROI settings saved in UserSettings
        angle_rotation = UserSettingsGUIvalue.angle_rotation;
        image_rot = imrotate(image,angle_rotation,'bilinear','crop');
        ROI = UserSettingsGUIvalue.ROI;
        image_ROI = imcrop(image_rot,ROI);      
        image = image_ROI; %Needed for domainSelection
        figure
        imshow(image,[])
    end
    
    if(isempty(UserSettingsGUIvalue.LIF_ROI))
        %Need a LIF_ROI if LIF is on, even if manualROI selection is off.
        if(lightIntensityFilter || synchronizationdef )
            figure
            imshow(image,[]) %Load the start frame
            uiwait(msgbox('Select the reference region for Linght Intensity Filtering'));
            [~,rect] = imcrop; %Output rect are the cooardinates for LIF ROI
        end
    else %Use the old ROI settings saved in UserSettings
        rect = UserSettingsGUIvalue.LIF_ROI;
    end
end

%% Main domain selection
if(Ndomains>0)
    domainROI = zeros(Ndomains,4); %domainROI(i,:) contains the imcrop coordinates for domain i
    uiwait(msgbox(['Select ' num2str(Ndomains) ' domains for image analysis']));
    for i=1:Ndomains      
        [~,domainROI(i,:)] = imcrop;
        hold on
        rectangle('Position',domainROI(i,:),'EdgeColor','r','LineWidth',2);        
    end
else
   domainROI = [0 0 size(image,2) size(image,1)]; 
end

%Create a domain Mask
domainMask = zeros([size(image(:,:,1)) size(domainROI,1)]);
for i=1:size(domainROI,1)
    xyCoord = bbox2points(domainROI(i,:));
    domainMask(:,:,i) = uint8(roipoly(image,xyCoord(:,1),xyCoord(:,2)));

end

%% Additional domain selection
GAdomains = zeros([size(image(:,:,1)) NGAdomains]);
if(NGAdomains>0)
    uiwait(msgbox(['Select ' num2str(NGAdomains) ' additional domains for general analysis']));
    for i=1:NGAdomains      
        [GAdomains(:,:,i),x,y] = roipoly; %Save Mask in GAdomains
        AA = [x'; y'];
        GAvertices{i} = reshape(AA,1,2*size(AA,2)); %Save vertices in form [x1 y1 x2 y2..]
        hold on
        plot(x,y,'m','LineWidth',2)
    end
else
   GAvertices = [0 0 size(image,2) size(image,1)]; 
end

close(gcf) % Close open figure



end

