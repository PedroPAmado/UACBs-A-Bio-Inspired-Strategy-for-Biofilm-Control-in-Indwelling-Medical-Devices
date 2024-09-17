function [ image ] = plotROIs( image,UserSettingsGUIvalue)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%Retreive Variables from UserSettingsGUI
manualROI_Selection = UserSettingsGUIvalue.roiSelection;
lightIntensityFilter = UserSettingsGUIvalue.LIF;
synchronizationdef=UserSettingsGUIvalue.synchronization;
Ndomains = UserSettingsGUIvalue.nDomainsEdit;
try %Safety mechanism for code versions older than release 2.2
    NGAdomains = UserSettingsGUIvalue.GAdomainsEdit;
    GAvertices = UserSettingsGUIvalue.GAvertices;
catch
    NGAdomains = 0;
end

Rotation = UserSettingsGUIvalue.Rotation;
angle_rotation = UserSettingsGUIvalue.angle_rotation;
rect = UserSettingsGUIvalue.LIF_ROI;
ROI = UserSettingsGUIvalue.ROI;
domainROI = UserSettingsGUIvalue.domainROI;


if(Rotation)
    image = imrotate(image,angle_rotation,'bilinear','crop');
end
imageOriginal =cat(3, image, image, image);
%figure
image = imcrop(image,ROI);

col1 = ROI(1);
col2 = col1 + ROI(3);
row1 = ROI(2);
row2 = row1 + ROI(4);

if(lightIntensityFilter || synchronizationdef )
    image = insertShape(image,'Rectangle',rect,'Color','yellow','LineWidth',2);
end
if(Ndomains>0)
    for i = 1:Ndomains
        image = insertShape(image,'Rectangle',domainROI(i,:),'Color','red','LineWidth',2);
        image = insertText(image,domainROI(i,1:2), sprintf('%d',i),...
            'BoxOpacity',0.4,'BoxColor','red','FontSize',15,'Font','LucidaBrightDemiBold','TextColor','blue');
    end
end
if(NGAdomains>0)
    for i=1:NGAdomains
        image = insertShape(image,'Polygon',GAvertices{i},'Color','magenta','LineWidth',2);
        image = insertText(image,GAvertices{i}(1:2), sprintf('%d',(i+Ndomains)),...
            'BoxOpacity',0.4,'BoxColor','magenta','FontSize',15,'Font','LucidaBrightDemiBold','TextColor','blue');
    end
end
imageOriginal(row1:row2,col1:col2,:) = image;
if(manualROI_Selection || ~isempty(ROI))
   image = insertShape(imageOriginal,'Rectangle',ROI,'Color','blue','LineWidth',2);
end


end

