function [FrameMat, ProcessingSteps] = imageBinarization(FrameMat, UserSettings, ProcessingSteps)
%IMAGEBINARIZATION: This script binarizes the frames either manually by a
%threshold set in USerSettings or automized using the function
%Autothresholder from the vision package which is based on Otsu's Method

%Authors: ARTORG-Center, CVE-Group;
% - Damiano Schirinzi

if(UserSettings.autoBinarization)
    %Initialize auto-binarizator
    binarizator = vision.Autothresholder('ThresholdOutputPort',true,'EffectivenessOutputPort',true);
    %Allocate memory:
        %Threshold used by binarizator
        TH = zeros(1,size(FrameMat,3));
        %Efficiency metric: 0 attained by image with just one grayscale
        %level and 1 for perfectly binary images.
        EMETRIC = zeros(1,size(FrameMat,3));
    %Loop through all frames
    for i=1:size(FrameMat,3)
        %Autobinarization
        [FrameMat(:,:,i), TH(i), EMETRIC(i)] = step(binarizator, FrameMat(:,:,i));     
    end
    %Saving processing step
    ProcessingSteps.binarization = true;
    %Save thresholds and efficiency metrics to ProcessingSteps
    ProcessingSteps.autothreshold = TH;
    ProcessingSteps.emetric = EMETRIC;
else %Manual binarization
    %Get binarization threshold
    threshold_BW = UserSettings.threshold_BW;
    for i=1:size(FrameMat,3)
        %use working copy of frame
        b = FrameMat(:,:,i);
        %Binarization using threshold
        b(b>threshold_BW)=255;
        b(b<=threshold_BW)=0;
        %exchange binarized frame into FrameMat
        FrameMat(:,:,i) = b;
    end
    %Saving processing script
    ProcessingSteps.binarization = true;
end

end

