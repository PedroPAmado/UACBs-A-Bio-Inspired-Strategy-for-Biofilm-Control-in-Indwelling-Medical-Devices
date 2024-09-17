function [FrameMat, ProcessingSteps] = noiseFiltering(FrameMat, UserSettings, ProcessingSteps  )
%NOISEFILTERING Performs noise filtering on all frames based on
%mathematical morphology.
%   The noise filter will perform an erosion step on all areas with a
%   filtering Element specified in UserSettings and subsequently will
%   reopen the remaining areas either by reconstructing them to their
%   original size or by dilating them with the same filtering element as in
%   the first step.

%Authors: ARTORG-Center, CVE-Group;
% - Damiano Schirinzi

%Waitbar
h = waitbar(0,'Progress NF:');

%Loop through all Frames
for i=1:size(FrameMat,3)
    SE = UserSettings.NFelement;
    if(UserSettings.reconstruction)   
        image_eroded = imerode( FrameMat(:,:,i) ,SE); %erosion step
        %dilate image by reconstruction
         FrameMat(:,:,i) = imreconstruct(image_eroded, FrameMat(:,:,i) );
    else 
        %dilate image with Element SE
         FrameMat(:,:,i)  = imopen( FrameMat(:,:,i) , SE);
    end
    
    
    %Update waitbar
    waitbar(i/size(FrameMat,3),h,sprintf('Progress NF: %2.0f %%', i/size(FrameMat,3)*100));
end
%Saving processing step
ProcessingSteps.NF = true;

delete(h) %Delete waitbar
end

