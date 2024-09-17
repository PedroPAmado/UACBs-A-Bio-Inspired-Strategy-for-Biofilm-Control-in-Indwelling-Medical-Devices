function [ FrameMat,ProcessingSteps ] = backgroundRemoval( FrameMat, UserSettings, ProcessingSteps, background )
%BACKGROUNDREMOVAL: Removes the background from each frame depending on the
%background removal method specified in UserSettings. The different methods
%are listed below:
%Different background removal methods:
%   - first     :subtract the first frame as bg from all subsequent ones
%   - min       :subtract minimum background calculated over all frames
%   - avg       :subtract the average of all frames as background
%   - run_avg   :subtract the running average
%   - gausMix   :use gaussian mixture models for foreground detection

%Authors: ARTORG-Center, CVE-Group;
% - Damiano Schirinzi

switch UserSettings.method_bg_rem
    case {'first','min','avg'}
        %Loop thorugh all frames
        for i=1:size(FrameMat,3)
            %Subtract background from frame
            FrameMat(:,:,i) =  FrameMat(:,:,i)-background;
        end
        %Saving processing step
        ProcessingSteps.backgroundRemoval = true;
    case 'run_avg'
        %Loop thorugh all frames
        for i=1:size(FrameMat,3)
            %Subtract background from frame
            FrameMat(:,:,i) =  FrameMat(:,:,i)-background(:,:,i);
        end
        %Saving processing step
        ProcessingSteps.backgroundRemoval = true;
        
    case 'gausMix'
        %Foreground is the opposite of the backgroundmask
        foreground = uint8(~background);
        if(UserSettings.binarization) %Binarization step follows later
            %Loop thorugh all frames
            for i=1:size(FrameMat,3)
                %Multiply frame times binary foreground mask
                FrameMat(:,:,i) =  FrameMat(:,:,i).*foreground(:,:,i);
            end
        else %Take foreground mask from gaussian mixture model as already 
            %binarized frame.
            %Loop thorugh all frames
            for i=1:size(FrameMat,3)
                %Set frame equal to foreground mask
                FrameMat(:,:,i) = foreground(:,:,i);
            end
        end
        %Saving processing step
        ProcessingSteps.backgroundRemoval = true;
        
    otherwise
         error('ERROR: invalid backgorund removal method selected!')
         
end
end

