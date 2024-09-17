function [ FrameMat, ProcessingSteps ] = lightIntensityFiltering( FrameMat, UserSettings, ProcessingSteps, background )
%LIGHTINTENSITYFILTER This function filters out light intensity
%oscillations in the frames by pulling the intensity to the one of the
%background. In detail:
%   The light Intensity filter computes the difference in mean light 
%   intensity between the LIF_ROI of the background and the LIF_ROI
%   in each frame and then adds this difference to the whole frame.

%Authors: ARTORG-Center, CVE-Group;
% - Damiano Schirinzi


diff_lightInt = zeros(1,size(FrameMat,3)); 
if isempty(UserSettings.LIF_ROI)
    fprintf(2,'\nNo ROI for LIF selected! No light intensity filtering done\n')
    return
end

switch UserSettings.method_bg_rem
    case {'first','min','avg'}
        %Mean light intensity of ROI in backgorund 
        IL_ref = imcrop(background,UserSettings.LIF_ROI);
        mean_ref = mean2(IL_ref);
        %Loop thorugh frames
        for i=1:size(FrameMat,3)
            %Mean light intensity of ROI in current frame
            IL = imcrop(FrameMat(:,:,i),UserSettings.LIF_ROI);
            meanImg=mean2(IL);
            %Difference in intensity to background
            diff_lightInt(i) = meanImg - mean_ref;
            %Subtraction of difference in the whole frame
            FrameMat(:,:,i) = FrameMat(:,:,i) - diff_lightInt(i);
        end
        %Saving processing step
        ProcessingSteps.LIF = true;
    case 'run_avg'
        for i=1:size(FrameMat,3)
            %Mean light intensity of ROI in backgorund
            IL_ref = imcrop(background(:,:,i),UserSettings.LIF_ROI);
            mean_ref = mean2(IL_ref);
            %Mean light intensity of ROI in current frame
            IL = imcrop(FrameMat(:,:,i),UserSettings.LIF_ROI);
            meanImg=mean2(IL);
            %Difference in intensity to background
            diff_lightInt(i) = meanImg - mean_ref;
            %Subtraction of difference in the whole frame
            FrameMat(:,:,i) = FrameMat(:,:,i) - diff_lightInt(i);
        end
        %Saving processing step
        ProcessingSteps.LIF = true;
    case 'gausMix'
        %No reference background frame available to do light intensity filtering on
        disp('Alert: No ligth intensity filtering for gaussian Mixture Model applied.')
    otherwise
        error('ERROR: invalid backgorund removal method selected!')
end

disp('LIF statistics:')
fprintf('max intensity difference: %3.0f \n',max(diff_lightInt))
fprintf('min intensity difference: %3.0f \n',min(diff_lightInt))

end

