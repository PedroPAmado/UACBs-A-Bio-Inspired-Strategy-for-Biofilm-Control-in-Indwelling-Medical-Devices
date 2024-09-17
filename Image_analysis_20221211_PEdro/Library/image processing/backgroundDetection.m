function [ background] = backgroundDetection( FrameMat, UserSettings)
%BACKGROUNDDETECTION: This function detects the background of the frames
%using the method specified in UserSettings. The different methods are
%listed below:
%Different background removal methods:
%   - first     :subtract the first frame as bg from all subsequent ones
%   - min       :subtract minimum background calculated over all frames
%   - avg       :subtract the average of all frames as background
%   - run_avg   :subtract the running average
%   - gausMix   :use gaussian mixture models for foreground detection

%Authors: ARTORG-Center, CVE-Group;
% - Damiano Schirinzi
% - Aurelia Bucciarelli

%Progress bar
h = waitbar(0,'Progress BG Detection:');

switch UserSettings.method_bg_rem
    case 'first'
       background = FrameMat(:,:,1); %First frame
       
    case 'min'
       for i=1:UserSettings.frame_step:size(FrameMat,3)
           if i==1
               background = FrameMat(:,:,i);
           else
               background = min(FrameMat(:,:,i),background); %Takes minimum intesity of every pixel
           end
           waitbar(i/size(FrameMat,3),h,sprintf('Progress BG Detection: %2.0f %%', i/size(FrameMat,3)*100));
       end
       
    case 'avg'
        %Convert images to double precision to do arithmetics
        if(UserSettings.frame_step == 1)
            background = sum(double(FrameMat),3);
            background = uint8(background/size(FrameMat,3));
        else
            for i=1:UserSettings.frame_step:size(FrameMat,3)
                if i==1
                    background = double(FrameMat(:,:,i));
                else
                    background = double(FrameMat(:,:,i))+background;
                end
                waitbar(i/size(FrameMat,3),h,sprintf('Progress BG Detection: %2.0f %%', i/size(FrameMat,3)*100));
            end
            %Convert back to precision uint8 
            background = uint8(background/length(1:UserSettings.frame_step:size(FrameMat,3)));
        end
        
    case 'run_avg'
        %Conversion of precision for averaging
        FrameMat = double(FrameMat);
        %Memory allocation
        background = zeros(size(FrameMat));
        run_sum= zeros(size(FrameMat));
        %Averaging of first n-frames where n is the averaging frame number
        run_sum(:,:,1) = FrameMat(:,:,1);
        for i=2:UserSettings.runAvgFrameRange
            run_sum(:,:,i) = FrameMat(:,:,i)+ run_sum(:,:,(i-1));
        end
        run_sum(:,:,(UserSettings.runAvgFrameRange+1)) = run_sum(:,:,UserSettings.runAvgFrameRange);
        %Running through all frames continuing to just average n-frames
        for i=1:size(FrameMat,3)
            if(i<=UserSettings.runAvgFrameRange+1)
                background(:,:,i)=run_sum(:,:,(UserSettings.runAvgFrameRange+1))/UserSettings.runAvgFrameRange;
            else
                run_sum(:,:,i) = run_sum(:,:,(i-1))+FrameMat(:,:,(i-1))-FrameMat(:,:,(i-1-UserSettings.runAvgFrameRange));
                background(:,:,i) = run_sum(:,:,i)/UserSettings.runAvgFrameRange;
            end
            waitbar(i/size(FrameMat,3),h,sprintf('Progress BG Detection: %2.0f %%', i/size(FrameMat,3)*100));
        end
        FrameMat = uint8(FrameMat);
        background = uint8(background); 
       
    case 'gausMix'
        %Initialization of foregroundDetector:
        % - Number of training frame could be changed
        % - Number of gaussians as well if there are more modes in the background
        foregroundDetector = vision.ForegroundDetector('NumGaussians', 2,'NumTrainingFrames', 150);
        %Memory allocation
        foreground = zeros(size(FrameMat));
        for i=1:size(FrameMat,3)
            %foreground is a logical mask: true for foreground pixels
            foreground(:,:,i) = step(foregroundDetector,FrameMat(:,:,i));
            waitbar(i/size(FrameMat,3),h,sprintf('Progress  BG Detection: %2.0f %%', i/size(FrameMat,3)*100));
        end
        background = uint8(~foreground);
        
    otherwise
         error('ERROR: invalid backgorund removal method selected!')
         
end
delete(h) %Delete progress bar
end

