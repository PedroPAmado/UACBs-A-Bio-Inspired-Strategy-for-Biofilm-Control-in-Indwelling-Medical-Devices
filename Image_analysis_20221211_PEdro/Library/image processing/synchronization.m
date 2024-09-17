function [Flash,Flash_info,ProcessingSteps] = synchronization(FrameMat,startframe,endframe,FrameRate,inputVideoPath,inputVideoName,UserSettings,ProcessingSteps)


if isempty(UserSettings.LIF_ROI)
    fprintf(2,'\nNo ROI for synchronization selected! No synchronization done\n')
    return
end

try
    load([inputVideoPath inputVideoName(1:end-4) '_Pressure_pump_data.mat']);
    Nikon= ImportNikon([inputVideoPath inputVideoName(1:end-4) '_Nikon.xlsx'], 1, [startframe+2 endframe+2]);
    Nikon_info = importNikon_info([inputVideoPath inputVideoName(1:end-4) '_Nikon.xlsx'],1,[1 1]);

catch
     fprintf(2,'\nNo Info form pressure pump and nikon found! No synchronization done \n')
    return
end 
%Loop thorugh frames
for i=1:size(FrameMat,3)
    %Mean light intensity of ROI in current frame
    Im = imcrop(FrameMat(:,:,i),UserSettings.LIF_ROI);
    meanImg(i,1)=mean2(Im);
end

[~,locs]=findpeaks( meanImg,'MinPeakProminence',4);
n=0;
%%
figure 
plot(meanImg)
hold on
[pk,locs]=findpeaks( meanImg,'MinPeakProminence',4);
plot(locs,pk,'x','MarkerSize',12)
ylim([44 54])
ylabel('Mean intensity','interpreter', 'Latex')
xlabel('Frame','interpreter', 'Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',30);
set(gcf, 'Position', [25 100 1000 750] );
exportgraphics(gcf,[inputVideoPath inputVideoName(1:end-4) '_peak.png'],'Resolution',1000) 
close all 
%%
for i=40:-1:-20
   difference(1+n,1)= meanImg(locs-i,1)-meanImg(locs-i-1,1);
   difference(1+n,2)=locs-i;  
   n=n+1;
end
[~,maxind]=max(difference(:,1));
[~,minind]=min(difference(:,1));
framemicroscope=difference(maxind,2);
minframe=difference(minind,2);
timeflashon=(minframe-framemicroscope)/FrameRate;
timepumpflash=Experiment_data(end,1);
timemicroscopestart=Nikon(1,1);
timemicroscopeflash=Nikon(framemicroscope,1);
%Saving processing step
Flash=[timemicroscopestart,timemicroscopeflash,framemicroscope,timeflashon,timepumpflash];
Flash_info={'Timestart microscope','Flash time microscpe','Frame flash microscope','Flash on time microscope (40ms)','Frash pressure pump'};
ProcessingSteps.synchronization = true;
   

end