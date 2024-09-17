function [ ] = haematocrit( inputVideoPath, inputVideoNames, UserSettingsGUI)
%% Haematocrit.m
%@brief: This script estimates the haematocrit from 2D microscope images of
%capillaries with known dimensions.//
%This script analyzes automaticlly the results of different experiments
%(e.g. different feed Haematocrits) in an outer loop j and different
%flow velocities for each experiment in an inner loop k
%In the end Plots will be made for
%-each experiment and each flow velocity containing:
% Cell area Histogram; Number of cells per frame;
% Number of particles detected by PTV vs Haematocrit script; Haematocrit
%   - if PTV matrix not available for the analyzed experiment the 3rd
%   subplot will be replaced by the histogram of the equiv. cell diameters
%-for each experiment we will get two plots:
% The mean Haematocrit over the different velocities
% The percentage of Particles detected by the PTV
%   This second plot is done only if PTV data was available
%-At each run it will plot all PTV performances of each experiment over all
% velocities in a single plot.
%Output:
%this script will Output a matrix HtcResults if saveOutputMat is true. It will
%contain all the important information about the analysis done.

%Authors: ARTORG-Center, CVE-Group;
% - Damiano Schirinzi
% - Aurelia Bucciarelli

%To Add:
%- Watershedding to separate cells mathematically


%% USER INPUT %%%%%%%
%Images
closeSingleImages = false;
%Folder Path
FolderPath = inputVideoPath;
%input Mat filename
inputMatName = inputVideoNames;


%% GUI Input %%%%%%%%%%%

UserSettingsGUIvalue = UserSettingsGUI{2};

saveOutputMat = UserSettingsGUIvalue.saveMatHaemCheckBox;
saveImages = UserSettingsGUIvalue.saveFiguresHaemCheckBox;
Hf = UserSettingsGUIvalue.feedHaemEdit/100; %Feed Haematocrits
channel_width = UserSettingsGUIvalue.channelWidthEdit; %[micro meter]
channel_heigth = UserSettingsGUIvalue.channelHeightEdit; %[micro meter]
useDomainWidth = UserSettingsGUIvalue.domainWidthCheckBox;
domainROI = UserSettingsGUIvalue.domainROI;
Ndomains = size(domainROI,1);
calibration =UserSettingsGUIvalue.calibrationEdit; %[micrometer/pixel]
bloodtipe=UserSettingsGUIvalue.TipeofBloodDropDown;
Diameter=UserSettingsGUIvalue.DiametrerdefinitonDropDown;

%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
%Relative Haematocrit formula for tubes from Pries_1990
Ht = @(Hd,D) Hd*(Hd + (1-Hd)*(1+1.7*exp(-0.415*D)-0.6*exp(-0.011*D))); %H in [vol fraction] and D in [micro meter]

%Average cell volume
switch bloodtipe
    case 'Human'
        cell_Vol = 85; %[micro meter ^3]/cell %Taken from: tomaiuolo_2012
    case 'Pig (Namdee)'
        cell_Vol = 61; %[micro meter ^3]/cell %Taken from: 17_Namdee_2015
    case 'Pig (Salvioli)'
        cell_Vol = 66; %[micro meter ^3]/cell %Taken from bookmark: https://link.springer.com/article/10.1007%2FBF02537121 (Salvioli)
end
%% Calculations

for i=1:length(inputVideoNames)
    %the Frames are saved in a .mat file with the same same as the video
    %without the extension .avi.
    FrameMatList{i} = [inputVideoNames{i}(1:(end-4)) '.mat'];
end

for k=1:length(FrameMatList) %Loop over different flow velocities
    disp(['Haematocrit Calculations of: ' FrameMatList{k}])
    %Load Matrix with processed Frames
    FrameMat = load([FolderPath FrameMatList{k}]);
    Frames = FrameMat.FrameMat;
    %Allocate memory
    trueMeanA_vec = zeros(Ndomains,1); %Mean Area of cells for each domain
    trueMeanD_vec = zeros(Ndomains,1); %Mean Equivalent Diameter of cells for domain
    haematocrit_mat = zeros(Ndomains,size(Frames,3)); %Haematocrit at each frame for each domain
    haematocrit2D_mat = zeros(Ndomains,size(Frames,3)); %2D haematocrit at each frame for each domain
    Ncells_mat = zeros(Ndomains,size(Frames,3)); %Number of cells per frame for each domain
    Hpries_vec = zeros(Ndomains,1); %Pries Haematocrit for each domain
    for j=1:Ndomains
        %Reload Frames for next domain analysis
        Frames = FrameMat.FrameMat;
        string =['            Domain ' num2str(j)];    
        FramesCut = [];
        for i1=1:size(Frames,3)
            FramesCut(:,:,i1) = imcrop(Frames(:,:,i1),domainROI(j,:));
        end
        Frames = FramesCut;
        
        %Determine remaning channel dimensions
        if(useDomainWidth)
            channel_width = domainROI(j,4)*calibration;%[micro meter]
            string = [string '; channel width = ' num2str(channel_width) ' \mum'];
        end
        disp(string)
        channel_length = size(Frames,2)*calibration;   %[micro meter]
        channel_Vol = channel_width*channel_length*channel_heigth; %[micro meter ^3]
        %2D channel size in [px]
        channel_Area_p = size(Frames,2)*channel_width/calibration; %[px]
        %%%%%%%%%%%%%%%%%%%%
        
        %Equivalent tube diameter
        switch Diameter
            case 'Based on equality'
                %Based on equality of cross sectional area
                area_cross = channel_width*channel_heigth;
                D_equ = 2*sqrt(area_cross/pi);
            case 'Hidraulic'
                %Hidraulic diameter for rectangular duct (actually only accurate for turbulent flows!)
                D_equ = 2*channel_width*channel_heigth/(channel_width+channel_heigth);
                %Complex variable analysis by cross section mapping
                % In principle, any duct cross section can be solved analytically for the laminar flow
                % velocity distribution, volume flow, and friction factor. This is because any cross section
                % can be mapped onto a circle by the methods of complex variables
        end
        %Preallocate space
        area = [];
        diameter = [];
        meanAreaframe = zeros(1,size(Frames,3));
        
        %Loop through all frames
        for img=1:size(Frames,3)
            
            %read frame and binarize to logical values [0,1]
            frame = Frames(:,:,img) >100;
            %Calculate properties of all areas in the frame
            stats = regionprops('table',frame,'Area','Centroid','EquivDiameter');
            %Calculate the cell aerea in the current frame
            cellAreaframe(img) = sum(stats.Area);%[px]
            %Gather single areas of bloobs detected in the frames
            area = [area ; stats.Area];
            %Gather single diameters of bloobs detected in the frames
            diameter = [diameter; stats.EquivDiameter];
            
        end
        
        %Get true single cell area mean by looking at histogram of cell areas
        %instead of taking mathematical mean of the cell areas.
        [Nhist,Xhist]=hist(area,0:10:1000);
        index = find(Nhist == max(Nhist)); %true mean is where histogram peaks
        trueMeanA = Xhist(index);
        
        %Get true single cell perimeter mean by looking at histogram of cell perimeters
        %instead of taking mathematical mean of the perimeters.
        [Nhist,Xhist]=hist(diameter,0:1:35);
        index = find(Nhist == max(Nhist)); %true mean is where histogram peaks
        trueMeanD = Xhist(index);

        %Haematocrit calculation
        Ncells = cellAreaframe/trueMeanA;
        haematocrit = Ncells*cell_Vol/channel_Vol*100;
        haematocrit2D = cellAreaframe/channel_Area_p*100;
        
        %% Plots
        
        figure
        
        subplot(2,2,1)
        hist(area,10:10:1000)
        title('Histogram of cell area')
        xlabel('cell area')
        ylabel('')
        
        subplot(2,2,2)
        hold on
        plot(1:size(Frames,3),cellAreaframe/trueMeanA)
        plot(1:size(Frames,3),mean(cellAreaframe/trueMeanA)*ones(1,size(Frames,3)))
        legend('Htc Cal.','mean Htc Cal.')
        title('Number of cells per frame')
        xlabel('frame')
        ylabel('cells')
        
        
        subplot(2,2,3)
        hist(diameter,0:1:35)
        title('Histogram of cell diameters')
        xlabel('cell diameter')
        ylabel('')
        
        subplot(2,2,4)
        plot(1:size(Frames,3),haematocrit)
        hold on
        plot(1:size(Frames,3),mean(haematocrit)*ones(1,size(Frames,3)))
        plot(1:size(Frames,3),haematocrit2D)
        plot(1:size(Frames,3),mean(haematocrit2D)*ones(1,size(Frames,3)))
        plot(1:size(Frames,3),Ht(Hf,D_equ)*100*ones(1,size(Frames,3)))
        legend('Htc_{exp}','mean(Htc_{exp})','Htc_{2D}','mean(Htc_{2D})','Htc_{Pries} with corrected D')
        title('Haematocrit per frame')
        xlabel('frame')
        ylabel('Haematocrit [%]')
        sgtitle(['Video: ' FrameMatList{k}(1:(end-4)) '; Domain ' num2str(j)])
        
        if(saveImages)
            savefig([FolderPath 'Htc_' FrameMatList{k}(1:(end-4)) '.fig'])
            disp(['saved Htc figure to: Htc_D' num2str(j) '_' FrameMatList{k}(1:(end-4)) '.fig'])
        end
        %Save Domain results to vectors:
        trueMeanA_vec(j,1) = trueMeanA;
        trueMeanD_vec(j,1) = trueMeanD;
        haematocrit_mat(j,:) = haematocrit; 
        haematocrit2D_mat(j,:) = haematocrit2D; 
        Ncells_mat(j,:) = Ncells;
        Hpries_vec(j,1) = Ht(Hf,D_equ)*100;
    end
    if(closeSingleImages)
        close gcf
    end
    if(saveOutputMat)
        HtcResults = struct('Ncells',Ncells_mat,'haematocrit',haematocrit_mat,...
            'haematocrit2D',haematocrit2D_mat,'Area',trueMeanA_vec,'Diameter',...
            trueMeanD_vec,'Hf',Hf*100,'Hpries',Hpries_vec,'channel_Vol',channel_Vol,...
            'channel_length',channel_length);
        
        save([FolderPath 'Htc_' FrameMatList{k}],'HtcResults')
        disp(['saved Htc results to: Htc_' FrameMatList{k} ])
    end
    disp('-------------------------------------')
end
end