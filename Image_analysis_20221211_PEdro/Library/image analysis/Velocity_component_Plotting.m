function [Stats] = Velocity_component_Plotting(inputVideoPath,inputVideoNames,UserSettingsGUI,f)
%Authors: ARTORG-Center, CVE-Group;
% - Alberto Mantegazza
% - Damiano Schirinzi

% This script is an adaptation of OptC_AM and finalD_AM
% When can I use it? You can use this script if you have to analyze a single channel and you want to compute
% and display the two components of particles velocity. I could have added some lines in the aforementioned
% codes but would be confusionary because every time you wanted to use it you had to enable/disable some lines.
% Therefore, the idea to create a different code.

% Requirements: You need to run first the following scripts: roiDefA, findparticleB, countD.
% ATTENTION: CountD recalls the function CountC and here you must enable 2 lines (storing of Vx and Vy)
% if you want to be able this code.

%% Plot for each domain:
load([inputVideoPath 'PTV_' inputVideoNames(1:(end-4)) '.mat']);
Ndomains = size(f,2);
UserSettingsGUIstring = UserSettingsGUI{1};
UserSettingsGUIvalue = UserSettingsGUI{2};
Stats ={};
%Calibration: retreive frame rate and microscope calibration   
space_cal=UserSettingsGUIvalue.calibrationEdit; % calibration: [�m/pixel]
fps=UserSettingsGUIvalue.framerateEdit; %[fps]
time_cal=1/fps; % Time calibration: 1 frame = 1/298 s
%Calculate domain lengths for Line Density calculations
nNormalDomains = size(UserSettingsGUIvalue.domainROI,1); %since Ndomains = nNormalDomains+nGAdomains
for j=1:Ndomains
    if(j<=nNormalDomains) %Domain is rectangular
        length_domain(j)=UserSettingsGUIvalue.domainROI(j,3)*space_cal; % [�m]
    else %Domain is a Polygon (manually drawn rectangle)
        %Euclidian distance
        %verteces have to be the verteces of a rectangle more or less!
        vertices = UserSettingsGUIvalue.GAvertices{j-nNormalDomains}; %Format: [x1 y1 x2 y2..]
        la = sqrt((vertices(1)-vertices(3))^2+(vertices(2)-vertices(4))^2);
        lb = sqrt((vertices(5)-vertices(3))^2+(vertices(6)-vertices(4))^2); 
        length_domain(j)= max(la,lb)*space_cal; %longer side of the rectangle [�m]
    end
end

for j=1:Ndomains
    %% Velocity Matrix and RBC Matrix creation
    
    C1=f(j).domains(2:4,:); % It is a velocity matrix. 1st row: V; 2nd row: Vx; 3rd row: Vy.
    %old piece of code deleting entries:
    %D1=find(C1(1,:)>0); % Find when the velocity is different than zero (or NaN)
    %new code line:
    if(any(C1(1,:)<=0))
        disp('Warning: Velocity_component_Plotting.m; zero or negative absolute velocities found!')
    end
    C1(:,C1(1,:)<=0) = NaN; %Convert zero velocity and negative velocities to NaN
    nanInd = isnan(C1(1,:));
    if(any(isnan(C1(1,:)))) %If any NaN values in absolute velocity then interpolate
        %Interpolate around NaN values (linar interpolation from point before to point after)
        %[C1,indexInterpolated] = fillmissing(C1,'linear','3'); %only works in r2016b and newer
        xx = 1:size(C1,2);
        ss = [];
        % xx(~nanInd) contains the indices of the non-NaN values
        % C1(~nanInd) contains the non-NaN values
        % xx(nanInd) contains the indices of the NaN values, and thus the points at
        % which we would like to query the spline interpolator
        ss(1,:) = interp1(xx(~nanInd),C1(1,~nanInd),xx(nanInd));
        ss(2,:) = interp1(xx(~nanInd),C1(2,~nanInd),xx(nanInd));
        ss(3,:) = interp1(xx(~nanInd),C1(3,~nanInd),xx(nanInd));
        %replace NaN values with interpolated values:
        C1(:,nanInd) = ss;
    end
    %old piece of code deleting entries:
    %D2=C1(1:3,D1); % Saving V, Vx, Vy (Row-wise)
    
    D2=C1;
    
    RBC=f(j).domains(1,:);
    RBC(isnan(RBC)) =0;
    if(sum(nanInd)~=sum(RBC==0))
        disp('Warning:Velocity_component_Plotting.m; possible mismatch in U and RBC NaN entries!')
        %Mismatch given as we set to NaN also zero velocity entries--< but
        %velocity zero is an error of the PTV, cells always move in our
        %experiments!-> Set RBC entries also to zero.
        RBC(nanInd)=0;
    end
    %old piece of code deleting entries:
    %S2=RBC(1,D1);
    S2=RBC;
    
    %% Signal Filtering
    filtwindow = min(size(S2,2),199);
    filtwindow = filtwindow-((mod(filtwindow,2)+1)*1); %filtwindow must be odd
    
    RBCn1f = sgolayfilt(S2(1,1:end),2,filtwindow);
    U1f = sgolayfilt(D2(1,1:end),2,filtwindow);
    U1_xf = sgolayfilt(D2(2,1:end),2,filtwindow);
    U1_yf = sgolayfilt(D2(3,1:end),2,filtwindow);
    %% Finde were error is big
  Lineerr=zeros(1,length(error_per));
    for i=1:length(error_per)
        if error_per(i)>0.5
            Lineerr(i)=i;
        else 
            Lineerr(i)=0;
        end
    end
%Linefin = nonzeros(Lineerr);
    
    
Linencell=zeros(1,length(Ncells));
meanncells=mean(Ncells);
 for i=1:length(Ncells)
      if abs((Ncells(1,i)-meanncells)/meanncells)>0.5
            Linencell(i)=i;           
      else 
            Linencell(i)=0;
      end
 end   

    Linefin = nonzeros(Linencell);
    %% read begin and end
    try
        infoVideo=load ([inputVideoPath 'Info_Realigned.mat' ]);
        sx=size(infoVideo);
        for a=1:sx(1)
            if char(inputVideoNames)==infoVideo(a).Video
                ii=a;
            end
        end   
        beginrel=infoVideo(ii).beginning;
        endrel=infoVideo(ii).ending;
        MeanbeginU1f=mean(U1f(1,1:beginrel));
        MeanbeginU1_xf=mean(U1_xf(1,1:beginrel));
        MeanbeginU1_yf=mean(U1_yf(1,1:beginrel));
        x=1000;
        MeanendU1f=mean(U1f(1,(endrel+x):end));
        MeanendU1_xf=mean(U1_xf(1,(endrel+x):end));
        MeanendU1_yf=mean(U1_yf(1,(endrel+x):end));
        %% Save info to plot
        save([inputVideoPath inputVideoNames(1:(end-4)) '_InfoPlot_Domain_' num2str(j) '.mat'],'RBCn1f',...
          'Lineerr','Linencell','Linefin','beginrel','endrel','U1f','U1_xf','U1_yf','MeanbeginU1f',...
          'MeanbeginU1_xf','MeanbeginU1_yf','MeanendU1f','MeanendU1_xf', 'MeanendU1_yf');
        infovideo=true;
    catch 
        infovideo=false;
         save([inputVideoPath inputVideoNames(1:(end-4)) '_InfoPlot_Domain_' num2str(j) '.mat'],'RBCn1f',...
          'Lineerr','Linencell','Linefin','U1f','U1_xf','U1_yf');
    end    
        
 
    
    %% Plotting Number of RBC\Mean Velocity\Mean Velocity Components
    
    figure
    subplot(3,1,1);
    plot(1:length(RBCn1f), RBCn1f(1,:),'r','LineWidth',1); grid on;
    for Li=1:length(Linefin)
       xline(Linefin(Li),'-','k');
    end
    if infovideo
        xline(beginrel,'--','m');xline(endrel,'--','g');
    end
    xlabel('Frames'); ylabel('RBCn - Sgolayfilt -'); legend('RBCnf');
    title('Number of RBC (t)','FontSize',14)
    
    subplot(3,1,2);
    plot(1:length(U1f), U1f(1,:),'r','LineWidth',2);
    hold on
    for Li=1:length(Linefin)
        xline(Linefin(Li),'-','k');
    end
    if infovideo
        plot(1:beginrel, MeanbeginU1f*ones(1,beginrel),'k');
        plot((endrel+x):length(U1f),MeanendU1f *ones(1,(length(U1f)-(endrel+x)+1)),'k');
        xline(beginrel,'--','m');
        xline(endrel,'--','g');
    end
    xlabel('Frames');
    ylabel('Velocity [ Pixel / Frame ] - Sgolayfilt -');
    legend('U'); 
    grid on
    title('Mean Velocity (t)','FontSize',14)
    hold off
    
    subplot(3,1,3);
    plot(1:length(U1_xf), U1_xf(1,:),'r',1:length(U1_yf), U1_yf(1,:),'b','LineWidth',1);
    hold on
    for Li=1:length(Linefin)
        xline(Linefin(Li),'-','k');
    end
    if infovideo
        plot(1:beginrel, MeanbeginU1_xf*ones(1,beginrel),'k');
        plot(1:beginrel, MeanbeginU1_yf*ones(1,beginrel),'k');
        plot((endrel+x):length(U1_xf),MeanendU1_xf *ones(1,(length(U1_xf)-(endrel+x)+1)),'k');
        plot((endrel+x):length(U1_yf),MeanendU1_yf *ones(1,(length(U1_yf)-(endrel+x)+1)),'k');
        xline(beginrel,'--','m');
        xline(endrel,'--','g');
    end
    xlabel('Frames'); 
    ylabel('Velocity [ Pixel / Frame ] - Sgolayfilt -'); 
    legend('Ux','Uy');
    grid on
    title('Mean Velocity - V_x(t) and V_y(t)','FontSize',14)
    hold off
    suptitle(['Video: ' inputVideoNames '; Domain ' num2str(j)])
    savefig([inputVideoPath inputVideoNames(1:(end-4)) '_Domain_' num2str(j) '.fig'])
       
    
    %% Plotting    
    U1=((U1f.*space_cal)./time_cal)./1000; %[mm/s]
    
%     length_domain=UserSettingsGUIvalue.domainROI(j,3)*space_cal; or  length_domain(j)= max(la,lb)*space_cal % [�m]
    Line_den1=(RBCn1f(1,:)./length_domain(j))*1000; %[NRBC/mm]
    
%     RBC_flux1=U1(:).*Line_den1(:);
    RBC_flux1=U1.*Line_den1;
    
    Stats{j}= struct('U',U1,'Line_den',Line_den1,'RBC_flux',RBC_flux1,'RBCnf',RBCn1f);
    
    
end