function [particleInfo] = PTVlab_artorg(FrameMat, UserSettings, ProcessingSteps, PTVsettings)

    resultslistptv=cell(0);
    resultslist=cell(0);%clear old results
    MAXID = 0; %Particle ID
    particleInfo = []; %Output with all particles
   
    h = waitbar(0, {'PTV analysis: 00 %', 'Time left: N/A' });
    if(PTVsettings.visualPTVControl)
        %initialize visual check
        figCheck = figure;
        stoppedByUser = false;
        hold on
        %toggle button for stopping of video
        stopButton = uicontrol(figCheck, 'Style', 'togglebutton', 'String', 'Close');
        positionStop = get(stopButton,'Position');
        positionPause = [(positionStop(1)+positionStop(3)+10) positionStop(2:4)];
        positionPlay = [(positionPause(1)+positionPause(3)+10) positionStop(2:4)];
        %toggle button for pausing of video
        pauseButton = uicontrol(figCheck, 'Style', 'togglebutton', 'String', 'Pause','Position',positionPause );
        %toggle button for resuming of video
        playButton = uicontrol(figCheck, 'Style', 'togglebutton', 'String', 'Play','Position',positionPlay);
    end
    
    
    for i=1:2:(size(FrameMat,3)*2-2) %Loops through every second frame
        if i==1
           tic
        end
            image1= FrameMat(:,:,(i+1)/2); %/Replace mit FrameMat()
            image2= FrameMat(:,:,(i+1)/2+1); %/Replace mit FrameMat()
            if size(image1,3)>1
                image1=uint8(mean(image1,3));
                image2=uint8(mean(image2,3));
                disp('Warning: To optimize speed, your images are being set to grayscale, 8 bit!')
            end

            switch PTVsettings.ParticleDetectionAlgorithm
                case 'Gussian Mask'
                    gaussdetecmark= true;
                    dynadetecmark= false;
                case 'Dynamic Threshold'
                    gaussdetecmark= false;
                    dynadetecmark= true;
            end
            
            corrthreval= PTVsettings.corrthreval;
            sigmasize=  PTVsettings.sigmasize;
            intthreval=  PTVsettings.intthreval;
            
            % Run the particle detection first
            [image1, row1, col1]= PTVlab_detection_artorg (image1,gaussdetecmark,corrthreval,sigmasize,intthreval,dynadetecmark);
            [image2, row2, col2]= PTVlab_detection_artorg (image2,gaussdetecmark,corrthreval,sigmasize,intthreval,dynadetecmark);
                                       
            if size(resultslistptv,2)<(i+1)/2-1 || (i+1)/2==1 %Check if the results of the previous pair has been calculated
                indicator=1;
                prev_dis_result=[];
            else
                %             isempty(resultslistptv{1,(i+1)/2-1})
                if isempty(resultslistptv{1,(i+1)/2-1})==1
                    indicator=1;
                    prev_dis_result=[];
                elseif isempty(resultslistptv{1,(i+1)/2-1})==0
                    indicator=0;
                    prev_dis_result=[];
                    prev_dis_result(:,1)=resultslistptv{1,((i+1)/2-1)};
                    prev_dis_result(:,2)=resultslistptv{2,((i+1)/2-1)};
                    prev_dis_result(:,3)=resultslistptv{3,((i+1)/2-1)};
                    prev_dis_result(:,4)=resultslistptv{4,((i+1)/2-1)};
                    prev_dis_result(:,5)=resultslistptv{5,((i+1)/2-1)};
                    prev_dis_result(:,6)=resultslistptv{6,((i+1)/2-1)};
                end
            end
            
            switch PTVsettings.PTValgorithm 
                case 'Cross-correlation (cc)'
                    ccmark=true;%/ Cross correlation method
                    rmmark=false;%/ relaxation method 
                    hymark=false;%/ hybrid method 
                case 'Relaxation'
                    ccmark=false;%/ Cross correlation method
                    rmmark=true;%/ relaxation method 
                    hymark=false;%/ hybrid method 
                case 'Hybrid'
                    ccmark=false;%/ Cross correlation method
                    rmmark=false;%/ relaxation method 
                    hymark=true;%/ hybrid method mark 
            end
            
            switch PTVsettings.Crosscorrlparameters
               case 'By interrogation area'
                    det_nummark=false;%/ Relaxation with n particles ...%get(handles.det_num,'value');
                    det_areamark=true;%/ Relaxation with detection Area ...%get(handles.det_area,'value');
                case 'By min. # of particles'
                    det_nummark=true;%/ Relaxation with n particles ...%get(handles.det_num,'value');
                    det_areamark=false;%/ Relaxation with detection Area ...%get(handles.det_area,'value');
                
            end
    

            num_part = PTVsettings.minParticles; %/Relaxation Model nParticles
            area_size = PTVsettings.interrArea; %/ lenght area
            corrcc = PTVsettings.minCorr;
            percentcc = PTVsettings.simNeigh;
            tn=37; %str2double(get(handles.tn,'string'));
            tq=1 ;%str2double(get(handles.tq,'string'));
            minneifrm=1; %get(handles.minneifrm,'string');
            tqfrm1=88; %str2double(get(handles.tqfrm1,'string'));
            minprob=0; %str2double(get(handles.minprob,'string'));
            tqfcc=80;
            epsilon=0.01;
            percentrm=70;
            ninit=1;
            nframe=(i+1)/2;
            
            %run the ptv algorithm and save the results in resultslistptv
            [dis_result,indicator,MAXID]=ptv_CCRM_artorg(image1,image2,num_part,tn,tq,tqfcc,tqfrm1,percentcc,percentrm,epsilon,...
                corrcc,minprob,ccmark,rmmark,hymark,minneifrm,indicator,det_nummark,det_areamark,area_size,MAXID,...
                row1,col1,row2,col2,prev_dis_result,ninit,nframe);
            
            if ~isempty(dis_result)
%                 if isempty(roirect)==1
                    resultslistptv{1,(i+1)/2}=dis_result(:,1);
                    resultslistptv{2,(i+1)/2}=dis_result(:,2);
                    resultslistptv{3,(i+1)/2}=dis_result(:,3);
                    resultslistptv{4,(i+1)/2}=dis_result(:,4);
                    resultslistptv{5,(i+1)/2}=dis_result(:,5);
                    resultslistptv{6,(i+1)/2}=dis_result(:,6);
    
                try
                    x=resultslistptv{2,(i+1)/2};
                    y=resultslistptv{1,(i+1)/2};
                    typevector=resultslistptv{5,(i+1)/2};
                    u=resultslistptv{4,(i+1)/2}-resultslistptv{2,(i+1)/2};
                    v=resultslistptv{3,(i+1)/2}-resultslistptv{1,(i+1)/2};
                    id=resultslistptv{6,(i+1)/2};
                    typevector=resultslistptv{5,(i+1)/2};
                    
                    % Position and Velocities
                    particleInfo = [particleInfo; nframe*ones(size(x,1),1) x y u v id];
                                       
                    %make cluster of points. idx is the index of each cluster
                    RadiusCluster=80; %in pixel
                    idx=ncluster(x(typevector==1),y(typevector==1),RadiusCluster);
                    
                    %Give the matrix X, Y, U  V (can be improved)
                    meshsize=10; %(in pixel)
                    currentimage=imread(FrameMat(:,:,(i+1)/2));
                    [X, Y, U, V] = ptv2grid_artorg(x,y,u,v,currentimage,meshsize,idx);
                    
                    %save it in resultlist
                    resultslist{1,(i+1)/2}=X;
                    resultslist{2,(i+1)/2}=Y;
                    resultslist{3,(i+1)/2}=U;
                    resultslist{4,(i+1)/2}=V;
                catch
                    resultslist{1,(i+1)/2}=[];
                    resultslist{2,(i+1)/2}=[];
                    resultslist{3,(i+1)/2}=[];
                    resultslist{4,(i+1)/2}=[];
                end
                
                
            else
                resultslistptv{1,(i+1)/2}=nan;
                resultslistptv{2,(i+1)/2}=nan;
                resultslistptv{3,(i+1)/2}=nan;
                resultslistptv{4,(i+1)/2}=nan;
                resultslistptv{5,(i+1)/2}=nan;
                resultslistptv{6,(i+1)/2}=nan; %Added while debugging
                
                resultslist{1,(i+1)/2}=nan;
                resultslist{2,(i+1)/2}=nan;
                resultslist{3,(i+1)/2}=nan;
                resultslist{4,(i+1)/2}=nan;
                resultslist{5,(i+1)/2}=nan;
                resultslist{6,(i+1)/2}=nan; %Added while debugging
            end
            
            %visual check plot
            if(PTVsettings.visualPTVControl && stoppedByUser~=true)
                %Get centroids from PTV results
                indexFrame = find(particleInfo(:,1)==nframe);
                centroids = particleInfo(indexFrame,2:3);
                radii = 2*ones(1,size(centroids,1)); %[px]
                %Plot to check bloob detection with Centroids found
                %plot frame
                imshow(FrameMat(:,:,nframe),[])
                %plot centroids
                viscircles(centroids,radii);
                drawnow
                if (get(stopButton, 'Value')==1)
                    stoppedByUser = true;
                    close(figCheck)
                elseif (get(pauseButton, 'Value')==1)
                    waitfor(playButton,'Value',1)
                    set(playButton,'Value',0)
                    set(pauseButton,'Value',0)
                else
                    %continue to next frame
                end
            else
                %continue without plotting
            end
            
            zeit=toc;
            done=(i+1)/2;
            tocome=((size(FrameMat,3)-1))-done;
            zeit=zeit/done*tocome;
            hrs=zeit/60^2;
            mins=(hrs-floor(hrs))*60;
            secs=(mins-floor(mins))*60;
            hrs=floor(hrs);
            mins=floor(mins);
            secs=floor(secs);
            waitbar((i+1)/2/size(FrameMat,3),h,{sprintf('PTV analysis: %2.0f %%', (i+1)/2/size(FrameMat,3)*100), ['Time left: ' sprintf('%2.2d', hrs) 'h ' sprintf('%2.2d', mins) 'm ' sprintf('%2.2d', secs) 's'] });
    end
    waitbar((i+1)/2/size(FrameMat,3),h,{sprintf('PTV analysis: %2.0f %%', (i+1)/2/size(FrameMat,3)*100), 'Time left: N/A' });
    delete(h)
end
