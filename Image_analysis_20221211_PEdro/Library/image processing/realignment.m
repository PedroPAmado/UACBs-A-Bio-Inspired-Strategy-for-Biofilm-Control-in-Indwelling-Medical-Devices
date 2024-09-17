function [FrameMatReal,ProcessingSteps] = realignment( FrameMat,startframe,endframe,inputVideoPath, inputVideoNames,ProcessingSteps )
% clear all 

A=FrameMat;
AL=length(A);

original =A(:,:,1);
sizemat=size(A);
%% crop for y direction
try 
    load([inputVideoPath inputVideoNames(1:end-4) '_Info_Realigned.mat']);
    exist transformatrix;
    eq(endframe,endframesaved);
    eq(startframe,startframesaved);
catch
    figure
        imshow(A(:,:,1))
        uiwait(msgbox('Select ROI for perycte y direction'));
        [image_realignedy,Realignedy] = imcrop; %Save manually selected ROI in variable ROi
       close all
    aay=size(image_realignedy);
    

    B1y=uint8.empty;
    for ii=1:AL
        B1y(:,:,ii)=imcrop(A(:,:,ii),Realignedy); 
    end

    yfoto1=uint8.empty;
    for ii=1:AL
        yfoto1(:,ii)=B1y(:,aay(2),ii) ;
    end

     BW1y = im2bw(yfoto1,0.5);
  %% Where is the translation    
        figure 
        subplot(2, 1,1);
        imshow(yfoto1)
        axis on
        subplot(2, 1,2);
        imshow( BW1y)
        axis on

        starting= str2double( inputdlg('Wehere is the translation happening ') );
     
        crop1=starting-100;
        crop2=starting+200;
       if startframe>crop1
           crop1=startframe;
       end
       if endframe<crop2
           crop2=endframe;
       end
       
        
    yfotocropped=yfoto1(:,crop1:crop2);
    BW1cropped=BW1y(:,crop1:crop2);
   
        figure 
        subplot(2, 1,1);
        imshow(yfotocropped)
        axis on
        subplot(2, 1,2);
        imshow(BW1cropped)
        axis on

        be= str2double( inputdlg('Insert begin translation') );
        en = str2double( inputdlg('Insert end transletion') );

        beginning= be+crop1;
        ending = en+crop1;
 
        close all
    
    ytransin=zeros(1,AL);
     ytransfin=zeros(1,AL);
    for i=1:AL   
     [row,col]= find(BW1y(:,i)==0);
     ytransin(1,i)=row(1);
     ytransfin(1,i)=ytransin(1,i)-ytransin(1,1);      
    end   

  %% Translation in y
  
  
    transformatrixy=zeros(3,3,AL);
    for ii=1: beginning
      transformatrixy(:,:,ii)=[1 0 0; 0 1 0; 0 0 1];  
    end
    for ii=(beginning+1):ending
       transformatrixy(:,:,ii)=[1 0 0; 0 1 0; 0 -ytransfin(1,ii) 1]; 
     end

    meanendy=floor(mean(ytransfin(1,ending:end)));
    for ii=(ending+1):AL
       transformatrixy(:,:,ii)=[1 0 0; 0 1 0; 0 -meanendy 1]; 
    end
  
    FrameMatRealy = uint8.empty;

    for i=1:AL
    Trans_finy=affine2d(transformatrixy(:,:,i));
    outputViewy = imref2d(size(original));
    FrameMatRealy(:,:,i) = imwarp(A(:,:,i),Trans_finy,'OutputView',outputViewy);
    end
%% crop for x direction  

        figure
        imshow(A(:,:,1))
        uiwait(msgbox('Select ROI for perycte x direction'));
        [image_realignedx,Realignedx] = imcrop; %Save manually selected ROI in variable ROi
        close all
    aax=size(image_realignedx);

    B1x=uint8.empty;
    for ii=1:AL
        B1x(:,:,ii)=imcrop(FrameMatRealy(:,:,ii),Realignedx); 
    end

    xfoto1=uint8.empty;
    for ii=1:AL
     xfoto1(ii,:)=B1x(aax(1),:,ii) ;
    end

     BW1x = im2bw(xfoto1,0.5);
     xtransin=zeros(1,AL);
    xtransfin=zeros(1,AL);
    for i=1:AL   
         [row,col]= find(BW1x(i,:)==0);
         xtransin(1,i)=col(1);
        xtransfin(1,i)=xtransin(1,i)-xtransin(1,1);      
    end 
    %% Translation total
  
  
    transformatrix=zeros(3,3,AL);
    for ii=1: beginning
      transformatrix(:,:,ii)=[1 0 0; 0 1 0; 0 0 1];  
    end
    for ii=(beginning+1):ending
       transformatrix(:,:,ii)=[1 0 0; 0 1 0; -xtransfin(1,ii) -ytransfin(1,ii)  1]; 
    end
    meanendx=floor(mean(xtransfin(1,ending:end)));

    for ii=(ending+1):AL
       transformatrix(:,:,ii)=[1 0 0; 0 1 0; -meanendx -meanendy  1]; 
    end
  
   
   
end
 FrameMatReal =uint8.empty;
for i=1:AL
    Trans_fin=affine2d(transformatrix(:,:,i));
    outputView = imref2d(size(original));
    FrameMatReal(:,:,i) = imwarp(A(:,:,i),Trans_fin,'OutputView',outputView);
end
  
%%  figure for x direction 

B2x=uint8.empty;
for ii=1:AL
    B2x(:,:,ii)=imcrop(FrameMatReal(:,:,ii),Realignedx); 
end
aax=size(B2x);
xfoto2=uint8.empty;
for ii=1:AL
   xfoto2(ii,:)=B2x(aax(1),:,ii);
end

beginningfin=beginning-500;
endingfin=ending+500;

 if startframe>beginningfin
           beginningfin=startframe;
 end
 if endframe<endingfin
           endingfin=endframe;
 end
    
xfoto11=xfoto1((beginningfin):(endingfin),:);

xfoto22=xfoto2((beginningfin):(endingfin),:);

figure 
subplot(1,2,1);
    imshow(xfoto11)
subplot(1, 2,2);
  imshow(xfoto22)
 savefig([inputVideoPath inputVideoNames(1:end-4) '_Before_after_x_direction.fig']) ;
close all
%%  figure for y direction 
B2y=uint8.empty;

for ii=1:AL
    B2y(:,:,ii)=imcrop(FrameMatReal(:,:,ii),Realignedy); 
end
aay=size(B2y);
yfoto2=uint8.empty;
for ii=1:AL
   yfoto2(:,ii)=B2y(:,aay(2),ii);
end

yfoto11=yfoto1(:,(beginning-100):(ending+100));

yfoto22=yfoto2(:,(beginning-100):(ending+100));

figure 
subplot(2, 1,1);
    imshow(yfoto11)
subplot(2, 1,2);
  imshow(yfoto22)
 savefig([inputVideoPath inputVideoNames(1:end-4) '_Before_after_y_direction.fig']) ;
close all

endframesaved=endframe; %save the endframe
startframesaved=startframe; %save the startframe

%% saved realign
save([inputVideoPath inputVideoNames(1:end-4) '_Info_Realigned.mat'],'transformatrix','FrameMatReal','beginning','ending','yfoto1','yfoto2','xfoto1','xfoto2','Realignedy','Realignedx','endframesaved','startframesaved');
mat2avi(FrameMatReal,106,[inputVideoPath inputVideoNames(1:end-4) '_realigned.avi'],false);
ProcessingSteps.realignment = true;
end