function Image_Angle=Rotate(framefigure)
imshow(framefigure)
uiwait(msgbox('Select ROI for the rotation.'));
[framefigure,ROI] = imcrop;
Image= framefigure;
ImgTL = Image; 

ImgInv=ImgTL;%imcomplement(ImgTL);
ImgEdge = edge(ImgInv,'Sobel',[],'both','nothinning');
ImgMed = medfilt2(ImgEdge);


[H,T,R] = hough(ImgMed,'RhoResolution',1 );
%Create a figure (subplot) with two stacked plots
%one for images one for hough space plot
figure
subplot(2,1,1); 
imshowpair(ImgMed, ImgEdge, 'montage');
title('H.png');
subplot(2,1,2);

imshow(imadjust(rescale(H)),'XData',T,'YData',R,...
     'InitialMagnification','fit');
title('Hough transform of wells.jpg');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(gca,hot);

%Extract hough peaks, from the peaks find lines and plot them on the image.
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));

lines = houghlines(ImgMed,T,R,P,'FillGap',5,'MinLength',100);

%Rotate image based on average theta value
      lines_t=struct2table(lines); 
      lines_a=table2array(lines_t);
      for ii=1:size(lines_a,1)
          if lines_a(ii,5)<0
              lines_a(ii,5)=90+lines_a(ii,5);
          else
              lines_a(ii,5)=-90+lines_a(ii,5);
          end
      end
      
      Image_Angle=mean(lines_a(:,5)); %extract theta value from array

      
      theta=5;
      
    if abs(Image_Angle)<theta %only roate if angle is smaller than theta
        
    else 

            Image_Angle=0;

    end

 close all 

