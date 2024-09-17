function [X Y Usup Vsup] = ptv2grid_artorg(x,y,u,v,currentimage,res,idx)
%PTV2GRID Lay down the vectors from PTV on a regular grid

%   Copyright, Antoine Patalano Jul 04, 2012
%   antoine.patalano@gmail.com

%Altered by ARTORG Center: CVE-Group
% Author: Damiano Schirinzi
% this function was not really working from the beginning, defy from using
% it.

%   INPUT:
%   x and y data: a cell in the form xydata = {XvaluesOfCoordinates YvaluesOfCoordinates}.
%   u and v: a cell in the form uvdata= {UvaluesOfCVelocity YvaluesOfVelocity}.
%   currentimage ishould be currentimage=imread(filepath{selected})
%   roirect: region of interest
%   res: size of a element. ex:16 (in pixels)
%   OUTPUT:
%   X and Y: a 2D matrix with X and Y coordinates
%   Usup and Vsup: a 2D matrix with U and V velocities

%%
%define the limit of the ROI

    bordmaxx=size(currentimage,2);
    bordmaxy=size(currentimage,1);
    bordminx=0;
    bordminy=0;

%define the grid with the size of the element
xgrid=bordminx:res:bordmaxx;
ygrid=bordminy:res:bordmaxy;


% figure
[X,Y]=meshgrid(xgrid,ygrid);
Usup=nan*X;
Vsup=Usup;
ClustersAsgd=unique(idx);
for i=1:length(ClustersAsgd)
    if length(find(idx==ClustersAsgd(i)))>2 % check if there are more than
        %2 points in one cluster, if not the they won't be interpolated on
        %the grid because you need at least 3 points
        try % for the newest version of MATLAB
        FU=TriScatteredInterp(x(idx==ClustersAsgd(i)),y(idx==ClustersAsgd(i)),u(idx==ClustersAsgd(i)));
        FV=TriScatteredInterp(x(idx==ClustersAsgd(i)),y(idx==ClustersAsgd(i)),v(idx==ClustersAsgd(i)));
        U=FU(X,Y);
        V=FV(X,Y);
        catch    
        U = griddata(x(idx==ClustersAsgd(i)),y(idx==ClustersAsgd(i)),u(idx==ClustersAsgd(i)),X,Y);
        V = griddata(x(idx==ClustersAsgd(i)),y(idx==ClustersAsgd(i)),v(idx==ClustersAsgd(i)),X,Y);
        end
        Usup(isnan(U)==0)=U(isnan(U)==0);% Usup and Vsup are the matrix U and V of each clusters superposed together
        Vsup(isnan(V)==0)=V(isnan(V)==0);
        
    end
end

end