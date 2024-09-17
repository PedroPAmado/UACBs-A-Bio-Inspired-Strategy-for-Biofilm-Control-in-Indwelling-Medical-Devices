
function [f]=countD(A, particleInfo,Ndomains)
 
%the output f is as struct which has as output on rows the number of
%domains and inside of each domains there is a matrix with in the first row
% the number of partcicles per frame per domain and in the second row the
% mean velocity.

field = 'domains';

for i=1:Ndomains

    for j=1:max(A(i,:))
           D=find(A(i,:)==j);
           A1(1,j)=length(D);
           A1(2,j)=mean(particleInfo(D,7)); % Mean velocity of all particles spotted in the frame j

           % Extra code lines: switch on/off depending on the need
           A1(3,j)=mean(particleInfo(D,4)); % Mean Vx
           A1(4,j)=mean(particleInfo(D,5)); % Mean Vy
           clear D
    end
   value(1,i) = {A1};

end
 f= struct(field,value);
end

    
           
       

    