function [A] = findparticleB_AM(particleInfo,s,Ndomains)

% Compared to the standard script I added and "if" statement to check
% whether the x or y coordinate is null. If this statement is true the
% A(i,j) = 0 meaning to say I exclude that particles from the ROI. In the
% other cases the function operates as usual.

A=zeros(Ndomains,length(particleInfo)) ;
for j=1:length(particleInfo)  
    for i=1:Ndomains
        if particleInfo(j,3)==0 || particleInfo(j,2)==0
            A(i,j)=0;
            
        elseif s(i).domains(ceil(particleInfo(j,3)),ceil(particleInfo(j,2)))>0 
                % Adapt accordingly to
                A(i,j)=particleInfo(j,1);                                      % particles main direction
        else
            A(i,j)=0;   
            
        end
    end
end

end