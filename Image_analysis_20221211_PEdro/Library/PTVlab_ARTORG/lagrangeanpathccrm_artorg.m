function [dis_result,MAXID]=lagrangeanpathccrm_artorg(dis_result,prev_dis_result,MAXID)

%% use this function to track a particle with it's ID But have to add a
%% new column to dis_result and prev_dis_result . Maybe a 6th coulmn ????

maxID=max(prev_dis_result(:,6));

if maxID<MAXID
    maxID=MAXID;
end
MAXID=maxID;

sizedis=size(dis_result,1);

counter=0;
 for i=1:sizedis

     posrow=find( fix((prev_dis_result(:,3))*1000)/1000==(fix(dis_result(i,1)*1000)/1000));
     poscol=find(fix((prev_dis_result(posrow,4))*1000)/1000==(fix(dis_result(i,2)*1000)/1000));

    
    %%% In case of a new particle a new ID is created
    if length(poscol)==0
        counter=counter+1;
        dis_result(i,6)=maxID+counter;
    end
    
    %%% In case of the same particle the ID is kept
    if length(poscol)~=0
        dis_result(i,6)=prev_dis_result(posrow(poscol(1)),6);
    end
 end

end
    