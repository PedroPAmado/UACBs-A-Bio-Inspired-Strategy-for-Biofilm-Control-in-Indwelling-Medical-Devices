function [ ] = generalAnalysis( inputVideoPath, inputVideoNames, UserSettingsGUI)
%General Analysis:
%@brief: calculates average RBC velocity / RBC line density / RBC flux


%% roiDefA - Calculation of s struct from domainMask
UserSettingsGUIvalue = UserSettingsGUI{2};
UserSettingsGUIstring = UserSettingsGUI{1};

domainMask = UserSettingsGUIvalue.domainMask;
Ndomains = UserSettingsGUIvalue.nDomainsEdit;
for i=1:Ndomains
    BW=domainMask(:,:,i);
    value(1,i) = {BW};
end
GAmask = UserSettingsGUIvalue.GAdomains;
nGAdomains = UserSettingsGUIvalue.GAdomainsEdit;
for i=1:nGAdomains
    BW=GAmask(:,:,i);
    value(1,(i+Ndomains)) = {BW};
end
s = struct('domains',value);

%% findparticleB_AM & countD - generate the matrix A and the struct f from particleInfo

for i=1:length(inputVideoNames)
        disp(['General Analysis of: ' inputVideoNames{i}(1:(end-4))])
        %Load PTV mat file from corresponding input video
        try
            load([inputVideoPath 'PTV_' inputVideoNames{i}(1:(end-4)) '.mat'])
        catch
            disp(['Error: Could not load "' inputVideoPath 'PTV_' inputVideoNames{k}(1:(end-4)) '.mat"'])
            return
        end
        %script generating Matrix A
        A = findparticleB_AM(particleInfo,s,(Ndomains+nGAdomains));
        %CountParticles, generate struct f
        f =countD(A, particleInfo,(Ndomains+nGAdomains));
        %Velocity_component_Plotting: filters velocities, plots it and
        %calculates RBC velocity / RBC line density / RBC flux for all
        %domains --> are saved in struct Stasts of size(1xNdomains),
        %containing the files U,Line_den,RBC_flux for each domain.
        %Acces domain stats via Stats{i}.RBC_flux ect..
        [Stats] = Velocity_component_Plotting(inputVideoPath,inputVideoNames{i},UserSettingsGUI,f);
        %Intrinsic viscosity and bulk velocity calculations        
        if(Ndomains>0) %Intrinisc vicosity calculation only possible with haematocrit
            [b,uBulk,Kt,rfitstats] = intrinsicViscosity(inputVideoPath,inputVideoNames{i},UserSettingsGUI,f);
        end
        %Save calculations to matrix file
        try
            save([ inputVideoPath 'PTV_' inputVideoNames{i}(1:(end-4)) '.mat'],...
                'particleInfo','PTVsettings','error','Ncells','error_per','CRID','s','f','A','Stats','u','st_u','b','uBulk','Kt')
        catch
            save([ inputVideoPath 'PTV_' inputVideoNames{i}(1:(end-4)) '.mat'],...
                'particleInfo','PTVsettings','error','error_per','Ncells','CRID','s','f','A','Stats','b','uBulk','Kt','rfitstats')
        end
        disp(['Added s,A,f,Stats,uBulk,b,rfitstats and Kt to PTV results in: ' 'PTV_' inputVideoNames{i}(1:(end-4)) '.mat'])
        disp('-------------------------------------------------------')
                
end

end