function [b, uBulk, Kt, rfitstats ] = intrinsicViscosity(inputVideoPath, inputVideoName, UserSettingsGUI,f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

UserSettingsGUIstring = UserSettingsGUI{1};
UserSettingsGUIvalue = UserSettingsGUI{2};
%Retreive number of analysis domain (exclude general analysis domains)
domainMask = UserSettingsGUIvalue.domainMask;
Ndomains = size(domainMask,3);
%Retreive Usersettigns for unit transformation
calibration=UserSettingsGUIvalue.calibrationEdit; % calibration: [µm/pixel]
fps=UserSettingsGUIvalue.framerateEdit; %[fps]
%Load in data from PTV and haematocrit calculations
load([inputVideoPath 'Htc_' inputVideoName(1:(end-4)) '.mat'])

b_mat = zeros(2,Ndomains);
Kt_cellarray = {};
uBulk_cellarray = {};

for j=1:Ndomains
disp(['Domain ' num2str(j)]); 
%% Calculation of bulk velocity u0
%Try to find frames were there is just one cell in the domain and get the
%mean velocity of those single cells as estimate of the velocity of the
%plasma.
index1cell = find(f(j).domains(1,:)==1);
u1 = f(j).domains(2,index1cell); %[px/frame]

%Do a linear regression
%Take a vector with number of cells in domain
%Another with the average RBC velocity as funciton of the number of cells
nCells = 1:(max(f(j).domains(1,:)));
numberMesPoints = zeros(length(nCells),1);
for i=1:(max(f(j).domains(1,:)))
    index1cell = find(f(j).domains(1,:)==i);
    numberMesPoints(i)=length(index1cell);
    u(i) = mean(f(j).domains(2,index1cell));
end

%Remove Nan entries with indexNan(That number of cells never happened in the domain)
%Remove that is not statistically relevant with indexLow
indexNan = find(isnan(u));
indexLow = find(numberMesPoints<200);
indexGood = find(numberMesPoints>=200);
if(length(indexLow) >= (length(numberMesPoints)-2))
    disp('Warning Linear Regression: not sufficient datapoints reached threshold of 200 frames.')
    disp('Did the regression on all available data:')
    index0 = find(numberMesPoints==0);
    nCells(index0) = [];
    u(index0) = [];
    numberMesPoints(index0) = [];
    T = table(nCells',numberMesPoints,'VariableNames',{'Cells','frames'})    
else
    nCells(indexLow) = [];
    u(indexLow) = [];
    disp('Regression Data:')
    T = table(nCells',numberMesPoints(indexGood),'VariableNames',{'Cells','frames'})
end


u = u*calibration*fps/1000; %[mm/s]

%Linear Regression to find intercept with y-Axis i.e. velocity in cell free
%channel u0:
NCells = [ones(length(nCells),1) nCells'];
b = NCells\u';
%Robust linear regression
[brob,rfitstats(j)] = robustfit(nCells,u);
%Output u0
u0 = brob(1); %[mm/s]
b_mat(:,j) = brob; %[mm/s]

%Plot linear regression:
figure
scatter(nCells,u)
hold on
yCalc2 = [1 0; NCells]*b;
yCalcRob = [1 0; NCells]*brob;
plot([0 nCells],yCalc2,'--','MarkerSize',6,'LineWidth',2)
plot([0 nCells],yCalcRob,'-.','MarkerSize',6,'LineWidth',2)
title({['Video: ' inputVideoName ' ; Domain ' num2str(j)],'Linear Regression u = m*Ncells+u_0'},'FontSize',14)
legend({'Data',['l.r.: (m= ' num2str(b(2),'%1.3e') '; u_0= ' num2str(b(1),'%2.3f') ')'],...
    ['r.l.r.: (m= ' num2str(brob(2),'%1.3e') '; u_0= ' num2str(brob(1),'%2.3f') ')']},'FontSize',14,'Location','best')
minU = min(min(u)*0.9,brob(1));
maxU = max(max(u)*1.1,brob(1));
ylim([minU maxU])
xlabel('Ncells','FontSize',14)
ylabel('u [mm/s]','FontSize',14)


%Calculate bulk velocity uBulk (velocity of the whole suspension) over the 
%haematocrit and RBC velocity following the formula:
% H_tube/H_discharge = uBulk/uRBC -> solve for uBulk 
%Do not consider frames where there are 0 particles
index=find(f(j).domains(1,:)>0);
Hfeed=UserSettingsGUIvalue.feedHaemEdit/100; %Fraction not in [%]
uRBC = f(j).domains(2,index); %[px/frame]
uRBCx = f(j).domains(3,index);%[px/frame]
haematocrit = HtcResults.haematocrit(1,index)./100; %Fraction not in [%]
uBulk = haematocrit.*uRBC./Hfeed; %[px/frame]
%unit transformation
uBulk = uBulk*calibration*fps/1000; %[mm/s]
uRBC = uRBC*calibration*fps/1000; %[mm/s]
uRBCx = uRBCx*calibration*fps/1000; %[mm/s]
%filter the velocity with Savitzky-Golay filter
filtF = round(length(uBulk)*2/80); %Weighting length of sgolay filter
filtF = filtF+(~mod(filtF,2)); %Ensuring weighting length is odd
uBulkf = sgolayfilt(uBulk,2,filtF);

uBulk_cellarray{j} = uBulk;
%% Calculation of intrinsic viscosity
%Formula:
% Kt = 1/H_tube*[u0/uRBC-1]
Kt = ((u0./uRBC) -1)./haematocrit;

Kt_cellarray{j}=Kt;
%% Plot
figure
plot(1:length(uBulk),uBulk)
hold on
plot(1:length(uBulkf),uBulkf)
title({['Video' inputVideoName '; Domain ' num2str(j)],'Bulk velocity u_{bulk} = u_{RBC}*H_{t}/H_{f}'})
xlabel('frames containing cells')
ylabel('u_{bulk} [mm/s]')
figure
plot(1:length(Kt),Kt)
title({['Video' inputVideoName '; Domain ' num2str(j)],'Intrinsic Viscosity K_t'})
xlabel('frames containing cells')
ylabel('K_t [-]')


end



end