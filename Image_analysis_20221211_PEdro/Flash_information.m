%% tabula rasa
clc
clear all
close all
%Dialogue for choosing a directory:
inputVideoPath = uigetdir('','Select the folder containing the information fo the Flash:');
if(inputVideoPath(end)~='\')
    inputVideoPath = [inputVideoPath '\'];
end
%Get the names of all Falsh files in the specific fokder:
inputVideoNames = dir([inputVideoPath '*_Flash.mat']); 
inputVideoNames = {inputVideoNames.name};

%%
for i=1:1%size(inputVideoNames,2)
    inputVideoName=inputVideoNames{i};
    load([inputVideoPath inputVideoName(1:end-10) '_Pressure_pump_data.mat']);
    load([inputVideoPath inputVideoNames{i}]);
    tstart_microscope(1,i)=Flash(1,1);
    tflash_microscope(1,i)=Flash(1,2);
    frameflash_microscope(1,i)=Flash(1,3);
    tfalsh_on(1,i)=Flash(1,4);
    tfalsh_pump(1,i)=Flash(1,5);
    difference_pump_miroscope(1,i)=tflash_microscope(1,i)-tfalsh_pump(1,i);
    %Flash_summary=[tstart_microscope; tflash_microscope;frameflash_microscope;tfalsh_on,tfalsh_pump;difference_pump_miroscope];
    time_pressurepump(:,i)=Experiment_data(:,1);
    setpressure_pressurepump(:,i)=Experiment_data(:,2);
    measuredpressure_pressurepump(:,i)=Experiment_data(:,3);
    for ii=1:(size(Experiment_data,1))
        difference_pressure_pressure_pump(ii,i)=Experiment_data(ii,2)-Experiment_data(ii,3);
    end
    for ii=1:(size(Experiment_data,1)-1)
        difference_time_pressure_pump(ii,i)=Experiment_data(ii+1,1)-Experiment_data(ii,1);
    end
   %clear Flash Experiment_data Experiment_setup Experiment_setup Experiment_data_info
    
end
maxstart=max(tstart_microscope);
minstart=min(tstart_microscope);
Flash_summary=[tstart_microscope; tflash_microscope; frameflash_microscope;tfalsh_on; tfalsh_pump; difference_pump_miroscope];

%% Plot
 tflash_microscope_name=cellstr(char(ones(size(tflash_microscope,1),1) * 'Time Microscope'));
 tflash_pump_name=cellstr(char(ones(size(tflash_microscope,1),1) * 'Time Pressure Pump'));
 Groups=[tflash_microscope_name; tflash_pump_name];
 Times=[tflash_microscope;tfalsh_pump];
 grouporder={'Time Microscope','Time Pressure Pump'};
 red=[0.6350 0.0780 0.1840];
 green=[0.4660 0.6740 0.1880];
 colors=[red; green];
points1=[1];
points2=[2];
r1=repmat(points1,1,size(tflash_microscope,2))';
r2=repmat(points2,1,size(tflash_microscope,2))';
points=[r1 r2];
figure
%boxplot([tflash_microscope,tfalsh_pump],'Labels',{'Time Microscope','Time Pressure pump'});
%violinplot(Times,  Groups,'Width',0.4,'ViolinColor', colors,'GroupOrder',grouporder);
for i=1:size(tflash_microscope,2)
plot(points(i,:),[tflash_microscope(1,i),tfalsh_pump(1,i)])
hold on
end


%%
inputVideoNames2 = dir([inputVideoPath '*_Pressure_pump_data.mat']); 
inputVideoNames2 = {inputVideoNames2.name};
%%
figure 
hold on
for i=28%:size(inputVideoNames2,2)
    inputVideoName=inputVideoNames2{i};
    load([inputVideoPath2 inputVideoName]);
plot(Experiment_data(:,1),Experiment_data(:,3))
end
%%
figure 
hold on
for i=12:21
    inputVideoName=inputVideoNames2{i};
    load([inputVideoPath inputVideoName]);
plot(Experiment_data(:,1),Experiment_data(:,3))
end

figure 
hold on
for i=22:31
    inputVideoName=inputVideoNames2{i};
    load([inputVideoPath inputVideoName]);
plot(Experiment_data(:,1),Experiment_data(:,3))
end
figure 
hold on
for i=32:41
    inputVideoName=inputVideoNames2{i};
    load([inputVideoPath inputVideoName]);
plot(Experiment_data(:,1),Experiment_data(:,3))
end