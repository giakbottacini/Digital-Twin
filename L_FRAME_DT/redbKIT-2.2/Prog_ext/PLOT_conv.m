%% PLOT graph for POMs convergence
clc
clear all
close all
pb_ID_f = '0_2';
start=['/Users/matteotorzoni/Desktop/Corigliano/Mensola_matteo/Dati/mensola_calore/Danneggiamento/damaged_',pb_ID_f];
%Prompt user for filename
[fname, pname] = uigetfile('*.mat','Choose input file',start);  
%Create fully-formed filename as a string
filename = fullfile(pname, fname);
%Read in the data
POM_conv= load(filename);
X=POM_conv.node_positions_conv(1,:);
X=[0 X];
time=0.001:0.001:0.4;
dofs=POM_conv.dofID;
%Estraggo i POM dalla struttura
V1=POM_conv.V1(:,1);
V2=POM_conv.V1(:,2);

%Ordino i POM
V1_R = reshape(V1, 494, 400);
V2_R = reshape(V2, 494, 400);
for j=1:size(V1_R,2)
    V1_p(1,j)=0;
    V2_p(1,j)=0;
    for i=1:length(dofs)
        V1_p(i+1,j)=V1_R(dofs(i),j);
        V2_p(i+1,j)=V2_R(dofs(i),j);
    end
end

[~, ~, ~] = mkdir([start,'/PLOT_conv']);

UNO=figure;
hold on
for j=1:size(V1_p,2)
    if V1_p(end,j)>=0
        plot3(X,time(j)*ones(size(X)),V1_p(:,j))
    else
        plot3(X,time(j)*ones(size(X)),-V1_p(:,j))
    end
end
set(gca,'fontname','times')
xlabel('x [m]','FontSize',16); ylabel('Training time [sec]','FontSize',16); zlabel('Amplitude [-]','FontSize',16)
box on
view(3)
saveas(UNO,fullfile([start,'/PLOT_conv'],'POM_1'),'fig')
saveas(UNO,fullfile([start,'/PLOT_conv'],'POM_1'),'pdf')

DUE=figure;
hold on
for j=1:size(V2_p,2)
    if V2_p(end,j)>=0    
        plot3(X,time(j)*ones(size(X)),V2_p(:,j))
    else
        plot3(X,time(j)*ones(size(X)),-V2_p(:,j))
    end
end
set(gca,'fontname','times')
xlabel('x [m]','FontSize',16); ylabel('Training time [sec]','FontSize',16); zlabel('Amplitude [-]','FontSize',16)
box on
view(3)
saveas(DUE,fullfile([start,'/PLOT_conv'],'POM_2'),'fig')
saveas(DUE,fullfile([start,'/PLOT_conv'],'POM_2'),'pdf')