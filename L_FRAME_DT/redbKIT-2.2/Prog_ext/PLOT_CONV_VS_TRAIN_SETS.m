%% PLOT_CONV_VS_TRAIN_SETS
clc
clear all
close all

% N_sets = [0, 5, 20, 40, 80, 240];
% N_modi = [0, 13, 21, 20, 21, 22];
% Accuracy=[0, 45, 65, 70, 76];
% 
% 
% figure
% hold on
% xp=1:1:300;
% vq2 = interp1(N_sets,N_modi,xp,'linear','extrap');
% plot(xp,vq2,':.');
% plot(N_sets(2:end),N_modi(2:end),'rs','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
% set(gca,'fontname','times');
% xlabel('Trainig sets [-]','FontSize',16); ylabel('POMs [-]','FontSize',16);
% xlim([0 300]); ylim([0 40]);
% box on
% grid on


N_sets = [0.1, 5, 20, 40, 80, 240];
N_modi = [0, 13, 21, 20, 21, 22];

fun = @(x)sseval(x,N_sets,N_modi);
N_fit=0.01:2:300;
x0 = rand(2,1);
bestx = fminsearch(fun,x0);
A = bestx(1);
lambda = bestx(2);
yfit = A*log10(lambda*N_fit);
figure
plot(N_fit,yfit,':.');
hold on
plot(N_sets(2:end),N_modi(2:end),'rs','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
set(gca,'fontname','times');
xlabel('Trainig sets [-]','FontSize',16); ylabel('POMs [-]','FontSize',16);
xlim([0 300]); ylim([0 30]);
legend('Fitted Curve','Tests')
box on
grid on

function sse=sseval(x,ascissa,ordinata)
    A=x(1);
    lambda=x(2);
    sse=sum((ordinata-A*log10(lambda*ascissa)).^2);
end

% 
% figure
% hold on
% xp=1:1:300
% vq2 = interp1(N_sets,Accuracy,xp,'poly3','extrap');
% plot(xp,vq2,':.');
% plot(N_sets(2:end),Accuracy(2:end),'rs','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);
% set(gca,'fontname','times');
% xlabel('Trainig sets [-]','FontSize',16); ylabel('Test Accuracy [%]','FontSize',16);
% xlim([0 300]); ylim([0 100]);
% box on
% grid on