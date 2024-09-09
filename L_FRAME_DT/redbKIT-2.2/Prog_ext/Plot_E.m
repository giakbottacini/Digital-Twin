%% Plot_E
clc
clear all
close all

E_15=210e+09
alpha=2.4e-4

T=0:0.1:40;
T_=T-15;
X_=0:1:15;
X_P=E_15.*ones(length(X_),1);
Y_=0:10000000:E_15;
Y_P=15.*ones(1,length(Y_));
E_correct=E_15*(1-alpha*T_);
figure
ylim([2.08e+11 2.11e+11])
xlim([0 40])
hold on
plot(T,E_correct)
plot(X_,X_P,'k--')
plot(Y_P,Y_,'k--')
xlabel('T [°C]','FontSize',16)
ylabel('E(T) [Pa]','FontSize',16)
set(gca,'fontname','times')