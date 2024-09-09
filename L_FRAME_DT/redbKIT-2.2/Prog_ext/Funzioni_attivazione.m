%% Funzioni di attivazione
clc
clear all
close all
syms z
logit = 1/(1+exp(-z));
tanh = 2/(1+exp(-2*z))-1;
dlogit = diff(logit);
dtanh = diff(tanh);
logit = matlabFunction(logit);
tanh = matlabFunction(tanh);
dlogit = matlabFunction(dlogit);
dtanh = matlabFunction(dtanh);
z=-5:0.01:5;
figure
hold on
plot(z,logit(z),'b','LineWidth',2)
plot(z,tanh(z),'k','LineWidth',2)
plot(z,max(0,z),'r','LineWidth',2)
set(gca,'fontname','times','FontSize',18);
xlim([-5 5]); ylim([-1.5 1.5]);
legend('Funzione logistica','Tangente iperbolica','Funzione rettificatore','Location','Southeast')
box on
grid on
figure
hold on
plot(z,dlogit(z),'b','LineWidth',2)
plot(z,dtanh(z),'k','LineWidth',2)
plot(z,max(sign(z),0),'r','LineWidth',2)
set(gca,'fontname','times','FontSize',18);
xlim([-5 5]); ylim([-0.5 1.5]);
legend('Funzione logistica','Tangente iperbolica','Funzione rettificatore','Location','Southeast')
box on
grid on
