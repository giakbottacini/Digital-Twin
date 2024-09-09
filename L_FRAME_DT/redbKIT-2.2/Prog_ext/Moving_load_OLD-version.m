%% Dynamical beam theory with moving_load

%% Plot deformata al variare della velocità con carico ad x

clc
clear all
close all

n_modi=3;
L = 15;
x = L/2;
rho = 3286;
A = 5.9065;
xi=0:L/100:L;
index=find(xi==x);%Poisizione carico
F=18500;
E = 20.4*10^9;
I_x = 1.6636;
v_ = 160; %160-260
v = v_/3.6;
for i = 1 : n_modi
    phi(i,:) = sin(i*pi*xi/L);
    w(i) = (i*pi/L)^2*sqrt((E*I_x)/(rho*A)); %omega
    T = L/v;
    t = 0 : T/100 : T;
    W(i) = i*pi/L*v; %Omega
    alpha(i) = W(i)/w(i);
    q(i,:) = (2*F/(rho*A*L)*(1/w(i)^2)*(1/(1-alpha(i)^2))).*(alpha(i).*sin(w(i).*t)-sin(W(i).*t));
end
for j = 1 : length(t)
spost(:,j) = phi'*q(:,j);
end

%soluz statica
vs = 0.01;
for i = 1 : n_modi
    Ts = L/vs;
    t = 0 : Ts/100 : Ts;
    Ws(i) = i*pi/L*vs; %Omega
    alphas(i) = Ws(i)/w(i);
    qs(i,:) = (2*F/(rho*A*L)*(1/w(i)^2)*(1/(1-alphas(i)^2))).*(alphas(i).*sin(w(i).*t)-sin(Ws(i).*t));
end
for j = 1 : length(t)
sposts(:,j) = phi'*qs(:,j);
end;

figure
set(gca,'fontname','times');
xlabel('Trainig sets [-]','FontSize',16); ylabel('POMs [-]','FontSize',16);
xlim([0 15]); ylim([-5e-5 0]);
box on
grid on
hold on
for j = 1 : length(t)
    a1(j)=plot(xi,spost(:,j),'bs-','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4,'HandleVisibility','off');drawnow
    a2(j)=plot(xi,sposts(:,j),'r','LineWidth',2,'HandleVisibility','off');drawnow
    if j>1
    delete(a1(j-1))
    delete(a2(j-1))
    end
    if j == length(t)
        plot(xi,spost(:,j),'bs-','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);drawnow
        plot(xi,sposts(:,j),'r','LineWidth',2);drawnow
%     a(i)=plot(data_lineB.X(i).*0.001,data_lineB.Y(i),'ko', 'MarkerFaceColor','k');drawnow
%     F(i) = getframe(gcf)
    end
end


% figure
% hold on
% plot(xi,spost,'bs-','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)
% plot(xi,sposts,'r','LineWidth',2)
% set(gca,'fontname','times');
% xlabel('Trainig sets [-]','FontSize',16); ylabel('POMs [-]','FontSize',16);
% xlim([0 15]); %ylim([0 30]);
% legend([num2str(v),'Km/h'],'Static')
% box on
% grid on

%% Plot spostamento in mezzeria al variare della velocità

clc
clear all
close all

n_modi=5;
L = 15;
x = L/2; %Posizione osservata
rho = 3286;
A = 5.9065;
xi=0:L/100:L;
index=find(xi==x);
F=18500;
E = 20.4*10^9;
I_x = 1.6636;
v_ = [160 180 200 220 240 610 0.1];
v = v_./3.6;
for i = 1 : n_modi
    phi(i,:) = sin(i*pi*xi/L);
    w(i) = (i*pi/L)^2*sqrt((E*I_x)/(rho*A)); %omega
    for j = 1 : length(v)
        T(j) = L/v(j);
        t(j,:) = 0 : T(j)/100 : T(j);
        W(i,j) = i*pi/L*v(j); %Omega
    end
end

for j = 1 : length(v)
    for i = 1 : n_modi
        alpha(i,j) = W(i,j)/w(i);
    end
end

for j = 1 : length(v)
    for i = 1 : n_modi
        q(i,j,:) = (2*F/(rho*A*L)*(1/w(i)^2)*(1/(1-alpha(i,j)^2))).*(alpha(i,j).*sin(w(i).*t(j,:))-sin(W(i,j).*t(j,:)));
    end
end
    
for j = 1 : length(v)
    for k = 1 :size(t,2)
        spost(j,k) = phi(:,index)'*q(:,j,k);
    end
end  

figure
hold on
plot(t(1,:)./T(1),spost(1,:),'bs-','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)
plot(t(2,:)./T(2),spost(2,:),'bs-','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4)
plot(t(3,:)./T(3),spost(3,:),'bs-','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',4)
plot(t(4,:)./T(4),spost(4,:),'bs-','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',4)
plot(t(5,:)./T(5),spost(5,:),'bs-','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',4)
plot(t(6,:)./T(6),spost(6,:),'bs-','MarkerEdgeColor','k','MarkerFaceColor','c','MarkerSize',4)
plot(t(7,:)./T(7),spost(7,:),'r','LineWidth',2)
set(gca,'fontname','times');
xlabel('Trainig sets [-]','FontSize',16); ylabel('POMs [-]','FontSize',16);
xlim([0 1]); %ylim([0 30]);
legend([num2str(v_(1)),'Km/h'],[num2str(v_(2)),'Km/h'],[num2str(v_(3)),'Km/h'],[num2str(v_(4)),'Km/h'],[num2str(v_(5)),'Km/h'],[num2str(v_(6)),'Km/h'],'Static')
box on
grid on