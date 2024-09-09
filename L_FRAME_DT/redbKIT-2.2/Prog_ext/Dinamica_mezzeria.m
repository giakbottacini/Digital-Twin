%% Plot spostamento in mezzeria

clc
clear all
close all

n_modi=10;
L = 15;
x = L/2; %Posizione osservata
rho = 3286;
A = 5.9065;
pi = 3.14159;
xi=0:L/100:L;
index=find(xi==x);
F=18500;
E = 20.4*10^9;
I_x = 1.6636;
v_ = [160 600 0.1];
v = v_./3.6;

for i = 1 : n_modi
    phi(i,:) = sin(i*pi*xi/L);
    phi2(i,:)= -(i*pi/L)^2.*sin(i*pi*xi/L);
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
plot([0 1],[0 0],'k--','LineWidth',1,'HandleVisibility','off')
hold on
plot(t(1,:)./T(1),spost(1,:),'b','LineWidth',2)
plot(t(2,:)./T(2),spost(2,:),'g','LineWidth',2)
plot(t(3,:)./T(3),spost(3,:),'r','LineWidth',2)
set(gca,'fontname','times');
xlabel('$\frac{t_{i}}{T_{i}}$ [-]','FontSize',16,'Interpreter','latex'); ylabel('w($\frac{L}{2}$,t) [m]','FontSize',16,'Interpreter','latex');
xlim([0 1]); %ylim([0 30]);
legend([num2str(v_(1)),'Km/h'],[num2str(v_(2)),'Km/h'],'Static','Location','SouthWest')
box on
grid on

%% Plot momento in mezzeria

xF = L.*ones(1,length(xi))-xi;
for i=1:length(xF)
    for j=1:length(xi)
momsA(i,j) = 0.001.*F.*(xF(i)./L).*xi(j);
momsB(i,j) = 0.001.*F.*((xF(i)./L).*xi(j)-(xi(j)-xi(i)));
    end
end
for i=1:length(xF)
    for j=1:length(xi)
        if j<=i
            momstatico(i,j)= momsA(i,j);
        else
            momstatico(i,j)= momsB(i,j);
        end
    end
end
asseF=xi./L;

for j = 1 : length(v)
    for k = 1 :size(t,2)
        mom(j,k) = E*I_x.*0.001.*(phi2(:,index)'*q(:,j,k));
    end
end  

figure
plot([0 1],[0 0],'k--','LineWidth',1,'HandleVisibility','off')
hold on
plot(t(1,:)./T(1),mom(1,:),'b','LineWidth',2)
plot(t(2,:)./T(2),mom(2,:),'g','LineWidth',2)
plot(asseF,momstatico(round(length(momstatico)/2),:),'r','LineWidth',2)
set(gca,'fontname','times','YDir','reverse');
xlabel('$\frac{t_{i}}{T_{i}}$ [-]','FontSize',16,'Interpreter','latex'); ylabel('M($\frac{L}{2}$,t) [kNm]','FontSize',16,'Interpreter','latex');
xlim([0 1]); ylim([-20 100]);
legend([num2str(v_(1)),'Km/h'],[num2str(v_(2)),'Km/h'],'Static','Location','SouthWest')
box on
grid on
