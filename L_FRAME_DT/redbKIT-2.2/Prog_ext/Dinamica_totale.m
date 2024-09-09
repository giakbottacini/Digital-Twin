%% Dynamical beam theory with moving_load

%% Plot deformata al variare della velocit? con carico ad x

clc
clear all
close all

n_modi=10;
L = 15;
x = L/2;
pi = 3.14159;
rho = 3286;
A = 5.9065;
xi=0:L/500:L;
F=185000;
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

E = 20.4*10^9;
I_x = 1.6636;
v_ = 160; %160-260
v = v_/3.6;
for i = 1 : n_modi
    phi(i,:) = sin(i*pi*xi/L);
    phi2(i,:)= -(i*pi/L)^2.*sin(i*pi*xi/L);
    w(i) = (i*pi/L)^2*sqrt((E*I_x)/(rho*A)); %omega
    T = L/v;
    t = 0 : T/500 : T;
    W(i) = i*pi/L*v; %Omega
    alpha(i) = W(i)/w(i);
    q(i,:) = (2*F/(rho*A*L)*(1/w(i)^2)*(1/(1-alpha(i)^2))).*(alpha(i).*sin(w(i).*t)-sin(W(i).*t));
end
for j = 1 : length(t)
spost(:,j) = phi'*q(:,j);
mom(:,j) = E*I_x*0.001.*(phi2'*q(:,j));
end

%soluz statica
vs = 0.001;
for i = 1 : n_modi
    Ts = L/vs;
    t = 0 : Ts/500 : Ts;
    Ws(i) = i*pi/L*vs; %Omega
    alphas(i) = Ws(i)/w(i);
    qs(i,:) = (2*F/(rho*A*L)*(1/w(i)^2)*(1/(1-alphas(i)^2))).*(alphas(i).*sin(w(i).*t)-sin(Ws(i).*t));
end
for j = 1 : length(t)
sposts(:,j) = phi'*qs(:,j);
end

%% displ
close(gcf)
figure
set(gca,'fontname','times');
xlabel('x [m]','FontSize',16); ylabel('w(x,t) [m]','FontSize',16);
xlim([0 15]); ylim([-5e-4 1e-4]);
box on
grid on
hold on
plot([0 L],[0 0],'k--','LineWidth',1,'HandleVisibility','off')
for j = 1 : length(t)
    if j==1
        a1(j)=plot(0,spost(1,j),'b','LineWidth',2);drawnow
        a2(j)=plot(0,sposts(1,j),'r','LineWidth',2);drawnow
        legend([num2str(v_),'Km/h'],'Static','Location','SouthEast')  
    end
    if j>1
    delete(a1(j-1));
    delete(a2(j-1));
    delete(a3(j-1));
    end
    a1(j)=plot(xi,spost(:,j),'b','LineWidth',2,'HandleVisibility','off');drawnow
    a2(j)=plot(xi,sposts(:,j),'r','LineWidth',2,'HandleVisibility','off');drawnow
    a3(j)=arrow([xi(j) 0.5e-4],[xi(j) 0]);drawnow
    M(j) = getframe(gcf);
end
all_valid = true;
flen = length(M);
for K = 1 : flen
  if isempty(M(K).cdata)
    all_valid = false;
    fprintf('Empty frame occurred at frame #%d of %d\n', K, flen);
  end
end
if ~all_valid
   error('Did not write movie because of empty frames')
end
  % create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
  writerObj.FrameRate = 100;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=2:length(M)
    % convert the image to a frame
    frame = M(i);    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

%% moment
close(gcf)
figure 
set(gca,'fontname','times','YDir','reverse');
xlabel('x [m]','FontSize',16); ylabel('M [kNm]','FontSize',16);
xlim([0 15]); ylim([-100 800]);
box on
grid on
hold on
plot([0 L],[0 0],'k--','LineWidth',1,'HandleVisibility','off')
for j = 1 : length(t)
    if j==1
        a1(j)=plot(0,mom(1,j),'b','LineWidth',2);drawnow
        a2(j)=plot(0,momstatico(j,1),'r','LineWidth',2);drawnow
        legend([num2str(v_),'Km/h'],'Static','Location','SouthEast')  
    end
    if j>1
    delete(a1(j-1))
    delete(a2(j-1))
    delete(a3(j-1))
    end
    a1(j)=plot(xi,mom(:,j),'b','LineWidth',2,'HandleVisibility','off');drawnow
    a2(j)=plot(xi,momstatico(j,:),'r','LineWidth',2,'HandleVisibility','off');drawnow
    a3(j)=arrow([xi(j) -80],[xi(j) 0]);drawnow
    M(j) = getframe(gcf);
end
all_valid = true;
flen = length(M);
for K = 1 : flen
  if isempty(M(K).cdata)
    all_valid = false;
    fprintf('Empty frame occurred at frame #%d of %d\n', K, flen);
  end
end
if ~all_valid
   error('Did not write movie because of empty frames')
end
  % create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
  writerObj.FrameRate = 100;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=2:length(M)
    % convert the image to a frame
    frame = M(i);    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);