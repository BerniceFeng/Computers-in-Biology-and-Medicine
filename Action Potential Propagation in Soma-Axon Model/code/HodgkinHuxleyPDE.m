% Solve the PDE
% C v_t + g(v-E) = r/(2*rho) v_xx
% Define parameters
clear
clc
close all

global check;
global t1p t2p ip; %parameters for the function izero(t)
in_HH
in_mhnv
r = 0.0238; % cm
rho = 0.0354; % mV/uA
Cs=1.0;     %(muF/cm^2)
Cd=1.0;     %(muF/cm^2)
gNabar=120; %((muA/mV)/cm^2) max possible Na+ conductance per unit area
gKbar=36;  %((muA/mV)/cm^2) max possible K+ conductance per unit area
gLbar=0.3;  %((muA/mV)/cm^2) leakage conductance per unit area
ENa = 45;   %(mV)
EK = -82;   %(mV)
EL = -59;   %(mV)
% Define the x grid and initial conditions
L = 4; % cm
A = 100; % cm^2
N = 101;
h = L/(N-1);
x = (0:N-1)'*h;
% Matrix for solving the linear system
M = sparse(N+1,N+1);
% Second derivative matrix (from last time)
L = zeros(N,N);
vd = -70*ones(N,1); 
m = alpham(vd)./(alpham(vd)+betam(vd));
n = alphan(vd)./(alphan(vd)+betan(vd));
jj = alphah(vd)./(alphah(vd)+betah(vd)); % CALL H GATES J GATES - avoids conflict with dx
dt = 1e-2;
tf = 8;
clockmax=tf/dt;
allt = zeros(clockmax);
allv = zeros(clockmax,N);
allm = allv;
alln = allv;
allh = allv;
Z = zeros(N+1,1);

for clock=1:clockmax
    t=clock*dt;                      %note time
    ms=snew(ms,alpham(vs),betam(vs),dt); %update m
    hs=snew(hs,alphah(vs),betah(vs),dt); %update h
    ns=snew(ns,alphan(vs),betan(vs),dt); %update n
    m=snew(m,alpham(vd),betam(vd),dt); %update m
    jj=snew(jj,alphah(vd),betah(vd),dt); %update h
    n=snew(n,alphan(vd),betan(vd),dt); %update n
    gNas=gNabar*(ms^3)*hs;    %sodium conductance
    gKs =gKbar*(ns^4);    %potassium conductance
    gs=gNas+gKs+gLbar;         %total conductance
    gEs=gNas*ENa+gKs*EK+gLbar*EL;         %gE=g*E
    gNa=gNabar.*(m.^3).*jj;%gNabar;%.*(m.^3).*j;    %sodium conductance
    gK =gKbar.*(n.^4);%gKbar;%.*(n.^4);    %potassium conductance
    g=gNa+gK+gLbar;         %total conductance
    gE=gNa*ENa+gK*EK+gLbar*EL;         %gE=g*E

    vvector = [vs; vd];
    vold = vvector;
    M(1,1) = Cs/dt + gs + pi*r^2/(rho*h*A);
    M(1,2) = -pi*r^2/(rho*h*A);
    Z(1) = Cs/dt*vvector(1) + gEs + izero(t);
    for j=2:N
        M(j,j) = Cd/dt + g(j-1) + r/(rho*h^2);
        M(j,j+1) = -r/(2*rho*h^2);
        M(j,j-1) = -r/(2*rho*h^2);
    end
    M(N+1,N) = -r/(rho*h^2);
    M(N+1,N+1) = Cd/dt + g(N) + r/(rho*h^2);
    for i=2:N+1
        Z(i) = Cd/dt*vvector(i) + gE(i-1);
    end
    vvector = M \ Z;
    vs = vvector(1);
    vd = vvector(2:N+1);
    
%     Check
      if(check)
        I0=-pi*r^2/rho*(vvector(2)-vvector(1))/h;
        RHS1=izero(t)-I0;
        LHS1=Cs*(vvector(1)-vold(1))/dt+gs*vvector(1)-gEs;
        LHS1-RHS1;
        L=zeros(N,N+1);
        for iRow=1:N-1
        L(iRow,iRow+1)=-2/h^2;
        L(iRow,iRow)=1/h^2;
        L(iRow,iRow+2)=1/h^2;
        end
        L(N,N+1)=-2/h^2;
        L(N,N)=2/h^2;
        LHS2=Cd*(vvector(2:end)-vold(2:end))/dt+g.*(vvector(2:end))-gE;
        RHS2=r/(2*rho)*(L*vvector);
        max(abs(RHS2-LHS2));     
      end

    % Plot results
    allt(clock)=clock;
    allv(clock,:)=vd;
    allh(clock,:)=j;
    allm(clock,:)=m;
    alln(clock,:)=n;
    figure(1)
    subplot(3,1,1),plot(x,vd)   
    ylim([-85 65])
    xlim([-0.5 4.3])
    title(strcat('t=',num2str(t),', Vs=',num2str(vs)))
    subplot(3,1,2)
    rsoma = 0.1;
    xx=[-rsoma 0:h:4]';
    scatter(xx,0*xx,[1000; 36*ones(N,1)],vvector,'filled');
    caxis([-80, 40])
    xlim([-0.5 4.3])
    subplot(3,1,3),plot(x,m,x,jj,x,n)
    ylim([0 1])
    legend('m','h','n')
    xlabel('t')
    xlim([-0.5 4.3])
    filename = 'action.gif';
    del = 0.1; % time between animation frames
    drawnow 
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if clock == 1
    imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
    else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
    end
end

% figure(2)
% plot(allv(:,1));
% hold on
% plot(allv(:,(N+1)/2))
% plot(allv(:,N))
% hold off
% legend('start','middle','end')
% title("Voltage of Axon with Area 10")
% % title("Voltage of Axon with Current 80 Applied to Soma")
% drawnow