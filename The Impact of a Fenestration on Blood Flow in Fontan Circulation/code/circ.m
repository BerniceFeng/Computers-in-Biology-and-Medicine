%filename:  circ.m
clear all % clear all variables
clf       % and figures
global T TS tauS tauD;
global G dt CHECK N;
% in_circ  %initialize
pmax = 5;
Nflows=6;
R_plot=zeros(1,pmax);
Qm_plot=zeros(Nflows,pmax);
for loop=1:pmax
    loop
    revised
    R_plot(loop)=Rpsv;
    t_all = [];
    for klok=1:floor(T/dt)
        t=klok*dt;
        t_all = [t_all; t];
    end
    t_all;
    fmax = 0;
    gT=(T/tao2)^m2;
    for i = 1:length(t_all)
        g1t = (t_all(i)/tao1)^m1;
        g2t = (t_all(i)/tao2)^m2;
        f = (g1t/(1+g1t))*(1/(1+g2t)-1/(1+gT));
        if f>fmax
            fmax=f;
        end
    end
    
    for klok=1:klokmax
      t=klok*dt;
    
      P_old=P;
      C_old=C;
      %find current values of left and right 
      %ventricular compliance and store each 
      %of them in the appropriate slot in the array C:
      C(iLV)=CV_now(t,fmax);
      %find self-consistent valve states and pressures:
      set_valves
      %store variables in arrays for future plotting:
      t_plot(klok)=t;
      C_plot(:,klok)=C;
      P_plot(:,klok)=P;
      V_plot(:,klok)=Vd+C.*P;
      Pdiff=P(iU)-P(iD); %pressure differences 
                         %for flows of interest:
      Q_plot(:,klok)=(Gf.*(Pdiff>0)+Gr.*(Pdiff<0)).*Pdiff;
     
      %(the net flow is computed in each case)
    end
    %plot results:
%     circ_out
    Qm_plot(:,loop)=mean(Q_plot(:,end-2000:end)');
    Pm_plot(:,loop)=mean(P_plot(:,end-2000:end)');
end
figure(1)
subplot(3,1,1),plot(R_plot,Qm_plot(2,:),'LineWidth',1)
title('Systemic Flow','FontSize',16,'FontName','Arial')
xlabel('Resistance of fenestration ($mmHg$ $min$ $L^{-1}$)','interpreter','latex','fontweight','bold','fontsize',14)
ylabel('($L$ $min^{-1}$)','interpreter','latex')
subplot(3,1,2),plot(R_plot,Qm_plot(3,:),'LineWidth',1)
title('Pulmonary Flow','FontSize',16)
xlabel('Resistance of fenestration ($mmHg$ $min$ $L^{-1}$)','interpreter','latex','fontweight','bold','fontsize',14)
ylabel('($L$ $min^{-1}$)','interpreter','latex')
subplot(3,1,3),plot(R_plot,Qm_plot(6,:),'LineWidth',1)
title('Fenestration Flow','FontSize',16)
xlabel('Resistance of fenestration ($mmHg$ $min$ $L^{-1}$)','interpreter','latex','fontweight','bold','fontsize',14)
ylabel('($L$ $min^{-1}$)','interpreter','latex')
figure(2)
subplot(2,1,1),plot(R_plot,Pm_plot(3,:),'LineWidth',1)
title('Systemic Veins pressure','FontSize',16,'FontName','Arial')
xlabel('Resistance of fenestration ($mmHg$ $min$ $L^{-1}$)','interpreter','latex','fontweight','bold','fontsize',14)
ylabel('($mmHg$)','interpreter','latex')
subplot(2,1,2),plot(R_plot,Pm_plot(5,:),'LineWidth',1)
title('pulmonary Veins pressure','FontSize',16,'FontName','Arial')
xlabel('Resistance of fenestration ($mmHg$ $min$ $L^{-1}$)','interpreter','latex','fontweight','bold','fontsize',14)
ylabel('($mmHg$)','interpreter','latex')
figure(3)
plot(R_plot,Pm_plot(1,:),'LineWidth',1)
figure(4)
plot(R_plot,Pm_plot(2,:),'LineWidth',1)
figure(5)
plot(R_plot,Pm_plot(4,:),'LineWidth',1)