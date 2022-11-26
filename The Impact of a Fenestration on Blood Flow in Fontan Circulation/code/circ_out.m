%filename: circ_out.m
%script to plot results of computer simulation 
%of the entire circulation.

figure(1)
subplot(2,1,1),plot(t_plot,P_plot([iLV,isa],:))
xlim([0.5 0.8])
title('Single ventricle and sysmetic arteries pressure','FontSize',16)
xlabel('Time ($minutes$)','interpreter','latex','fontweight','bold','fontsize',14)
ylabel('Pressure ($mmHg$)','interpreter','latex')
legend('Single ventricle','Sysmetic arteries')

subplot(2,1,2),plot(t_plot,P_plot(isv,:))
title('Systemic veins pressure','FontSize',16)
xlim([0.5 0.8])
xlabel('Time ($minutes$)','interpreter','latex','fontweight','bold','fontsize',14)
ylabel('Pressure ($mmHg$)','interpreter','latex')

%systemic and pulmonary flows:
figure(2)
subplot(3,1,1),plot(t_plot,Q_plot([jAo,js],:))
xlim([0.2 1])
subplot(3,1,2),plot(t_plot,Q_plot([jTr,jp],:))
xlim([0.2 1])
subplot(3,1,3),plot(t_plot,Q_plot([jMi,jpsv],:))
xlim([0.2 1])

figure(3)
plot(V_plot(iLV,:),P_plot(iLV,:))
title('Single ventricle PV loop','FontSize',16)
xlabel('Volume ($liters$)','interpreter','latex','fontweight','bold','fontsize',14)
ylabel('Pressure ($mmHg$)','interpreter','latex')

% fenestration flow
figure(4)
plot(t_plot,Q_plot(jpsv,:))
xlim([0.2 1])
