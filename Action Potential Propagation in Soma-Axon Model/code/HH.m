%filename HH.m
%numerical solution of the space-clamped Hodgkin-Huxley equations
clear all
clf
global check;
global t1p t2p ip; %parameters for the function izero(t)
in_HH
in_mhnv
for klok=1:klokmax
  t=klok*dt;                      %note time
  m=snew(m,alpham(vs),betam(vs),dt); %update m
  h=snew(h,alphah(vs),betah(vs),dt); %update h
  n=snew(n,alphan(vs),betan(vs),dt); %update n
  gNa=gNabar*(m^3)*h;    %sodium conductance
  gK =gKbar*(n^4);    %potassium conductance
  g=gNa+gK+gLbar;         %total conductance
  gE=gNa*ENa+gK*EK+gLbar*EL;         %gE=g*E
  %save old value of v for checking purposes:
  v_old=vs;
  %update v:
  vs=(vs+(dt/C)*(gE+izero(t)))/(1+(dt/C)*g);
  if(check)
    E=gE/g;
    chv=C*(vs-v_old)/dt+g*(vs-E)-izero(t)
  end
  %store results for future plotting:
  mhn_plot(:,klok)=[m h n]';
  v_plot(klok)=vs;
  t_plot(klok)=t;
end
subplot(2,1,1),plot(t_plot,v_plot)
subplot(2,1,2),plot(t_plot,mhn_plot)
