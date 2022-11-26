function CV=CV_now(t,fmax)
%filename: CV_now.m
global T;
m1=1.32;
m2=27.4;
tao1=0.269*T;
tao2=0.452*T;
Emin=79.52;
Emax=5232;
tc=rem(t,T); %tc=time in the current cycle, 
             %measured from start of systole.
%instead of using the piecewise function, we use a continuous function
g1=(tc/tao1)^m1;
g2=(tc/tao2)^m2;
gT=(T/tao2)^m2;

if fmax ~=0
    k = (Emax-Emin)/fmax;
    Ev = k*(g1/(1+g1))*(1/(1+g2)-1/(1+gT))+Emin;
    CV = 1/Ev;
else 
    CV = 0.01;
end