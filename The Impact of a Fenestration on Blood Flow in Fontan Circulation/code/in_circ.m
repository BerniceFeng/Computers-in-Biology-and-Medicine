% %filename: in_circ.m  (initialization for circ)
% T =0.0125    %Duration of heartbeat (minutes)
% TS=0.0050    %Duration of systole   (minutes)
% tauS=0.0025  %CLV time constant during systole (minutes)
% tauD=0.0075  %CLV time constant during diastole (minutes)
% Rs=17.5      %Systemic resistance  (mmHg/(liter/minute))
% Rp= 1.79     %Pulmonary resistance (mmHg/(liter/minute))
% %Unrealistic valve resistances,
% %Chosen small enough to be negligible.
% RMi=0.01  %mitral valve resistance (mmHg/(liter/minute))
% RAo=0.01  %aortic valve resistance (mmHg/(liter/minute))
% RTr=0.01  %tricuspid valve resistance (mmHg/(liter/minute))
% RPu=0.01  %pulmonic valve  resistance (mmHg/(liter/minute))
% %The following values of Csa and Cpa are approximate.
% %They will need adjustment to make the systemic 
% %blood pressure be roughly 120/80 mmHg
% %and to make the pulmonary 
% %blood pressure be roughly 25/8 mmHg.
% Csa=0.00175  %Systemic  arterial compliance (liters/mmHg)
% Cpa=0.00412  %Pulmonary arterial compliance (liters/mmHg)
% Csv=1.75     %Systemic  venous compliance (liters/mmHg)
% Cpv=0.08     %Pulmonary venous compliance (liters/mmHg)
% CLVS=0.00003 %Min (systolic)  value of CLV (liters/mmHg)
% CLVD=0.0146  %Max (diastolic) value of CLV (liters/mmHg)
% CRVS=0.0002  %Min (systolic)  value of CRV (liters/mmHg)
% CRVD=0.0365  %Max (diastolic) value of CRV (liters/mmHg)
% Vsad=0.825   %Systemic arterial volume at P=0 (liters)
% Vpad=0.0382  %Pulmonary arterial volume at P=0 (liters)
% Vsvd=0       %Systemic venous volume at P=0 (liters)
% Vpvd=0       %Pulmonary venous volume at P=0 (liters)
% VLVd=0.027   %Left ventricular volume at P=0 (liters)
% %VRVd=0.027   %Right ventricular volume at P=0 (liters)
% dt=0.01*T    %Time step duration (minutes)
% %This choice implies 100 timesteps per cardiac cycle.
% klokmax=60*T/dt %Total number of timesteps 
% %This choice implies simulation of 15 cardiac cycles.
% %Assign an index to each compliance vessel 
% %of the model circulation:
%filename: in_circ.m  (initialization for circ)



T =0.016    %Duration of heartbeat (minutes)
m1=1.32
m2=27.4
tao1=0.269*T
tao2=0.452*T
Emin=79.52
Emax=5232
% Vd=0.028
TS=0.0050    %Duration of systole   (minutes)
tauS=0.0025  %CLV time constant during systole (minutes)
tauD=0.0075  %CLV time constant during diastole (minutes)
Rs=20.78      %Systemic resistance  (mmHg/(liter/minute))
Rp= 0.5517     %Pulmonary resistance (mmHg/(liter/minute))
%Unrealistic valve resistances,
%Chosen small enough to be negligible.
RMi=0.01  %mitral valve resistance (mmHg/(liter/minute))
RAo=0.01  %aortic valve resistance (mmHg/(liter/minute))
RTr=0.01  %tricuspid valve resistance (mmHg/(liter/minute))
RPu=0.01  %pulmonic valve  resistance (mmHg/(liter/minute))
%The following values of Csa and Cpa are approximate.
%They will need adjustment to make the systemic 
%blood pressure be roughly 120/80 mmHg
%and to make the pulmonary 
%blood pressure be roughly 25/8 mmHg.
Csa=0.0007333  %Systemic  arterial compliance (liters/mmHg)
Cpa=0.00412  %Pulmonary arterial compliance (liters/mmHg)
Csv=0.0990     %Systemic  venous compliance (liters/mmHg)
Cpv=0.01     %Pulmonary venous compliance (liters/mmHg)
% CLVS=0.001 %Min (systolic)  value of CLV (liters/mmHg)
% CLVD=0.0146  %Max (diastolic) value of CLV (liters/mmHg)
% CRVS=0.0002  %Min (systolic)  value of CRV (liters/mmHg)
% CRVD=0.0365  %Max (diastolic) value of CRV (liters/mmHg)
Vsad=0.7051   %Systemic arterial volume at P=0 (liters)
Vpad=0.0930  %Pulmonary arterial volume at P=0 (liters)
Vsvd=2.869       %Systemic venous volume at P=0 (liters)
Vpvd=0.1475       %Pulmonary venous volume at P=0 (liters)
VLVd=0.027   %Left ventricular volume at P=0 (liters)
%VRVd=0.027   %Right ventricular volume at P=0 (liters)
dt=0.01*T    %Time step duration (minutes)
%This choice implies 100 timesteps per cardiac cycle.
klokmax=ceil(60*T/dt) %Total number of timesteps 
%This choice implies simulation of 15 cardiac cycles.
%Assign an index to each compliance vessel 
%of the model circulation:
iLV=1
isa=2
isv=3
ipa=4
ipv=5
N=5
%Enter parameters and initial values 
%into correct slots in arrays.
%Note that the code that follows is independent 
%of the specific numbering scheme chosen above.
%Compliance vector:
C=zeros(N,1);  
%This makes C a column vector of length N.
fmax=0
C(iLV)=CV_now(0,fmax);  %initial value
C(isa)=Csa;
C(isv)=Csv;
% C(iRV)=CV_now(0,CRVS,CRVD);  %initial value
C(ipa)=Cpa;
C(ipv)=Cpv;
C  %This writes the result on the screen.
%Pressure vector (initial values) at end of diastole:
P=zeros(N,1);  
%This makes P a column vector of length N.
P(iLV)= 9;
P(isa)=80;
P(isv)= 2;
% P(iRV)= 2;
P(ipa)= 8;
P(ipv)= 5;
P  %This writes the result on the screen.
%Vector of dead volumes (volume at zero pressure);
%Note: Vd is only needed for output purposes.  
%It drops out of the equations we solve for P, 
%but we need it if we want to output (e.g., plot)  
%the volume of any compliance vessel.
Vd=zeros(N,1);  
%This makes Vd a column vector of length N.
Vd(iLV)=VLVd;
Vd(isa)=Vsad;
Vd(isv)=Vsvd;
% Vd(iRV)=VRVd;
Vd(ipa)=Vpad;
Vd(ipv)=Vpvd;
Vd  
% Rig the volume
Vset = 5.0;
Vtotal = sum(Vd+C.*P);
P(isa) = P(isa) + (Vset-Vtotal)/C(isa);
newV = sum(Vd+C.*P)
%This writes the results on the screen.
%Conductance matrix:
G=zeros(N,N);  
%This makes G an NxN matrix filled with zeros.
%Any element of G that is not explicitly 
%made nonzero remains zero, 
%thus modeling an infinite resistance connection, 
%that is, no connection at all.
G(iLV,isa)=1/RAo;  %But G(isa,iLV)=0 (no leak)
G(isa,isv)=1/Rs;   %no valve
G(isv,isa)=1/Rs;   %no valve
G(isv,ipa)=1/RTr;  %But G(iRV,isv)=0; (no leak)
G(ipa,isv)=1/RTr;
% G(iRV,ipa)=1/RPu;  %But G(ipa,iRV)=0; (no leak)
G(ipa,ipv)=1/Rp;   %no valve
G(ipv,ipa)=1/Rp;   %no valve
G(ipv,iLV)=1/RMi;  %But G(iLV,ipv)=0; (no leak)


% RVSD = 0.01;
% G(iLV,iRV)=1/RVSD; % VSD
% G(iRV,iLV)=1/RVSD;


G  %This writes the result on the screen.
%Matrix of initial valve states:
S=zeros(N,N)   
%This makes S an NxN matrix filled with zeros
%(and writes it on the screen).
%Start with all valves closed.
%Valves will adjust to pressures 
%during first time step.
%Initialize arrays to store data for plotting:
t_plot=zeros(1,klokmax);
C_plot=zeros(N,klokmax);
P_plot=zeros(N,klokmax);
%Other variables that we might want to plot 
%can be found from these.
%For self-checking in P_new, set CHECK=1.
%To skip self-checking set CHECK=0.
%(should be much faster with CHECK=0)
CHECK=0;
%Initialize flow computation (for output purposes only)
%assign an index to each flow of interest:
jAo=1;
js =2;
jTr=3;
%jPu=4;
jp =4;
jMi=5;
% jVSD=6;
Nflows=5;
%note index of upstream and downstream chamber 
%for each flow:
iU=zeros(Nflows,1);
iD=zeros(Nflows,1);
iU(jAo)=iLV;
iD(jAo)=isa;
iU(js )=isa;
iD(js )=isv;
iU(jTr)=isv;
% iD(jTr)=iRV;
% iU(jPu)=iRV;
% iD(jPu)=ipa;
iD(jTr)=ipa;
iU(jp )=ipa;
iD(jp )=ipv;
iU(jMi)=ipv;
iD(jMi)=iLV;

% iU(jVSD)=iLV;
% iD(jVSD)=iRV;


%extract the conductances from the matrix G:
Gf=zeros(Nflows,1);
Gr=zeros(Nflows,1);
for j=1:Nflows
  Gf(j)=G(iU(j),iD(j)); %forward conductance
  Gr(j)=G(iD(j),iU(j)); %reverse conductance
end
%create arrays to store current pressure differences 
%and history over time of the net flows:
Pdiff=zeros(Nflows,1);
Q_plot=zeros(Nflows,klokmax);
