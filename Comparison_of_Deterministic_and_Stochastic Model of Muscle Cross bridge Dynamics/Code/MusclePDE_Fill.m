% Deterministic PDE for muscle mechanics
VPVals = [];
VPTRUE = [];
for v = 5:5:150
% v = 50;
% N = 3200;
N = 1600;
x0 = 1;
MinusInf = -10;
% MinusInf = 0;
dx = (x0-MinusInf)/(N-1);
xx = (0:N-1)'*dx+MinusInf;
tf = 2;
dt = 0.0025*dx;
nSteps = tf/dt;
alpha0 = 14;
beta0 = 126;
IntegrationWts = [1/2 ones(1,N-2) 1/2]';
% Initial conditions
U = alpha0/(alpha0+beta0);
ubc = alpha0*(1-U)/v; 
u = ubc*exp(-beta0/v+beta0*xx/v);
% Solve advection ODE
for iT=1:nSteps
% Find U and evaluate the BC
U = sum(u.*IntegrationWts)*dx;
ubc = alpha0*(1-U)/v;
u(end) = ubc;
for j=1:N-1
% u(j)= u(j) + dt*(v*(u(j+1)-u(j))/dx-beta0*exp(xx(j))*u(j));
u(j)= u(j) + dt*(v*(u(j+1)-u(j))/dx-beta0*u(j));
end
% Advection for the first N-1 terms
end


hold on
% plot(xx,u)
% utrue = alpha0*(beta0/(alpha0+beta0))/v*exp(-beta0/v+beta0*xx/v);
% plot(xx,utrue)
% max(abs(u-utrue))
Ptrue = alpha0*p1/(alpha0+beta0)*...
    ((exp(mu*x0)-1)-mu*v/beta0)/(1+mu*v/beta0);
% Plot P and force velocity curve
p1 = 4;
mu = 0.322;
P = p1*sum(u.*(exp(mu*xx)-1).*IntegrationWts)*dx;
VPVals = [VPVals; v P];
% VPTRUE = [VPTRUE; v Ptrue];
end
% 
plot(VPVals(:,2),VPVals(:,1))
xlabel('$P$ per CB')
ylabel('$v$ per sarcomere')
a=xlim;
xlim([0 a(2)])
% hold on