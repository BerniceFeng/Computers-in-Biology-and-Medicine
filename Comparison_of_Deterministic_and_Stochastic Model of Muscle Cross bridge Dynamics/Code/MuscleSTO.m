% clear all
% close all
global t_start v_zero;
t_start = 10;
v_zero = 50;
NumTrials = 5;
VPVals = [];
Nb = 1000;
alpha = 14;
beta = 126;
dt = 0.01 / (alpha+beta);
tmax = 30;
clockmax = ceil(tmax/dt);
x0 = 1;
x1 = 10;
p1 = 4;
mu = 0.322;
a = zeros(1,Nb);
% a(1)=1
x = zeros(1,Nb);
for iT=1:NumTrials
for clock = 1:clockmax
    x(find(a))=x(find(a))-v(clock*dt)*dt;
%     pc=(beta*exp(x)*dt).*a+(alpha*dt)*(1-a);
    pc=(beta*dt)*a+(alpha*dt)*(1-a);
    c=(rand(1,Nb)<pc)|(x>x1);
    a=xor(c,a);
    x(find(a&c))=x0;
    x(find(~a))=0;
    U=sum(a)/Nb;
    P=sum(p1*(exp(mu*x)-1))/Nb;
end
    binedges = 0:0.05:1;
    val(iT,:) = histcounts(x(find(a)),binedges)/Nb/0.05;
    Pvals4(iT) = P;
end

bincenters = (binedges(1:end-1)+binedges(2:end))/2;
errorbar(bincenters,mean(val),...
    2*std(val)/sqrt(NumTrials),'-o','LineWidth',2.0)
% plot(VPVals(:,2),VPVals(:,1))


% errorbar(mean(Pvals4),v_zero,...
% 2*std(Pvals4)/sqrt(NumTrials),'horizontal','-o','LineWidth',2.0)

% xlabel('$P$ per CB')
% ylabel('$v$ per sarcomere')
% a=xlim;
% xlim([0 a(2)])
% hold on