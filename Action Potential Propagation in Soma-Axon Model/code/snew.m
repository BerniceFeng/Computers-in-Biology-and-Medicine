function s=snew(s_old,alpha,beta,dt)
%filename: snew.m
s=(s_old+dt*alpha)./(1+dt*(alpha+beta));
end