%filename: in_mhnv.m
%this script simulates holding the
%voltage constant for a long time
%prior to t=0 at the value v=vhold,
%and then applying a current shock
%at t=0 to step the voltage
%to a new level, v=vstart.
%The gating variables m,h,n
%cannot change suddenly and therefore
%remain constant during the current shock.
%
%simulate time prior to t=0:
vs=vhold
%set m,h,n equal to their steady values
%under constant-v conditions:
ms=alpham(vs)/(alpham(vs)+betam(vs))
hs=alphah(vs)/(alphah(vs)+betah(vs))
ns=alphan(vs)/(alphan(vs)+betan(vs))
%
%now let voltage jump to its value
%just after t=0, without making
%any further change in m,h,n:
vs=vstart
