function a=alpham(v)
%filename: alpham.m
theta=(v+45)/10;
a = (theta==0)+1.0*theta./(1-exp(-theta)).*(theta~=0);