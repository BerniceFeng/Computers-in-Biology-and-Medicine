function a=alphan(v)
%filename: alphan.m
theta=(v+60)/10;
a = 0.1*(theta==0)+0.1*theta./(1-exp(-theta)).*(theta~=0);
end
