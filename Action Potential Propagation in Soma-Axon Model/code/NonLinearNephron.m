% Solve the equation 
% a = e^alpha (1-alpha) for a given
avals = linspace(0.01,0.99,99);
alphavals = 0*avals;
for ia=1:length(avals)
    x = 0.5; % Initial guess
    a = avals(ia);
    v = f(x,a);
    vp = fprime(x);
    while (abs(v) > 1e-10)
        x = x - v/vp;
        v = f(x,a);
        vp = fprime(x);
    end
    % Compare to matlab 
    x2 = fzero(@(s) exp(s).*(1-s)-a,0.5);
    if (max(abs(x2-x)) > 1e-10)
        keyboard
    end
    alphavals(ia)=x;
end
plot(avals,alphavals)

function v = f(x,a)
    v = exp(x).*(1-x)-a;
end

function vp = fprime(x)
    vp = exp(x).*(1-x)-exp(x);
end

