clc; close all; clear all; beep off;

% parameters
p.rho_m = 0.1;
p.rho_p = 1; % bigger than p.rho_m
p.h = 1.5;
p.L = 1;

% don't change
p.g = 9.8;
p.a = log(p.rho_p/p.rho_m)/p.h;






left = 1.5;
right = 5;

fun = @(mu) (4-p.a^2-mu.^2).*tanh(mu*p.h/2)+2*p.a*mu;

f_left = fun(left);
f_right = fun(right);

if f_left*f_right >= 0
    error('root not straddled');
end




while right-left > 1e-8
    mid = 0.5*(left+right);
    f_mid = fun(mid);
    if sign(f_mid)== sign(f_left)
        left = mid;
    elseif sign(f_mid) == sign(f_right)
        right = mid;

    elseif sign(f_mid) == 0
        break
    else
        error('problem with bisection method');
    end
end


mu = mid


lam2 = -(2*p.g*tanh(mu*p.h/2))/(mu-p.a*tanh(mu*p.h/2))



