clc; close all; clear all;

%
% system parameters
%

p.a = 0.5;

%
% directory management
%

curr_dir = cd; cd('../code');

% -------------------------------------------------------------------------
% guess for the periodic solution - shooting to find approximate periodic
% -------------------------------------------------------------------------

ode_options = odeset('RelTol',1e-10,'AbsTol',1e-10);
y0 = [1;0];
ode_fun = @(t,z)[z(2);2*p.a*z(2)*(z(1)^2 - 1) - z(1)];
sol = ode15s(ode_fun,[0,-600],y0,ode_options);

y0 = deval(sol,sol.x(end));
sol = ode15s(ode_fun,[0,8],y0,ode_options);
left_X = 4;
right_X = 8;
while right_X-left_X > 1e-12
   mid_X = (left_X+right_X)/2;
   temp = deval(sol,mid_X);
   if temp(1) < y0(1)
       left_X = mid_X;
   else
       right_X = mid_X;
   end
end

sol = ode15s(ode_fun,[0,mid_X],y0,ode_options);

% -------------------------------------------------------------------------
% Solve for the periodic solution
% -------------------------------------------------------------------------

bvp_options = bvpset('RelTol',1e-12,'AbsTol',1e-12,'Nmax', 20000);

ode_fun = @(t,y,X)ode_periodic(t,y,X,p);
ode_fun_per = ode_fun;
perp = [0.4,0.1];
bc_fun = @(ya,yb,X)[ya-yb;perp*[ya(2);2*p.a*ya(2)*(ya(1)^2 - 1) - ya(1)]];

% solve the periodic with the initial guess coming from the shooting method
x_dom = linspace(0,mid_X,31);
pre_guess = @(t)deval(sol,t);
solinit = bvpinit(x_dom,pre_guess,1);
sol_per = bvp5c(ode_fun,bc_fun,solinit,bvp_options);

% use the bvp solution as an initial guess with domain [0,1] to get period
x_dom = linspace(0,1,31);
pre_guess = @(t)deval(sol_per,t*sol_per.x(end));
solinit = bvpinit(x_dom,pre_guess,sol_per.x(end));
sol_per = bvp5c(ode_fun,bc_fun,solinit,bvp_options);

% use the bvp solution as initial guess with the actual domain [0,period]
x_dom = linspace(0,sol_per.parameters,31);
pre_guess = @(t)deval(sol_per,t/sol_per.parameters);
solinit = bvpinit(x_dom,pre_guess,1);
sol_per = bvp5c(ode_fun,bc_fun,solinit,bvp_options);

hold on;
plot(sol_per.y(1,:),sol_per.y(2,:),'-b','LineWidth',2);
h = xlabel('x');
set(h,'FontSize',18);
h = ylabel('y');
set(h,'FontSize',18);
h = zlabel('Action');
set(h,'FontSize',18);
drawnow;

cd('../data');
save(replace(['profile_a_',num2str(p.a)],'.','P'),'sol_per');
cd(curr_dir);