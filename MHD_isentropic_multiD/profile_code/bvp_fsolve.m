function sol = bvp_fsolve(ode,bc,solinit,options)

% left and right end points of the interval [a,b]
sol.a = solinit.boundaries(1);
sol.b = solinit.boundaries(end);

% ODE
sol.ode = ode;
% BCs

sol.bc = bc;
% number of interpolation nodes for Chebyshev polynomial
sol.N = options.N;
[x0, DM] = chebdif(sol.N, 1); 
sol.xtilde = x0; % xtilde \in [-1,1]
sol.x = 0.5*(sol.a+sol.b)+0.5*(sol.a-sol.b)*sol.xtilde; % x \in [a,b]
sol.Tx = (2/(sol.a-sol.b))*DM(:,:,1); % derivative of interpolated function

% dimension of the first order ODE system
sol.dim = length(solinit.guess(sol.a));

% initial guess of solution
y0 = zeros(sol.dim,sol.N);
for j = 1:sol.N
   y0(:,j) = solinit.guess(sol.x(j)); 
end

if isfield(solinit,'params')
    params = solinit.params;
    if size(params,2) > size(params,1)
       params = params.'; 
    end
else
   params = []; 
end

ind = 0:1:sol.N-1;
Tcf = (2/(sol.N-1))*cos((pi/(sol.N-1))*(ind.'*ind));
Tcf(1,:) = 0.5*Tcf(1,:);
Tcf(end,:) = 0.5*Tcf(end,:);
Tcf(:,1)  = 0.5*Tcf(:,1);
Tcf(:,end) = 0.5*Tcf(:,end);

sol.Tcf = Tcf;



u0 = [reshape(y0.',sol.N*sol.dim,1); params];

[u,f,exitflag] = fsolve(@(u) ode_fun(u,ode,bc,sol),u0,options.algorithm);

sol.exitflag = exitflag;
sol.y = real(reshape(u(1:sol.dim*sol.N),sol.N,sol.dim).');
sol.params = real(u(sol.dim*sol.N+1:end));

sol.cf = Tcf*sol.y.';
sol.der_cf = Tcf*(sol.Tx*sol.y.');


% functions a la semi-object oriented fashion
sol.deval = @(x)local_deval(x,sol);

sol.check_err = @(num_pts)check_err(num_pts,sol);

sol.plot = @(num_pts)plot_profile(num_pts,sol);

sol.solver = 'bvp_fsolve';

% -------------------------------------------------------------------------
% plot
% -------------------------------------------------------------------------
function plot_profile(num_pts,sol)

dom = linspace(sol.a,sol.b,num_pts);

figure;
hold on;
y = sol.deval(dom);
plot(dom,y,'LineWidth',2);
% h = legend('R','S','Nr','Ni','Location','Best');
% set(h,'FontSize',18);
% h = xlabel('x');
% set(h,'FontSize',18);
% h = gca;
% set(h,'FontSize',18);
drawnow;


% -------------------------------------------------------------------------
% check_err 
% -------------------------------------------------------------------------

function [ode_err,bc_err] = check_err(num_pts,sol)

% check residual error of the profile ODE and BCs
ode_err = 0;
dom = linspace(sol.a,sol.b,num_pts);
[f,f_der] = sol.deval(dom);

if isempty(sol.params)
    for j = 1:length(dom)
        ode_err = max(ode_err,max(abs(f_der(:,j)-sol.ode(dom(j),f(:,j)))));
    end   
    bc_err = max(abs(sol.bc(sol.deval)));
else
    for j = 1:length(dom)
        ode_err = max(ode_err,max(abs(f_der(:,j)-sol.ode(dom(j),f(:,j),sol.params))));
    end
    bc_err = max(abs(sol.bc(sol.deval,sol.params)));
end



 
% -------------------------------------------------------------------------
% eval 
% -------------------------------------------------------------------------

function [f,f_der] = local_deval(dom,sol)

% evaluate the profile and its derivative
dom_tilde = (2/(sol.a-sol.b))*(dom-0.5*(sol.a+sol.b));
theta_dom = acos(dom_tilde);
ind = 0:1:sol.N-1;
T = cos(theta_dom.'*ind);
f = (T*sol.cf).';
f_der = (T*sol.der_cf).';


% -------------------------------------------------------------------------
% ode_fun 
% -------------------------------------------------------------------------

function out = ode_fun(u,ode,bc,sol)
% ODE for fsolve

sol.y = real(reshape(u(1:sol.dim*sol.N),sol.N,sol.dim).');
sol.params = real(u(sol.dim*sol.N+1:end));
sol.cf = sol.Tcf*sol.y.';
sol.der_cf = sol.Tcf*(sol.Tx*sol.y.');
current_deval = @(x)local_deval(x,sol);

y = reshape(u(1:sol.dim*sol.N),sol.N,sol.dim).';
params = u(sol.dim*sol.N+1:end);

if isempty(params)
    %current_deval
    out_bc = bc(current_deval);
else
    out_bc = bc(current_deval,params);
end
out = zeros(length(out_bc)+sol.dim*sol.N,1);

ders = zeros(sol.dim,sol.N);
for j = 1:sol.N
    if isempty(params)
        ders(:,j) = ode(sol.x(j),y(:,j));
    else
        ders(:,j) = ode(sol.x(j),y(:,j),params);
    end
end

for j = 1:sol.dim
    out((j-1)*sol.N+1:j*sol.N) = sol.Tx*y(j,:).'-ders(j,:).';
end

out(sol.dim*sol.N+1:end) = out_bc;




