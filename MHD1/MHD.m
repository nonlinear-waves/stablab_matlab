clc; beep off; clear all; close all;

% -------------------------------------------------------------------------
% parameters
% -------------------------------------------------------------------------

p.B=2;
p.gamma=5/3;
p.mu0=1;
p.sigma=1;
p.vp=0.0001;
p.mu=1;
p.eta = -2*p.mu/3;

% -------------------------------------------------------------------------
% dependent parameters
% -------------------------------------------------------------------------

p.a = p.vp^p.gamma*((1-p.vp)/(1-p.vp^p.gamma));

% -------------------------------------------------------------------------
% profile solution
% -------------------------------------------------------------------------

s.F = @(x,y,s,p)((2*p.mu+p.eta)^(-1)*y*(y-1+p.a*(y^(-p.gamma)-1)));
s.Flinear = @(y,p)((2*p.mu+p.eta)^(-1)*(y-1+p.a*(y^(-p.gamma)-1)) ... 
    +(2*p.mu+p.eta)^(-1)*y*(1-p.a*p.gamma*y^(-p.gamma-1)));
%number of profile equations to integrate
s.n = 1;
s.order = 1;
s.phase = 0.5*(1+p.vp);
%end states
s.UL = 1;
s.UR = p.vp;

s.stats = 'on';
%tolerance at end states
s.tol = 1e-6;

[p,s] = profile_flux(p,s);


%
% structure variables
%

[s,e,m,c] = emcset(s,'front',[2,2],'default'); % default for MHD1 is reg_reg_polar
% [s,e,m,c] = emcset(s,'front',[2,2],'reg_adj_polar');
% [s,e,m,c] = emcset(s,'front',[2,2],'adj_reg_polar');
% [s,e,m,c] = emcset(s,'front',[2,2],'reg_reg_polar');
% [s,e,m,c] = emcset(s,'front',[2,2],'adj_reg_compound');
% [s,e,m,c] = emcset(s,'front',[2,2],'reg_adj_compound');

% display a waitbar
c.stats = 'off';
 
m.ode_fun = @ode15s;

%
% preimage contour
%

R=1;
circpnts=30; imagpnts=30; spread=4; zerodist=10^(-4); 
preimage=semicirc(circpnts,imagpnts,c.ksteps,R,spread,zerodist);

%
% compute Evans function 
%

halfw=contour(c,s,p,m,e,preimage);
w = [halfw fliplr(conj(halfw))];

% 
% process and display data
%

wnd=winding_number(w);
fprintf('Winding Number: %1.1d\n',wnd);
plot_evans(w);
