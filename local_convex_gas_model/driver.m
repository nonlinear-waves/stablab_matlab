% driver for balanced flux formulation

% for full gas in Lagrangian coordinates
beep off; clc; clear all; close all;

% parameters

p.S_neg = - 1;

p.none = 1;
p.mu = 1;
p.kappa = 1;
p.S0 = 0;
p.tau0 = 1;

% plus and minus infinity end states
 p = RH_local_model(p);

% phase condition
s.phase = [0.5*(p.tau_plus+p.tau_neg),0.5*(p.S_plus+p.S_neg),];

%order in which to apply phase conditions
s.order = [1,2];

%profile ode
s.F = @F_local_model;

%Jacobian file
s.Flinear = @Flinear_local_model;

%number of profile equations to integrate
s.n = 2;

%end states
s.UL = [p.tau_neg; p.S_neg ];
s.UR = [p.tau_plus; p.S_plus ];

s.stats = 'on';
%tolerance at end states
s.tol = 1e-6;

[p,s] = profile_flux(p,s);
s_old = s;

S_neg_vals = linspace(-1,-5,30);
for j = 2:length(S_neg_vals)
    clc;
    disp(j);
    p.S_neg = S_neg_vals(j);
    % plus and minus infinity end states
    p = RH_local_model(p);
    % phase condition
    s.phase = [0.5*(p.tau_plus+p.tau_neg),0.5*(p.S_plus+p.S_neg),];
    s.UL = [p.tau_neg; p.S_neg; ];
    s.UR = [p.tau_plus; p.S_plus; ];
    [p,s] = profile_flux(p,s,s_old);
    s_old = s;
end

plot_profile_local_model(p,s);

% Evans matrix
Amat = @A_local_model_balflux;

%
% structure variables
%

[s,e,m,c] = emcset(s,'front',LdimRdim(Amat,s,p),'default',Amat);

% refine the Evans function computation to achieve set relative error
c.refine = 'off';

% display a waitbar
c.stats = 'print';

% use drury_no_radil; MU is too large in the method of continuous
% orthogonalization
m.method = @drury_no_radial;

c.max_R = 1000;
c.Rtol = 0.2;
R = radius(c,s,p,m,e);

circpnts = 10; imagpnts = 30; spread = 4; zerodist = 10^(-4);
preimage = semicirc(circpnts,imagpnts,c.ksteps,R,spread,zerodist);

% 
% compute Evans function
%

halfw = contour(c,s,p,m,e,preimage);
w = [halfw fliplr(conj(halfw))];

%
% plot the Evans function
%

plot_evans(w);

wnd = winding_number(w);

fprintf('\n\nWinding Number: %4.4g\n\n',wnd);






