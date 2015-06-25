beep off; clc; clear all; close all;

%
% parameters
%

p.gamma=5/3;
p.vp = 1e-5;

%
% dependpent variables
% 

p.a=-(1-p.vp)/(1-p.vp^-p.gamma);

%
% profile
%

s.n = 1; % this is the dimension of the profile ode
% we divide the domain in half to deal with the 
% non-uniqueness caused by translational invariance
% s.side = 1 means we are solving the profile on the interval [0,X]
s.side=1; 
s.F=@F; % F is the profile ode
s.Flinear = @Flinear; % Flinear is the profile ode Jacobian
s.UL =1; % These are the endstates of the profile and its derivative at x = -infty
s.UR = p.vp; % These are the endstates of the profile and its derivative at x = +infty
s.phase = 0.5*(s.UL+s.UR); % this is the phase condition for the profile at x = 0
s.order = 1; % this indicates to which componenet the phase conditions is applied
s.stats = 'on'; % this prints data and plots the profile as it is solved
[p,s] = profile_flux(p,s);

evan_mat = @A;
% [s,e,m,c] = emcset(s,'front',[1,2],'default',evan_mat); % default for isentropic is reg_adj_polar
% [s,e,m,c] = emcset(s,'front',[1,2],'reg_reg_polar',evan_mat);
% [s,e,m,c] = emcset(s,'front',[1,2],'adj_reg_polar',evan_mat);
[s,e,m,c] = emcset(s,'front',[1,2],'reg_adj_polar',evan_mat);

% display a waitbar
c.stats = 'print';

R=(sqrt(p.gamma)+0.5)^2;
circpnts=20; imagpnts=20; spread=4; zerodist=1e-4;
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
