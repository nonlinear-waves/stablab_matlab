
%{
INFO

Check out the following article for more details about
this system.

Ghazaryan, Anna, Stephane Lafortune, and Vahagn Manukian. 
"Spectral Analysis of Fronts in a Marangoni-Driven Thin 
Liquid Film Flow Down a Slope." SIAM Journal on Applied 
Mathematics 80, no. 1 (2020): 95-118.

%}

beep off; clc; clear all; close all;

%
% parameters
%

p.hL=4;
p.hR=1;
p.D=1;
p.be=.2;
p.C=.01;
p.maxG=2;

%
% dependent parameters
%

p.ss=1/3*(p.hL^2+p.hL*p.hR+p.hR^2);
p.K1=-1/3*p.hL*p.hR*(p.hL+p.hR);

%
% solve for the profile
%

s.n = 4; % this is the dimension of the profile ode
% we divide the domain in half to deal with the 
% non-uniqueness caused by translational invariance
% s.side = 1 means we are solving the profile on the interval [0,X]
s.side=1; 
s.F=@F; % F is the profile ode
s.Flinear = @J; % J is the profile ode Jacobian
s.UL = [p.hL;0;0;0]; % These are the endstates of the profile and its derivative at x = -infty
s.UR = [p.hR;0;0;0]; % These are the endstates of the profile and its derivative at x = +infty
s.phase = [-3*p.K1/p.ss; 0 ;0; p.maxG]; % this is the phase condition for the profile at x = 0
s.order = [1; 4]; % this indicates to which componenent the phase conditions is applied
s.stats = 'on'; % this prints data and plots the profile as it is solved

s.guess = @(x)profile_guess(x,s,p); % guess for the profile solution

% there are some other options you specify. You can look in profile_flux to see them
[p,s] = profile_flux(p,s); % solve profile for first time

%
% plot the profile
%

plot_profile = 'on';
if strcmp(plot_profile,'on')
    % plot final solution
    figure;
    x = linspace(s.L,s.R,1001);
    y = zeros(4,length(x));
    for j = 1:length(x)
        y(:,j) = real(soln(x(j),s));
    end
    plot(x,y,'LineWidth',2);
    h = xlabel('x');
    set(h,'FontSize',22);
    h2 = gca;
    set(h2,'FontSize',22);
    drawnow;
end


%
% structure variables
%

[s,e,m,c] = emcset(s,'front',[3,3],'reg_reg_polar');

% display a waitbar
c.stats = 'print'; % 'on', 'print', or 'off' (parallel matlab does not work with 'on' option)
c.ksteps = 2^8;
c.lambda_steps=100;

%
% preimage contour - semi-circle
%

% This is a semi circle. You can also do a semi annulus or a rectangle (look in bin_main)
circpnts=30; imagpnts=30; R=1; spread=2; zerodist = 1e-4;
preimage=semicirc(circpnts,imagpnts,c.ksteps,R,spread,zerodist,c.lambda_steps);

%
% compute Evans function
%
c.refine = 'on';
c.tol = 0.1;

tstart = tic;
[halfw pp2]=contour(c,s,p,m,e,preimage);
w = halfw/halfw(1); % We normalize the Evans function.
w = [w fliplr(conj(w))]; % We compute the Evans function on half of contour then reflect
tstop = toc(tstart);

fprintf('\nRun time = %4.4g seconds.\n',tstop);

% 
% process and display data
%

wnd=winding_number(w); % determine the number of roots inside the contour
fprintf('Winding Number: %1.1d\n',wnd);

% plot the Evans function (normalized)
plot_evans(w)






