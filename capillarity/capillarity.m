beep off; clc; clear all; close all;

% INFO

% Check out 
% http://math.byu.edu/~jeffh/publications/papers/cap.pdf
% for details about this system


%
% parameters
%

p.gamma=1.4;
p.vp=.15;
% this example solves the profile for $d = -0.45$ and then uses continuation thereafter
d_vals = -(0.45:0.001:0.451); 

%
% dependent parameters
%

p.a=-(1-p.vp)/(1-p.vp^-p.gamma);

% solve profile. Use continuation as an example
for j = 1:length(d_vals)

    p.d = d_vals(j);
    disp(p);
    
    %
    % profile
    %

    s.n = 2; % this is the dimension of the profile ode
    % we divide the domain in half to deal with the 
    % non-uniqueness caused by translational invariance
    % s.side = 1 means we are solving the profile on the interval [0,X]
    s.side=1; 
    s.F=@F; % F is the profile ode
    s.Flinear = @J; % J is the profile ode Jacobian
    s.UL = [1;0]; % These are the endstates of the profile and its derivative at x = -infty
    s.UR = [p.vp;0]; % These are the endstates of the profile and its derivative at x = +infty
    s.phase = 0.5*(s.UL+s.UR); % this is the phase condition for the profile at x = 0
    s.order = 1; % this indicates to which componenet the phase conditions is applied
    s.stats = 'on'; % this prints data and plots the profile as it is solved

    if j == 1
        % there are some other options you specify. You can look in profile_flux to see them
        [p,s] = profile_flux(p,s); % solve profile for first time
        s_old = s;
    else
        % this time we are using continuation
        [p,s] = profile_flux(p,s,s_old); % solve profile
    end
    
end


% plot final solution
x = linspace(s.L,s.R,200);
y = zeros(2,length(x));
for j = 1:length(x)
    y(:,j) = soln(x(j),s);
end
plot(x,y)
drawnow;


%
% structure variables
%

% Here you can choose the method you use for the Evans function, or you can set the option
% to default and let it choose. [2,2] is the size of the manifold evolving from - and + infy in the 
% Evans solver. 'front' indicates you are solving for a traveling wave front and not a periodic solution
% [s,e,m,c] = emcset(s,'front',LdimRdim(@A,s,p),'default'); % default for capillarity is reg_reg_polar
% [s,e,m,c] = emcset(s,'front',[2,2],'reg_adj_polar');
% [s,e,m,c] = emcset(s,'front',[2,2],'adj_reg_polar');
% [s,e,m,c] = emcset(s,'front',[2,2],'reg_reg_polar');

% This choice solves the right hand side via exterior products
[s,e,m,c] = emcset(s,'front',[2,2],'adj_reg_compound',@A,@Ak);

% display a waitbar
c.stats = 'print'; % 'on', 'print', or 'off' (parallel matlab does not work with 'on' option)
c.ksteps = 2^8;

% Here you can choose differnt ode solvers to evolve the Evans function manifolds you can see
% which works best for your system

% m.ode_fun = @ode45;

m.ode_fun = @ode15s;

% m.options = odeset('RelTol',1e-6,'AbsTol',1e-8,'Refine',1,'Stats','off','MaxOrder',2);
% m.ode_fun = @ode15s;
 


%
% preimage contour
%

% This is a semi circle. You can also do a semi annulus or a rectangle (look in bin_main)
circpnts=30; imagpnts=30; R=10; spread=2; zerodist=10^(-2);
preimage=semicirc(circpnts,imagpnts,c.ksteps,R,spread,zerodist);

%
% compute Evans function
%

tstart = tic;
halfw=contour(c,s,p,m,e,preimage);
tstop = toc(tstart);
fprintf('\nRun time = %4.4g seconds.\n',tstop);
w = halfw/halfw(1);
w = [w fliplr(conj(w))]; % We compute the Evans function on half of contour then reflect

% 
% process and display data
%

wnd=winding_number(w); % determine the number of roots inside the contour
fprintf('Winding Number: %1.1d\n',wnd);

% plot the Evans function (normalized)
plot_evans(w)



