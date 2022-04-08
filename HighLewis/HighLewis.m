beep off; clc; clear all; close all;

% INFO

% Computation of the Evans function for the High Lewis Number Combustion
% model as present in the book "An introduction to traveling waves: fronts,
% pulses, and wavetrains" written by
% Ghazaryan, Anna, Stephane Lafortune, and Vahagn Manukian.


%
% parameters
%

addpath('WaveComputation');
p.be=1;
p.c = integrated_find_c(p.be);





% solve profile. 

    s.n = 2; % this is the dimension of the profile ode
    % we divide the domain in half to deal with the 
    % non-uniqueness caused by translational invariance
    % s.side = 1 means we are solving the profile on the interval [0,X]
    s.side=1; 
    s.F=@F; % F is the profile ode
    s.Flinear = @J; % J is the profile ode Jacobian
    s.UL = [1/p.be;0]; % These are the endstates of the profile and its derivative at x = -infty
    s.UR = [10^(-9);1]; % These are the endstates of the profile and its derivative at x = +infty
    s.phase = [1/2/p.be; 1/2]; % this is the phase condition for the profile at x = 0
    s.order = [2]; % this indicates to which componenent the phase conditions is applied
    s.stats = 'on'; % this prints data and plots the profile as it is solved

   
        % there are some other options you specify. You can look in profile_flux to see them
        [p,s] = profile_flux(p,s); % solve profile for first time
        s_old = s;
    


%plot final solution
figure(66);
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
[s,e,m,c] = emcset(s,'front',[2,1],'reg_reg_polar');

% This choice solves the right hand side via exterior products
% [s,e,m,c] = emcset(s,'front',[3,3],'adj_reg_compound',@A,@Ak);

% display a waitbar
c.stats = 'print'; % 'on', 'print', or 'off' (parallel matlab does not work with 'on' option)
%c.ksteps = 2^8;
%c.lambda_steps=100;

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
circpnts=100; imagpnts=100; r=0.1; spread=2; zerodist=0;
preimage=semicirc(circpnts,imagpnts,c.ksteps,r,spread,zerodist,c.lambda_steps)-0.0001;

%circpnts=100; imagpnts=100; innerpnts = 50; r=1; 
%spread=4; inner_radius=0.01;
%preimage=semicirc2(circpnts,imagpnts,innerpnts,c.ksteps,r,spread,inner_radius,c.lambda_steps);
%The countour to be used to pinpoint the eigenvalue at lambda=0.1765i when beta = 6.572
%preimage = 0.1765*1i+0.05*exp(2*pi*1i*linspace(0,1,points+(points-1)*c.ksteps));
%preimage = 0.1765*1i-0.1*1i+1i*linspace(0,.5,points+(points-1)*c.ksteps);
%preimage = 0.16+0.05*exp(2*pi*1i*linspace(0,1,points+(points-1)*c.ksteps));

%Below, conditions are checked to make sure the Evans function stays in a region where it is analytic  
for cc = 1:length(preimage)
    ll=preimage(cc);
    if or(real(-p.c/2-1/2*sqrt(p.c^2+4*ll))>real(ll)/p.c,real(ll)<-p.c^2/4)
        disp('ALERT at + infinity')
    end;
    if or(real((ll+p.be*exp(-p.be))/p.c)<real(-p.c/2-1/2*sqrt(p.c^2+4*ll)),real(ll)<-p.c^2/4)
        disp('ALERT at - infinity')
    end;
end


figure(99);
plot(preimage);

%
% compute Evans function
%
c.refine = 'on'
c.tol=0.1;

tstart = tic;
[halfw pp2]=contour(c,s,p,m,e,preimage);



tstop = toc(tstart);
fprintf('\nRun time = %4.4g seconds.\n',tstop);
w = halfw/halfw(1);
%w=halfw;

w = [w fliplr(conj(w))]; % We compute the Evans function on half of contour then reflect


% 
% process and display data
%

wnd=winding_number(w); % determine the number of roots inside the contour
fprintf('Winding Number: %1.1d\n',wnd);

% plot the Evans function (normalized)
plot_evans(w)


figure(77);
plot(w);





