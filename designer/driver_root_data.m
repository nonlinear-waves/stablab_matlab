% driver_root_data
%
% Find roots for the designer system,
%
%   u_t + (A(v)u)_x = (B(v)u_x)_x
%   v_t + ((1/2)u \cdot A'(v)u + (1/2)v^2)_x = v_{xx}
%
%   (\bar u, \bar v) = (0, \bar v(x))
%
%   A(v) = R_{\theta(v)}A_m R_{-\theta(v)}
%   B(v) = Id
%
%   R_{\theta} = \begin{pmatrix}{ \cos(\theta) & - \sin(\theta)) \\ \sin(\theta)&
%   \cos(\theta) \end{pmatrix}
%
%   Am = [1 0; 0 -1];
%  
%   \theta(v) = M \pi v
%
%   \bar v_{\gamma}(x) = -\gamma tanh(\gamma x/2)

beep off; clc; clear all; close all;

%
% parameters
%

p.M = 3.2836;
p.gamma = 0.66;

% Evans function matrix
Amat = @A_coord_change;

s.I = 10/p.gamma;
s.R = s.I;
s.L = -s.I;

% stablab structures
[s,e,m,c] = emcset(s,'front',[2,2],'reg_reg_polar',Amat);

% Evans options
m.options = odeset('RelTol',1e-6,'AbsTol',1e-8,'Refine',1,'Stats','off');

c.tol = 0.2; 
c.lambda_steps = 0;
% refine the Evans function computation to achieve set relative error
c.refine = 'on';
c.best_refine = 'on';
% display a waitbar
c.stats = 'off';
% display squares during root solving
c.pic_stats = 'on';
% Root solving
c.lambda_steps = 0;
c.stats = 'off';
c.ksteps = 2^8;
c.moments = 'off';
c.tol = 0.2;

c.root_fun = @contour;
tol = 0.05;
rt_ep = 1e-4;
R = 0.11;
box = [1e-4,rt_ep,R,R];
rts1 = root_solver1(box,tol,p,s,e,m,c);

box = [1e-4,-1e-4,R];
rts2 = root_solver2(box,tol,p,s,e,m,c);

rts = [rts1,conj(rts1),rts2];

plot(conj(rts1),'.k','MarkerSize',18);

disp(rts)


