% Burgers
beep off; clc; clear all; close all

tstart = tic;

%
% parameters
%

p.ul=1;
p.ur=0;
p.integrated = 'off';
 
%
% numerical infinity
%

s.I=12;
s.R=s.I;
s.L=-s.I;

%
% set STABLAB structures to local default values
%

[s,e,m,c] = emcset(s,'front',[1,1],'default'); % default for Burgers is reg_reg_polar

%
% preimage contour
%

circpnts=20; imagpnts=20; innerpnts = 5; r=10; spread=4; zerodist=10^(-2);
preimage=semicirc2(circpnts,imagpnts,innerpnts,c.ksteps,r,spread,zerodist,c.lambda_steps);

%
% compute the Evans function
%

halfw=contour(c,s,p,m,e,preimage);
w = [halfw fliplr(conj(halfw))];
wnd=winding_number(w);

%
% display Evans function output and statistics
%

plot(w,'-*k');
fprintf('Time: %1.2d seconds\n',toc(tstart));
fprintf('Winding Number: %1.1d\n',wnd);

