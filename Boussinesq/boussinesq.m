% boussinesq driver

clc; clear all; close all; beep off;

% parameters
p.S = 0.4;

%
% profile
%

s.I = 8;
s.R=s.I;
s.L=-s.I;

%
% STABLAB structures
%

% [s,e,m,c] = emcset(s,'front',[2,2],'default'); % default for boussinesq is reg_reg_polar
% [s,e,m,c] = emcset(s,'front',[2,2],'reg_adj_polar');
% [s,e,m,c] = emcset(s,'front',[2,2],'adj_reg_polar');
% [s,e,m,c] = emcset(s,'front',[2,2],'reg_reg_polar');
[s,e,m,c] = emcset(s,'front',[2,2],'reg_adj_compound'); 
% [s,e,m,c] = emcset(s,'front',[2,2],'adj_reg_compound');

% display a waitbar
c.stats = 'print';

%
% preimage
%

points = 50;
preimage = 0.16+0.05*exp(2*pi*1i*linspace(0,0.5,points+(points-1)*c.ksteps));

% 
% compute Evans function
%

halfw=contour(c,s,p,m,e,preimage);
w = [halfw(1:end-1) fliplr(conj(halfw))];

% 
% process and display data
%

wnd=winding_number(w);
plot(w/w(1),'.-k');
fprintf('Winding Number: %1.1d\n',wnd);


