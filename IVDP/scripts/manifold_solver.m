clc; close all; clear all;

%
% system parameters
%

p.a = 0.5;
p.eps = 0;
p.sigma1 = 1;
p.sgima2 = 1;

%
% directory management
%

curr_dir = cd; cd('../code');

% Get the coefficients for the unstable manifold at negative infinity
N = 40; % The number of terms for the manifold
scale = 6; % Scale the eigenvectors in the parameterization method
UM = manifold_infty_coeffs(p,N,scale);

cd('../data');
save(replace(['manifold_a_',num2str(p.a),'_eps_',num2str(p.eps)],'.','P'),'UM');
cd('../scripts');