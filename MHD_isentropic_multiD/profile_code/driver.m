clc; clear all; 

% add DMSUITE to path
try
    eval(['addpath ',pwd,'/../../../dmsuite']);
catch me
    error(['This code requires DMSUITE by J.A.C Weideman, ', ...
        'which is free from the MATLAB File Exchange. Place ', ...
        'DMSUITE at the same file level as stablab_matlab.']);
end

% Independent parameters
p.gamma = 5/3;
p.u1_p = 0.9;

phi = 0.1;
p.nu = phi;
p.mu = phi;
p.eta = -(2/3)*p.mu;
p.kappa = phi;

% solve for the profile
end_point_tol = 1e-8;
[s,p] = get_profile(p,end_point_tol);

% plot the profile
x = linspace(s.L,s.R,1000);
y = s.sol.deval(x);
plot(x,y,'-k','LineWidth',2);

