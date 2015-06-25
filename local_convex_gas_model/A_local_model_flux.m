function out = A_local_model_flux(x,lambda,s,p)
% Unbalanced flux Evans matrix for the system named local_model.

% wave speed
spd = p.spd;

% parameters
v_neg = p.v_neg;
tau_neg = p.tau_neg;
S_neg = p.S_neg;
m1 = p.m1;
none = p.none;
mu = p.mu;
kappa = p.kappa;


% interpolate profile solution
temp = soln(x,s);

%profile values
tau = temp(1);
S = temp(2);

v = v_neg - spd * (tau - tau_neg);
v_x = -spd * tau_x;
T_x = S_x * exp(S) / tau - exp(S) / tau ^ 2 * tau_x;


% Evans matrix
out = [lambda / spd 0 0 -lambda / spd 0; 0 0 0 lambda 0; lambda * (-exp(S) / tau ^ 2 + tau - (exp(S) / tau + 0.1e1) / tau) / spd 0 0 lambda * (-(-exp(S) / tau ^ 2 + tau - (exp(S) / tau + 0.1e1) / tau) / spd + v) lambda * (exp(S) / tau + 0.1e1) / exp(S) * tau; tau / mu * (v_x * mu / tau ^ 2 - 0.3e1 * exp(S) / tau ^ 3 - 0.1e1) / spd tau / mu 0 -tau / mu * ((v_x * mu / tau ^ 2 - 0.3e1 * exp(S) / tau ^ 3 - 0.1e1) / spd + spd) 0.1e1 / mu; (-v * tau / kappa * (v_x * mu / tau ^ 2 - 0.3e1 * exp(S) / tau ^ 3 - 0.1e1) + tau / kappa * (-spd * (-exp(S) / tau ^ 2 + tau - (exp(S) / tau + 0.1e1) / tau) + v_x * mu * v / tau ^ 2 + T_x * kappa / tau ^ 2 + v * (-0.3e1 * exp(S) / tau ^ 3 - 0.1e1))) / spd -v * tau / kappa tau / kappa v * tau / kappa * ((v_x * mu / tau ^ 2 - 0.3e1 * exp(S) / tau ^ 3 - 0.1e1) / spd + spd) - tau / kappa * ((-spd * (-exp(S) / tau ^ 2 + tau - (exp(S) / tau + 0.1e1) / tau) + v_x * mu * v / tau ^ 2 + T_x * kappa / tau ^ 2 + v * (-0.3e1 * exp(S) / tau ^ 3 - 0.1e1)) / spd + spd * v + v_x * mu / tau - exp(S) / tau ^ 2 + tau) -v / kappa - tau / kappa * (spd * (exp(S) / tau + 0.1e1) / exp(S) * tau - v / tau);];


