function out = A_local_model_balflux(x,lambda,s,p)
% Balanced flux Evans matrix for the system named local_model.

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

ynot = [tau,S,];

temp2 = F_local_model(x,ynot,s,p);

tau_x = temp2(1);
S_x = temp2(2);

v = v_neg - spd * (tau - tau_neg);
v_x = -spd * tau_x;
T_x = S_x * exp(S) / tau - exp(S) / tau ^ 2 * tau_x;


% Evans matrix
out = [lambda / spd 0 0 -0.1e1 / spd 0; 0 0 0 1 0; lambda * (-exp(S) / tau ^ 2 + tau - (exp(S) / tau + 0.1e1) / tau) / spd 0 0 -(-exp(S) / tau ^ 2 + tau - (exp(S) / tau + 0.1e1) / tau) / spd + v (exp(S) / tau + 0.1e1) / exp(S) * tau; lambda * tau / mu * (v_x * mu / tau ^ 2 - 0.3e1 * exp(S) / tau ^ 3 - 0.1e1) / spd lambda * tau / mu 0 -tau / mu * ((v_x * mu / tau ^ 2 - 0.3e1 * exp(S) / tau ^ 3 - 0.1e1) / spd + spd) 0.1e1 / mu; lambda * (-v * tau / kappa * (v_x * mu / tau ^ 2 - 0.3e1 * exp(S) / tau ^ 3 - 0.1e1) + tau / kappa * (-spd * (-exp(S) / tau ^ 2 + tau - (exp(S) / tau + 0.1e1) / tau) + v_x * mu * v / tau ^ 2 + T_x * kappa / tau ^ 2 + v * (-0.3e1 * exp(S) / tau ^ 3 - 0.1e1))) / spd -lambda * v * tau / kappa lambda * tau / kappa v * tau / kappa * ((v_x * mu / tau ^ 2 - 0.3e1 * exp(S) / tau ^ 3 - 0.1e1) / spd + spd) - tau / kappa * ((-spd * (-exp(S) / tau ^ 2 + tau - (exp(S) / tau + 0.1e1) / tau) + v_x * mu * v / tau ^ 2 + T_x * kappa / tau ^ 2 + v * (-0.3e1 * exp(S) / tau ^ 3 - 0.1e1)) / spd + spd * v + v_x * mu / tau - exp(S) / tau ^ 2 + tau) -v / kappa - tau / kappa * (spd * (exp(S) / tau + 0.1e1) / exp(S) * tau - v / tau);];


