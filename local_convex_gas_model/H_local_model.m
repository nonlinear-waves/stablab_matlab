function out = H_local_model(p,tau,S)
% Hugoniot curve equation


S0 = p.S0;
tau0 = p.tau0;
none = p.none;
mu = p.mu;
kappa = p.kappa;


out = exp(S) / tau + S + tau ^ 2 / 0.2e1 - exp(S0) / tau0 - S0 - tau0 ^ 2 / 0.2e1 + (exp(S) / tau ^ 2 - tau + exp(S0) / tau0 ^ 2 - tau0) * (tau - tau0) / 0.2e1;
