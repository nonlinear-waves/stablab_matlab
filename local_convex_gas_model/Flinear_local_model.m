function out = Flinear_local_model(y,p)


tau = y(1);
S = y(2);

% wave speed
spd = p.spd;

% values at - infinity
v_neg = p.v_neg;
tau_neg = p.tau_neg;
S_neg = p.S_neg;

% parameters
m1 = p.m1;
none = p.none;
mu = p.mu;
kappa = p.kappa;

%Jacobian
out = [-0.1e1 / spd / mu * (spd ^ 2 * (tau - tau_neg) + exp(S) / tau ^ 2 - tau - exp(S_neg) / tau_neg ^ 2 + tau_neg) - tau / spd / mu * (spd ^ 2 - 0.2e1 * exp(S) / tau ^ 3 - 0.1e1) -0.1e1 / tau / spd / mu * exp(S); (-0.1e1 / tau / spd / mu - 0.1e1 / exp(S) * tau * (v_neg - spd * (tau - tau_neg)) / kappa) * (spd ^ 2 * (tau - tau_neg) + exp(S) / tau ^ 2 - tau - exp(S_neg) / tau_neg ^ 2 + tau_neg) + 0.1e1 / exp(S) * tau / kappa * (-spd * (exp(S) / tau + S + tau ^ 2 / 0.2e1 + (v_neg - spd * (tau - tau_neg)) ^ 2 / 0.2e1 - exp(S_neg) / tau_neg - S_neg - tau_neg ^ 2 / 0.2e1 - v_neg ^ 2 / 0.2e1) + (v_neg - spd * (tau - tau_neg)) * (exp(S) / tau ^ 2 - tau) - v_neg * (exp(S_neg) / tau_neg ^ 2 - tau_neg)) + tau * ((0.1e1 / tau ^ 2 / spd / mu - 0.1e1 / exp(S) * (v_neg - spd * (tau - tau_neg)) / kappa + 0.1e1 / exp(S) * tau * spd / kappa) * (spd ^ 2 * (tau - tau_neg) + exp(S) / tau ^ 2 - tau - exp(S_neg) / tau_neg ^ 2 + tau_neg) + (-0.1e1 / tau / spd / mu - 0.1e1 / exp(S) * tau * (v_neg - spd * (tau - tau_neg)) / kappa) * (spd ^ 2 - 0.2e1 * exp(S) / tau ^ 3 - 0.1e1) + 0.1e1 / exp(S) / kappa * (-spd * (exp(S) / tau + S + tau ^ 2 / 0.2e1 + (v_neg - spd * (tau - tau_neg)) ^ 2 / 0.2e1 - exp(S_neg) / tau_neg - S_neg - tau_neg ^ 2 / 0.2e1 - v_neg ^ 2 / 0.2e1) + (v_neg - spd * (tau - tau_neg)) * (exp(S) / tau ^ 2 - tau) - v_neg * (exp(S_neg) / tau_neg ^ 2 - tau_neg)) + 0.1e1 / exp(S) * tau / kappa * (-spd * (-exp(S) / tau ^ 2 + tau - (v_neg - spd * (tau - tau_neg)) * spd) - spd * (exp(S) / tau ^ 2 - tau) + (v_neg - spd * (tau - tau_neg)) * (-0.2e1 * exp(S) / tau ^ 3 - 0.1e1))) tau * (0.1e1 / exp(S) * tau * (v_neg - spd * (tau - tau_neg)) / kappa * (spd ^ 2 * (tau - tau_neg) + exp(S) / tau ^ 2 - tau - exp(S_neg) / tau_neg ^ 2 + tau_neg) + (-0.1e1 / tau / spd / mu - 0.1e1 / exp(S) * tau * (v_neg - spd * (tau - tau_neg)) / kappa) * exp(S) / tau ^ 2 - 0.1e1 / exp(S) * tau / kappa * (-spd * (exp(S) / tau + S + tau ^ 2 / 0.2e1 + (v_neg - spd * (tau - tau_neg)) ^ 2 / 0.2e1 - exp(S_neg) / tau_neg - S_neg - tau_neg ^ 2 / 0.2e1 - v_neg ^ 2 / 0.2e1) + (v_neg - spd * (tau - tau_neg)) * (exp(S) / tau ^ 2 - tau) - v_neg * (exp(S_neg) / tau_neg ^ 2 - tau_neg)) + 0.1e1 / exp(S) * tau / kappa * (-spd * (exp(S) / tau + 0.1e1) + (v_neg - spd * (tau - tau_neg)) * exp(S) / tau ^ 2));];
