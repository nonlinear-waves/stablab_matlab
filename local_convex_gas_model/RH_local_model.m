function p = RH_local_model(p)
% Rankine Hugoniot conditions


st.H_fun = @H_local_model;

    st.left_H = 1e-6;
    st.right_H = 1000;
    st.H_ind_var = 's';
    p.tau_neg = full_gas_Hugoniot(p,st,p.S_neg);


p.S_plus = p.S0;
p.tau_plus = p.tau0;

none = p.none;
mu = p.mu;
kappa = p.kappa;


S = p.S_neg;
tau = p.tau_neg;
press_neg = exp(S) / tau ^ 2 - tau;
T_neg = exp(S) / tau + 0.1e1;
p.T_neg = T_neg;

S = p.S_plus;
tau = p.tau_plus;
press_plus = exp(S) / tau ^ 2 - tau;
T_plus = exp(S) / tau + 0.1e1;
p.T_plus = T_plus;

p.spd = -sqrt((press_neg-press_plus)/(p.tau_plus-p.tau_neg));

p.m1 = 0;

p.v_plus = -p.spd*p.tau_plus-p.m1;
p.v_neg = -p.spd*p.tau_neg-p.m1;
