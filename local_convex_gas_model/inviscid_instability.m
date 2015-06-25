clc; close all; beep off; clear all;

%
% parametrs
%

p.S0 = 0; 
p.tau0 = 1;

S_vals = linspace(p.S0,-10,1000);

% 
% program
%

p.model = ' $\\bar e(s,\\tau) = e^s/\\tau+S+\\tau^2/2$ ';

%
% instantiate variables
%

ceros = zeros(length(S_vals),1);
tau_vals = ceros;
e = ceros;
press = ceros;
e_s = ceros;
e_tau_tau = ceros;
neg_e_s_tau = ceros;
neg_e_tau_tau_tau = ceros;
e_s_s = ceros;
e_s_tau = ceros;
conv = ceros;
delta = ceros;
H = ceros;
sigma = ceros;
v = ceros;

p.none = 0;
p.mu = 1;
p.kappa = 1;

S0 = p.S0;
tau0 = p.tau0;
for j = 1:length(S_vals)
    
    S = S_vals(j);    
    st.H_fun = @H_local_model;
    st.left_H = 1e-6;
    st.right_H = 1000;
    st.H_ind_var = 's';
    tau = full_gas_Hugoniot(p,st,S);
    
    tau_vals(j) = tau;
    
    delta(j) = (exp(S0) - tau0 ^ 3) / tau0 ^ 2 - exp(S) / tau ^ 2 + tau - (0.1e1 + sqrt(((exp(S0) - tau0 ^ 3) / tau0 ^ 2 - exp(S) / tau ^ 2 + tau) / (tau0 - tau) / (-0.2e1 * exp(S0) / tau0 ^ 3 - 0.1e1))) * (exp(S0) / tau0 + 0.1e1) * (0.2e1 * exp(S0) / tau0 ^ 3 + 0.1e1) / exp(S0) * tau0 ^ 2;
 
    e(j) = exp(S) / tau + S + tau ^ 2 / 0.2e1;
    press(j) =  exp(S) / tau ^ 2 - tau;
    sigma(j) =  -sqrt(-(exp(S) / tau ^ 2 - tau - (exp(S0) - tau0 ^ 3) / tau0 ^ 2) / (tau - tau0));
    
    v(j) = sqrt(-(exp(S) / tau ^ 2 - tau - (exp(S0) - tau0 ^ 3) / tau0 ^ 2) / (tau - tau0)) * (tau - tau0);
    
end

for j = 1:length(delta)-1
    if delta(j)*delta(j+1) < 0
        taustar = tau_vals(j);
        Sstar = S_vals(j);
    end
end

%------------------------------------------------------------
% plot results
%------------------------------------------------------------

% use to plot x-axis
xax = [S_vals(1),S_vals(end)];
yax = [0,0];
x = S_vals;

%
% Plot the Hugoniot curve
%

hold on
plot(S_vals,tau_vals,'-or','MarkerSize',4,'LineWidth',2);
h = xlabel('S');
set(h,'FontSize',22);
set(h,'FontSize',22);
set(h,'FontSize',22);
try
    plot(Sstar,taustar,'xk','MarkerSize',18,'LineWidth',2);
catch me
end
h = gca;
set(h,'FontSize',22);
set(h,'LineWidth',2);

