clc; close all; clear all; beep off;

disp('Plot of the largest eigenvalue as a function of k');

% parameters
p.rho_m = 0.1;
p.rho_p = 0.2; % bigger than p.rho_m
p.h = 2;
p.L = 1;

% fixed
p.g = 9.8;
p.rho_0_der = (p.rho_p-p.rho_m)/p.h;
rho_0 = p.rho_m +p.rho_0_der*(linspace(-p.h,0,1000)+p.h);
p.L0 = 1/max(p.rho_0_der./rho_0);

choices = 1; % 1 for top eigenvalue, 2 for second top eigenvalue
k_max = 686;


left = 0.01;
right = 2;
k_vals = (1:1:k_max)/p.L;
rts1 = zeros(size(k_vals));
for j = 1:length(k_vals)
    disp(j);
    p.k = k_vals(j);
    rts1(j) = get_eigenvalues(p,choices,left,right);
    left = rts1(j);
    if j>1
        right = 2*rts1(j)-rts1(j-1);
    else
        right = 2;
    end

end

save('rts1','rts1');

figure;
hold on;
plot(k_vals,rts1,'.k','MarkerSize',8);
hx = xlabel('k');
set(hx,'FontSize',18);
hy = ylabel('\lambda_1(k)');
set(hy,'FontSize',18);
h = gca;
set(h,'FontSize',18);
plot([k_vals(1),k_vals(end)],[sqrt(p.g/p.L0),sqrt(p.g/p.L0)],'--b','LineWidth',2);

axis([-2,k_vals(end),0,1.1*max(rts1)]);






