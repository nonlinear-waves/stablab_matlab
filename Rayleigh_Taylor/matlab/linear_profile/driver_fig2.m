clc; close all; clear all; beep off;

disp('Plot of the second largest eigenvalue as a function of k');

ld = load('rts1');
rts1 = ld.rts1;

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

choices = 2; % 1 for top eigenvalue, 2 for second top eigenvalue
k_max = 685;

left = 0.001;
right = 2;
k_vals = (1:1:k_max)/p.L;
rts2 = zeros(size(k_vals));

figure;
hold on;
plot(k_vals,rts1(1:end-1),'.r','MarkerSize',8);
hx = xlabel('k');
set(hx,'FontSize',18);
hy = ylabel('\lambda_1(k)');
set(hy,'FontSize',18);
h = gca;
set(h,'FontSize',18);
plot([k_vals(1),k_vals(end)],[sqrt(p.g/p.L0),sqrt(p.g/p.L0)],'--b','LineWidth',2);

for j = 1:length(k_vals)-1

    disp(j);
    p.k = k_vals(j);
    rts2(j) = get_eigenvalues(p,choices,left,right);
    left = rts2(j);

    if j > 1
        right = min(2*rts2(j)-rts2(j-1),rts1(j+1)-1e-6);
    else
        choices = 1;
        right = rts1(2)-1e-6;
    end
    

    plot(p.k,rts2(j),'.k','MarkerSize',18);


end


save('rts2','rts2');

figure;
hold on;
plot(k_vals,rts2,'.k','MarkerSize',8);
hx = xlabel('k');
set(hx,'FontSize',18);
hy = ylabel('\lambda_2(k)');
set(hy,'FontSize',18);
h = gca;
set(h,'FontSize',18);
plot([k_vals(1),k_vals(end)],[sqrt(p.g/p.L0),sqrt(p.g/p.L0)],'--b','LineWidth',2);

axis([-2,k_vals(end),0,1.1*max(rts2)]);






