clc; close all; clear all; beep off;

disp('Plot of the top eigenvalue as a function of k');

% parameters
p.rho_m = 0.1;
p.rho_p = 1; % bigger than p.rho_m
p.h = 1.5;
p.L = 1;

% don't change
p.g = 9.8;
p.a = log(p.rho_p/p.rho_m)/p.h;
profile_ratio = @(x)p.a;
s.options = odeset('RelTol',1e-8,'AbsTol',1e-8);

p.Atwood_number =  p.rho_p -p.rho_m/(p.rho_p + p.rho_m);
fprintf('\nAtwood number = %4.4g\n',p.Atwood_number);

bound = sqrt(p.a*p.g);

% compute evans function


num_k_indices = 200;
num_eig_indices = 1;

M = zeros(num_k_indices,num_eig_indices);
k_vals = zeros(1,num_k_indices);

for j = 1:size(M,1)
    for index = 1:size(M,2)
        k = j/p.L;
        k_vals(j) = k;
        M(j,index) = get_eigenvalues(index,k,p.a,p.h);
    end
end

figure;
hold on;
for j = 1:size(M,2)
    plot(k_vals,M(:,j),'.k','MarkerSize',8);
end
plot([0,num_k_indices],[bound,bound],'--b','LineWidth',2);

mx = max(max(M));
axis([-1,num_k_indices,-0.1*mx,1.1*mx]);
h = xlabel('k');
set(h,'FontSize',18);
h = ylabel('\lambda');
set(h,'FontSize',18);
h = gca;
set(h,'FontSize',18);
