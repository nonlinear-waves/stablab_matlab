clc; clear all; close all; beep off; curr_dir = cd;

% load the data
curr_dir = cd;
cd('../data');
ld = load('batch101.mat');
d = ld.d;
cd(curr_dir);

remove = 2;
beta = zeros(1,length(d)-remove);
u1_p = zeros(1,length(d)-remove);
for j = 1:length(d)-remove
    
    p = d{j}.p;
    u1_p(j) = p.u1_p;
    beta(j) = d{j}.beta_vals(end);
        
end

% load the data
curr_dir = cd;
cd('../data');
ld = load('batch102.mat');
d = ld.d;
cd(curr_dir);

remove = 2;
beta2 = zeros(1,length(d)-remove);
u1_p2 = zeros(1,length(d)-remove);
for j = 1:length(d)-remove
    
    p = d{j}.p;
    u1_p2(j) = p.u1_p;
    beta2(j) = d{j}.beta_vals(end);
        
end

beta = [fliplr(beta2),beta];
u1_p = [fliplr(u1_p2),u1_p];

hold on;
plot(1-u1_p,beta,'.-k','LineWidth',2,'MarkerSize',18);
h = xlabel('u_1^--u_1^+');
set(h,'FontSize',18);
h = ylabel('\sigma');
set(h,'FontSize',18);
h = gca;
set(h,'FontSize',18);

