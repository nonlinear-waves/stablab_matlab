
%{
    This is the main file. 
%}

clc; close all; clear all; beep off; curr_dir = cd; t_start = tic;

% -------------------------------------------------------------------------
% Controls
% -------------------------------------------------------------------------

nu_vals = linspace(0.28,3.2,30);

scale_R = 0.5;

min_x_R = 4; 

figure;
hold on;
h = xlabel('\nu');
set(h,'FontSize',18);
h = ylabel('\theta_0');
set(h,'FontSize',18);
h = gca;
set(h,'FontSize',18);



for j = 1:length(nu_vals)

    [x_R,theta_0,theta_R,r_R] = get_x_R(nu_vals(j),scale_R,min_x_R);

    plot(nu_vals(j),theta_0,'.k','MarkerSize',8);
    drawnow;

end











