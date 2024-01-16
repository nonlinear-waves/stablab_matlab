
%{
    This is the main file. 
%}

clc; close all; clear all; beep off; curr_dir = cd; t_start = tic;

% -------------------------------------------------------------------------
% Controls
% -------------------------------------------------------------------------



for j = 280:5:400


   % left bound on the interval in the parameter nu
   nu_left = ['0.',num2str(j)];
   disp(nu_left);

   % right bound on the interval in the parameter nu
   nu_right = ['0.',num2str(j+1)];

   prover(nu_left,nu_right);

end






















