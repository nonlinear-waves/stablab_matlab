clc; clear all; beep off;  close all; drawnow; curr_dir = cd;
file_name = mfilename;

eval(['addpath ',pwd,'/../profile_code']);

%{
Check the concavity along the transition curve.
%}

cd('../data');
ld = load('batch1.mat');
cd('../core_code');
code_dir = cd;

data = ld.data;

eps_vals = 1e-5;
delta = 1e-4;

cnt = 0;
for j = 1:length(data)
    
%     try
%         clc; 
        j
        p = data{j}.p;
        % parameters
        phi = 0.1;
        p.nu = phi;
        p.mu = phi;
        p.eta = -(2/3)*p.mu;
        p.kappa = phi;

        beta_vec = zeros(1,length(eps_vals));

        for k = 1:length(eps_vals)
            beta_vec(k) = verify_concavity2(p,eps_vals(k),delta);
        end
        cnt = cnt + 1
        d{cnt}.beta = beta_vec;
        d{cnt}.p = p;
        d{cnt}.eps_vals = eps_vals;
        d{cnt}.delta = delta;
        
        cd('../data');
        save(file_name,'d');
        cd('../core_code');

%     catch me
%        cd(code_dir); 
%     end
end








cd(curr_dir);