clc; clear all; close all; beep off; home_dir = cd;

% load the critical curve data
cd('../data');
data_dir = cd;
ld = load('batch41.mat');
stored_data = ld.d;
cd('../core_code');
code_dir = cd;

% file name that current data is saved to
file_name = mfilename;

eps_vals = fliplr(linspace(1e-6,1e-3,10));

end_point_tol = 1e-6;
for q = 1:length(stored_data)
    
    clc;
    fprintf('\n\nPercent done = %.3g\n\n',100*((q-1)/length(stored_data)));
    
    p = stored_data{q}.p;    
    p.hstar = p.h1;

    % solve for the profile
    end_point_tol = 1e-8;
    [s,p] = get_profile(p,end_point_tol);
    
    beta_vals = zeros(1,length(eps_vals));
    success_cnt = 0;
    try
        for j = 1:length(eps_vals)
    
        
            % finite difference controls
            eps = eps_vals(j);
            h = eps;

            % Evans function controls
            s.system = 'parallel';
            s.mat_type = 'mbfv';
            s.sys_fun = 'A_aux_mhd_alpha';
            evan_mat = @A_pseudo_mhd_alpha;
            s.xi = eps;

            [s,e,m,c] = emcset(s,'front',LdimRdim(evan_mat,s,p),'reg_reg_polar',evan_mat);
            m.method = @drury;
            c.refine = 'off';
            c.ksteps = 2^8;
            c.lambda_steps = 0;
            c.Lproj = @projection1;
            c.Rproj = @projection1;
            m.options = odeset('RelTol',1e-13,'AbsTol',1e-13,'Refine',1,'Stats','off');


            % get the basis on which to project - D(xi,lambda)
            s.xi = eps;
            lambda = 0; 
            [proj_L, Q_L] = c.Lproj(evan_mat(s.L,lambda,s,p),1,0);
            [proj_R, Q_R] = c.Rproj(evan_mat(s.R,lambda,s,p),-1,0);


            % compute D(h+eps,0)
            s.xi = h+eps;
            lambda = 0;
            [proj_L, Q1] = c.Lproj(evan_mat(s.L,lambda,s,p),1,0);
            [proj_R, Q2] = c.Rproj(evan_mat(s.R,lambda,s,p),-1,0);
            L_basis = proj_L*Q_L;
            R_basis = proj_R*Q_R;
            D1 = evans(L_basis,R_basis,lambda,s,p,m,e);

            % compute D(eps,0)
            s.xi = eps;
            lambda = 0;
            [proj_L, Q1] = c.Lproj(evan_mat(s.L,lambda,s,p),1,0);
            [proj_R, Q2] = c.Rproj(evan_mat(s.R,lambda,s,p),-1,0);
            L_basis = proj_L*Q_L;
            R_basis = proj_R*Q_R;
            D2 = evans(L_basis,R_basis,lambda,s,p,m,e);

            % compute D(eps,eps*h)
            s.xi = eps;
            lambda = eps*h;
            [proj_L, Q1] = c.Lproj(evan_mat(s.L,lambda,s,p),1,0);
            [proj_R, Q2] = c.Rproj(evan_mat(s.R,lambda,s,p),-1,0);
            L_basis = proj_L*Q_L;
            R_basis = proj_R*Q_R;
            D3 = evans(L_basis,R_basis,lambda,s,p,m,e);

            beta = -(D1-D2)/(D3-D2);
            beta_vals(j) = beta;
            success_cnt = success_cnt+1;
        end
    catch me
        cd(code_dir);
    end
        
        d{q}.p = p;
        d{q}.eps_vals = eps_vals(1:success_cnt);
        d{q}.beta_vals = beta_vals(1:success_cnt);
            
    

    
%     hold on;
%     plot((eps_vals),(real(beta_vals)),'.-b','LineWidth',2);
%     plot((eps_vals),(imag(beta_vals)),'.-r','LineWidth',2);
%     
%     h = xlabel('\epsilon');
%     set(h,'FontSize',18);
%     h = ylabel('\beta');
%     set(h,'FontSize',18);
    
    
    cd('../data');
    save(file_name,'d');
    cd('../core_code');

end





