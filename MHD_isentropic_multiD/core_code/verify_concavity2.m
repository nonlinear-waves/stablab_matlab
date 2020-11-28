function beta = verify_concavity2(p,eps,delta)

% profile
s.something = 0;
[s,p] = profile_solve_pseudo(s,p);

if nargin < 2
    eps = 1e-5;
end

% Evans function controls
s.system = 'parallel';
s.mat_type = 'mbfv';
s.sys_fun = 'A_aux_mhd_alpha';
evan_mat = @A_pseudo_mhd_alpha;
s.xi = 0;
[s,e,m,c] = emcset(s,'front',LdimRdim(evan_mat,s,p),'reg_adj_polar',evan_mat);
m.method = @drury;
c.refine = 'off';
c.ksteps = 0;
c.lambda_steps = 0;
c.Lproj = @projection1;
c.Rproj = @projection1;
m.options = odeset('RelTol',1e-13,'AbsTol',1e-13,'Refine',1,'Stats','off');

% basis at (lambda = 0, xi = eps)
s.xi = eps;
lambda = delta;
[proj_L, fixed_basis_L] = c.Lproj(c.LA(s.L,lambda,s,p),1,0);
[proj_R, fixed_basis_R] = c.Rproj(c.RA(s.R,lambda,s,p),-1,0);
yl = proj_L*fixed_basis_L;
yr = proj_R*fixed_basis_R;
D_0_eps = evans(yl,yr,lambda,s,p,m,e);


% basis at (lambda = eps, xi = eps)
s.xi = eps;
lambda = delta+eps;
[proj_L, unused] = c.Lproj(c.LA(s.L,lambda,s,p),1,0);
[proj_R, unused] = c.Rproj(c.RA(s.R,lambda,s,p),-1,0);
yl = proj_L*fixed_basis_L;
yr = proj_R*fixed_basis_R;
D_eps_eps = evans(yl,yr,lambda,s,p,m,e);

% basis at (lambda = 0, xi = 2*eps)
s.xi = 2*eps;
lambda = delta;
[proj_L, unused] = c.Lproj(c.LA(s.L,lambda,s,p),1,0);
[proj_R, unused] = c.Rproj(c.RA(s.R,lambda,s,p),-1,0);
yl = proj_L*fixed_basis_L;
yr = proj_R*fixed_basis_R;
D_0_2eps = evans(yl,yr,lambda,s,p,m,e);


% basis at (lambda = eps, xi = 0)
s.xi = 0;
lambda = delta + eps;
[proj_L, unused] = c.Lproj(c.LA(s.L,lambda,s,p),1,0);
[proj_R, unused] = c.Rproj(c.RA(s.R,lambda,s,p),-1,0);
yl = proj_L*fixed_basis_L;
yr = proj_R*fixed_basis_R;
D_eps_0 = evans(yl,yr,lambda,s,p,m,e);

% basis at (lambda = 0, xi = 0)
s.xi = 0;
lambda = delta;
[proj_L, unused] = c.Lproj(c.LA(s.L,lambda,s,p),1,0);
[proj_R, unused] = c.Rproj(c.RA(s.R,lambda,s,p),-1,0);
yl = proj_L*fixed_basis_L;
yr = proj_R*fixed_basis_R;
D_0_0 = evans(yl,yr,lambda,s,p,m,e);

beta = (D_0_2eps-2*D_0_eps+D_0_0)/(D_eps_eps-D_eps_0-D_0_eps+D_0_0);

