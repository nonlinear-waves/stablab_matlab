function [w,prew] = lop_for_transition(p,preimage)

%
% dependent parameter
%

p.u1_m = 1;
p.rho_m = 1;
p.rho_p = 1/p.u1_p;
p.a0 = p.u1_p^p.gamma*(1-p.u1_p)/(1-p.u1_p^p.gamma);

%
% initialize stablab structures
%

s.xi = 1;
s.R = 1;
s.L = -1;
evan_mat = @A_lop_alpha; % USE THIS ONE
[s,e,m,c] = emcset(s,'lopatinski',LdimRdim(evan_mat,s,p),'default',evan_mat);
e.jump = @lop_jump;

% refine the Evans function computation to achieve set relative error
c.refine = 'off';
% display a waitbar
c.stats = 'off';
c.debug = 'off';
c.ksteps = 0;
c.basisL = @analytic_basis_local;
c.basisR = @analytic_basis_local;
c.Lproj= @projection_local;
c.Rproj= @projection_local;

%
% preimage contour
%

[w,prew] = contour(c,s,p,m,e,preimage);
w = real(w/w(1));
% w = real(w);
% figure
% plot(prew,w,'.-k');
% drawnow;



