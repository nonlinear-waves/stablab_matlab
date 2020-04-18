function out = fd_jac(U_n,U_o,K,H,p)
%
%Jacobian of finite difference scheme



%Jacobian

u_o = U_o(1,2:end-1);
u_om = U_o(1,1:end-2);
u_op = U_o(1,3:end);
u_n = U_n(1,2:end-1);
u_nm = U_n(1,1:end-2);
u_np = U_n(1,3:end);
%
out{1}{1}{1} = -1./(2.*H.^2);
%
out{1}{1}{2} = -3.*u_n.^2 + 1 + 1./K + H.^(-2);
%
out{1}{1}{3} = -1./(2.*H.^2);
