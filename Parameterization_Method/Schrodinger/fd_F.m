function out = fd_F(U_n,U_o,K,H,p)
%
%Function of the finite difference scheme

mu = p.mu;
nu = p.nu;
u_o = U_o(1,2:end-1);
u_om = U_o(1,1:end-2);
u_op = U_o(1,3:end);
u_n = U_n(1,2:end-1);
u_nm = U_n(1,1:end-2);
u_np = U_n(1,3:end);
v_o = U_o(2,2:end-1);
v_om = U_o(2,1:end-2);
v_op = U_o(2,3:end);
v_n = U_n(2,2:end-1);
v_nm = U_n(2,1:end-2);
v_np = U_n(2,3:end);


out = [-mu.*(v_n./2 + v_o./2) + (v_n./2 + v_o./2).*((u_n./2 + u_o./2).^2 + (v_n./2 + v_o./2).^2) + (u_n - u_o)./K + (-2.*v_n + v_nm + v_np)./(2.*H.^2) + (-2.*v_o + v_om + v_op)./(2.*H.^2);
2.*nu.*(v_n./2 + v_o./2) + u_n./2 + u_o./2 - (u_n./2 + u_o./2).*((u_n./2 + u_o./2).^2 + (v_n./2 + v_o./2).^2) + (v_n - v_o)./K - (-2.*u_n + u_nm + u_np)./(2.*H.^2) - (-2.*u_o + u_om + u_op)./(2.*H.^2);
];

