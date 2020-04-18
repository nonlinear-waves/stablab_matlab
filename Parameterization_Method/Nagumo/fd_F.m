function out = fd_F(U_n,U_o,K,H,p)
%
%Function of the finite difference scheme

u_o = U_o(1,2:end-1);
u_om = U_o(1,1:end-2);
u_op = U_o(1,3:end);
u_n = U_n(1,2:end-1);
u_nm = U_n(1,1:end-2);
u_np = U_n(1,3:end);


out = [-u_n.^3 + u_n + (u_n - u_o)./K - (-2.*u_n + u_nm + u_np)./(2.*H.^2) - (-2.*u_o + u_om + u_op)./(2.*H.^2);
];

