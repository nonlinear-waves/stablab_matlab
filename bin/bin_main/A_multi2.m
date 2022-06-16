function out = A_multi2(x,lambda,s,p)

xi = s.xi;

fun = str2func(s.sys_fun);
[nhyp,nsys,A0,A1,Atilde,Einv,b_21,bj,bjk,A1_12_der,a,a_der] = fun(x,s,p);

q.nhyp = nhyp;
q.nsys = nsys;

% choose matrix type
if strcmp(s.mat_type,'mbfv') % modified balanced flux variables
    rfun = norm(xi)+lambda;
else
    rfun = 1;
end

% make transformation to sharp variables
xi = xi/rfun;
lambda = lambda/rfun;

% create Axi
Axi = zeros(size(Atilde{1}));
for j = 1:length(Atilde)
    Axi = Axi + xi(j)*Atilde{j};
end

% create bxi
bxi = zeros(nsys-nhyp,nsys);
for j = 1:length(Atilde)
    bxi = bxi + xi(j)*bj{j};
end
bxi_21 = bxi(:,1:nhyp);
bxi_22 = bxi(:,nhyp+1:nsys);

% create bxixi
bxixi = zeros(nsys-nhyp,nsys);
for j = 1:size(bjk,1)
    for k = 1:size(bjk,2)
        bxixi = bxixi + xi(j)*xi(k)*bjk{j}{k};
    end
end
bxixi_21 = bxixi(:,1:nhyp);
bxixi_22 = bxixi(:,nhyp+1:nsys);

Csharp = -lambda*m11(A0,q)*a-1i*m11(Axi,q)*a;

Dsharp = -lambda*m11(A0,q)*a*m12(A1,q)+lambda*m12(A0,q) -...
    1i*m11(Axi,q)*a*m12(A1,q)+1i*m12(Axi,q);

Hsharp = -lambda*m21(A0,q)*a-1i*m21(Axi,q)*a-rfun*bxixi_21*a;

Jsharp = -lambda*m21(A0,q)*a*m12(A1,q)+lambda*m22(A0,q)-...
    1i*m21(Axi,q)*a*m12(A1,q)+1i*m22(Axi,q)-rfun*bxixi_21*a*m12(A1,q)+...
    rfun*bxixi_22;

Fsharp = b_21*a_der+rfun*b_21*a*Csharp+1i*rfun*bxi_21*a-m21(A1,q)*a;

Gsharp = b_21*a_der*m12(A1,q)+b_21*a*A1_12_der-1i*rfun*bxi_22+m22(A1,q)+...
                rfun*b_21*a*Dsharp+1i*rfun*bxi_21*a*m12(A1,q)-m21(A1,q)*a*m12(A1,q);

 out = [
   rfun*Csharp,   zeros(nhyp,nsys-nhyp), Dsharp;
   rfun*Hsharp, zeros(nsys-nhyp,nsys-nhyp), Jsharp;
   rfun*Einv*Fsharp, rfun*Einv, Einv*Gsharp
 ];
 
%%%%%%%%%%%%%%%%%%%%%%%%

% (1:nhyp,1:nhyp)           A_11
% (1:nhyp,nhyp+1:nsys)       A_12
% (nhyp+1:nsys,1:nhyp)       A_21
% (nhyp+1:nsys,nhyp+1:nsys)   A_22

function out = m11(A,q)
    out = A(1:q.nhyp,1:q.nhyp);

function out = m12(A,q)
    out = A(1:q.nhyp,q.nhyp+1:q.nsys);

function out = m21(A,q)
    out = A(q.nhyp+1:q.nsys,1:q.nhyp);
    
function out = m22(A,q)
    out = A(q.nhyp+1:q.nsys,q.nhyp+1:q.nsys);
    
    
    
    
    
