function Bound = BoundComputation(be)
%Computation of the bound for the absolute value of positive eigenvalues

[c, front, sol] = integrated_find_c(be);

n=10000;
dxi=(sol.x(end)-sol.x(1))/n;

xi=linspace(sol.x(1),sol.x(end),n);

v=front(xi);
u=v(:,1);
y=v(:,2);

%Computing mu for each i and storing those in the vector mm
mm=zeros(length(xi),1);
for i=1:length(xi)
mm(i)=mu(i);
end;

%Checking that a denominator does not become too small for Matlab to handle
if min(mm(:))<10^(-150)
    indf=find(mm<10^(-150),1,'first');
else
    indf=length(xi);
end;
            

u=v(1:indf,1);
y=v(1:indf,2);
mm(indf+1:end)=[];

%Computing mu2 (see below) for each i and storing those in the vector mm2
mm2=zeros(indf,1);
for i=1:indf
mm2(i)=mu2(i);
end;

%Bound=c^2/4;

Bound=c^2/4+max(y./u.^2.*exp(-1./u))+sqrt(2)*be/c*sqrt(sum(exp(-2./u)./mm.^2.*mm2.^2)*dxi);

%Function that defines tilde mu used in the bound
function m = mu(i)
    us=u(i:end);
    ys=y(i:end);
    m=exp(be/c*sum(exp(-1./us))*dxi);
end

%Function that defines an expression part of the bound 
function m2 = mu2(i)
    us=u(i:end);
    ys=y(i:end);
    mms=mm(i:end);
    m2=sum(mms.^2.*ys.^2./us.^4.*exp(-2./us))*dxi;
end






end

