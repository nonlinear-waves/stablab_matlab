function out=moments_roots(domain,range)
% 
% out=moments_roots(domain, range)
%
% Uses the method of moments to find the roots of the analytic function
% F(lambda) evaluated at the points "domain" of the simple, positively
% oriented contour Gamma. Here range=F(domain). The results become less
% accurate generally as the number of roots contained in the contour
% increase. If the contour has winding number of zero, then the function
% throws an error. Generally, it is better to not use this when more than
% one root is inside the contour because of the numerical error.

% NOTE: Todd Kapitula's student studied the method of moments and came up
% with an improvement by integrating differently.

% Created: 17 Dec 2010
% Last updated: 1 Jan 2011

% make sure the dimensions are correct for the moments function.
[sx,sy]=size(domain);
if sx > sy 
   domain=domain.'; 
end
[sx,sy]=size(range);
if sx > sy
   range=range.'; 
end

% find the winding number of the domain contour
wnd=round(winding_number(range));

% if wnd < 4, solve algebrically for the roots using closed formulas
if wnd==1 % Case when there is one root
    out=(moments(domain,range,1,0)/(2*pi*1i)); 
    return
elseif wnd==2 % Case when there are two roots
    mu1=moments(domain,range,1,0)/(4*pi*1i);
    mu2=moments(domain,range,2,mu1)/(2*pi*1i);
    out(1)=mu1+sqrt(mu2/2);
    out(2)=mu1-sqrt(mu2/2);
    return
elseif wnd==3 % Case when there are three roots
    mu1=moments(domain,range,1,0)/(2*pi*1i); %a+b+c
    mu2=moments(domain,range,2,0)/(2*pi*1i); %a^2+b^2+c^2
    mu3=moments(domain,range,3,0)/(2*pi*1i); %a^3+b^3+c^3
    [r1,r2,r3]=CubicRoots(-mu1,.5*(mu1^2-mu2),(-1/6)*(mu1^3-3*mu1*mu2+2*mu3));
    out(1)=r1;
    out(2)=r2;
    out(3)=r3;
    return
end

% calculate the first n moments were n=wnd.
m=zeros(wnd,1);
for j=1:wnd
    m(j,1)=(1/(2*pi*1i))*moments(domain,range,j,0);
end

s=zeros(wnd,1);
s(1)=m(1);
for k=2:wnd
    sums=0;
    for count=1:k-1
        sums=sums+s(k-count)*m(count)*(-1)^mod(count+1,2);
    end
    sums=sums+m(k)*(-1)^mod(k+1,2);
    s(k)=sums/k;
end

% form the polynomial
p=zeros(wnd+1,1);
p(1)=1;
for j=1:wnd
    p(j+1)=(-1)^j*s(j);
end

out=roots(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out1,out2,out3]=CubicRoots(a,b,c)
%[out1,out2,out3]=CubicRoots(a,b,c)
%
% x^3 +ax^2 +bx +c
%Takes the coefficients of a polynomial and returns the roots

%See Abstract Algebra text by Dummit and Foot
p=(1/3)*(3*b-a^2);
q=(1/27)*(2*a^3-9*a*b+27*c);
D=-4*p^3-27*q^2;
%Need to make sure that the cube roots are chosen so that AB=-3p
A=((-27/2)*q+1.5*sqrt(3*D)*1i)^(1/3);
B=((-27/2)*q-1.5*sqrt(3*D)*1i)^(1/3);
%Find best fit of cube roots such that AB=-3p
testAB=zeros(3,3);
for k=1:3
    for j=1:3
       testAB(k,j)=(A*exp(k*2*pi/3*1i))*(B*exp(j*2*pi/3*1i)); 
    end
end
testAB=testAB+3*p;
mini=min(min(abs(testAB)));
[rval,cval]=find(abs(testAB)<mini+eps);
A=A*exp(rval(1)*2*pi/3*1i);
B=B*exp(cval(1)*2*pi/3*1i);
ro=-.5+.5*sqrt(3)*1i;
alpha=(A+B)/3;
beta=(ro^2*A+ro*B)/3;
gamma=(ro*A+ro^2*B)/3;
out1=alpha-a/3;
out2=beta-a/3;
out3=gamma-a/3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

