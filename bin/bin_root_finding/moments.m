function y = moments(z,w,pow,mu,order)
% y = moments(z,w,pow,mu,order)
%
% Integrates z^pow*f'(z)/f(z) along a contour given by z where w = f(z).
% The moments are calculated about the complex constant mu. The order of
% the derivative difference approximation may be specified as 2 or 4 and defaults to 4.
%

% NOTE: This is not state of the art code. Using analyticity of the Evans
% function, we can do much better by using Chebyshev interpolation

if nargin < 5
    order=4;
end

if nargin < 4
    mu = 0;
end

[sx,sy]=size(z);
if sx > sy 
   z=z.'; 
end
[sx,sy]=size(w);
if sx > sy
   w=w.';
end

% Remove any repeat points from contours
aind=1;
bind=2;
count=1;
while bind <= length(z)
    if abs(z(aind)-z(bind)) > 10^-14
        index(count)=aind;
        count=count+1;
        aind=bind;
        bind=bind+1;
    else
        bind=bind+1;
    end
end

z=z(index);
w=w(index);

% Make sure the contour has an odd number of points
if mod(length(w),2)==0
   z=[z(1:end-1) 0.5*(z(end-1)+z(end)) z(end)];
   w=[w(1:end-1) 0.5*(w(end-1)+w(end)) w(end)];
end

if order ==2   
    h2=diff(z);
    h1=[0 -1*diff(z(1:end-1))];

    n= length(z);
    intgrand=zeros(1,n);
    for k= 2:n-1           %Here we approximate and store the values of z^pow*f'(z)/f(z)
        intgrand(k)=((z(k)-mu)^pow/w(k))*(1/((h2(k)-h1(k))*h1(k)*h2(k)))*(h2(k)^2*w(k-1)-h1(k)^2*w(k+1)+(h1(k)^2-h2(k)^2)*w(k));
    end

    intgrand(1) = ((z(1)-mu)^pow/w(1))*((w(2)-w(1))/(z(2)-z(1)));
    intgrand(end)=((z(end)-mu)^pow/w(end))*((w(end)-w(end-1))/(z(end)-z(end-1)));    
elseif order == 4
    longer_z = [z(end-2:end-1) z z(2:3)];
    longer_w = [w(end-2:end-1) w w(2:4)];
    n= length(z);
    fprime=zeros(n,1);
    for j=1:n
        mid=j+2;
        h1=longer_z(mid+2)-longer_z(mid);
        h2=longer_z(mid+1)-longer_z(mid);
        h3=longer_z(mid-1)-longer_z(mid);
        h4=longer_z(mid-2)-longer_z(mid);
        vec = [-h3 * h2 * h4 / h1 / (-h4 + h1) / (-h3 + h1) / (-h2 + h1) h4 * h1 * h3 / (-h2 + h1) / h2 / (-h4 + h2) / (-h3 + h2) -h1...
            * h2 * h4 / (h3 ^ 2 - h1 * h3 - h2 * h3 + h1 * h2) / h3 / (-h4 + h3) h2 * h1 * h3 / h4 / (-h4 ^ 3 - h2 * h1 * h4 + h1 * h4 ^ 2 ...
            + h2 * h4 ^ 2 + h3 * h4 ^ 2 - h3 * h1 * h4 - h3 * h2 * h4 + h2 * h1 * h3)];
         fprime(j)= (vec(1)*longer_w(mid+2)+vec(2)*longer_w(mid+1)+vec(3)*longer_w(mid-1)+vec(4)*longer_w(mid-2) ...
                          -longer_w(mid)*(sum(vec)));
    end

    intgrand=zeros(1,n);
    for k= 1:n           %Here we approximate and store the values of z^pow*f'(z)/f(z)
        intgrand(k)=((z(k)-mu)^pow/w(k))*fprime(k);
    end    
end

%Simpson Integration
y = 0;
for j = 1:2:n-2  
x0=z(j);
x1=z(j+1); 
x2=z(j+2);

a = intgrand(j)/( (x0-x1)*(x0-x2) );
b = intgrand(j+1)/( (x1-x0)*(x1-x2) );
c = intgrand(j+2)/( (x2-x0)*(x2-x1) );
 
SUM = (x2^3- x0^3)*(a + b + c)/3 + (x0^2- x2^2)* ( a*(x1+x2) + b*(x0 + x2) + c*(x0+x1) )/2 + (x2-x0)*( a*x1*x2 + b*x0*x2 + c*x0*x1 ); 

y = y + SUM;
end
