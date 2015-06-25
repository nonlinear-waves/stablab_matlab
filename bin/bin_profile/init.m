function out=init(x,y,box,ode)
% global box ode
pass=1;
if y(1)<box(1)
    pass=0;
end
if y(1)>box(2)
    pass=0;
end
if y(2)<box(3)
    pass=0;
end
if y(2)>box(4)
    pass=0;
end

if pass==0
    out=[0;0];
else
    out=ode(x,y);
end
   
