%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=init_3D(x,y,p,ode,box)

out=ode(x,y,p);

if y(1) > box(1,2)
    out=0;
elseif y(1) < box(1,1)
    out=0;
end
        
if y(2) > box(2,2)
    out=0;
elseif y(2) < box(2,1)
    out=0;
end

if y(3) > box(3,2)
    out=0;
elseif y(3) < box(3,1)
    out=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%