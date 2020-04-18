% u = [U V z]'
function cust_ode = ep0integrated_ode(hL,hR,be,D)

  function du = ode(t,u)
    du = zeros(2,1);
    du(1) = 1/be/u(1)^3*(u(1)^3-3*s*u(1)-3*K1-3*u(1)*u(2)*(s*u(1)+3*K1)/(u(1)*u(2)+4*D));
    du(2) = 2*u(2)/u(1)*(s*u(1)+3*K1)/(u(1)*u(2)+4*D);
  end
    s=1/3*(hL^2+hL*hR+hR^2);
    K1=-1/3*hL*hR*(hL+hR);
  
  

    cust_ode = @ode;
end
