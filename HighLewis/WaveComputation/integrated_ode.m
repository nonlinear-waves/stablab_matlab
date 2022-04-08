% u = [U V]'
function cust_ode = integrated_ode(c,be)
  function du = ode(t,u)
    du = zeros(2,1);
    du(1) = c/be * (1-u(2) - be*u(1));
    du(2) = (be/c)*u(2)*exp(-1/u(1));
  end


    cust_ode = @ode;
end
