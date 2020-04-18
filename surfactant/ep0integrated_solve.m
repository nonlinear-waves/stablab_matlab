% depends on integrated_ode_eps0.m, integrated_ode.m
function sol = ep0integrated_solve(hL,hR,D,be,L,initial)

    options = odeset('RelTol',1e-13,'AbsTol',1e-13);

    x_values = [0 L];
  
    ode = ep0integrated_ode(hL,hR,be,D);
    sol = ode45(ode, x_values, initial, options);
 
end
