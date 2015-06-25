function out = ode_to_bvp_guess(x,sol,s)

out = [ deval(sol,(s.R/s.I)*x); deval(sol,(s.L/s.I)*x)];



