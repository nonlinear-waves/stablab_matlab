function [p,s] = get_profile(p,s)

speed_guess = -0.0854708;

s.I = 300;
s.scale = 0.1;
guess_fun = @(x)guess(x,s);
s = solve_profile(s,p,speed_guess,guess_fun);
p.c = s.sol.parameters;