function [out,Dout] = P(rho,p)
% Pressure function for isentropic system in Eulerian coordinates

out = p.a*rho^p.gamma;

Dout = p.a*p.gamma*rho^(p.gamma-1);