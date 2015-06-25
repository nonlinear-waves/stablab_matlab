function out = double_F(x,y,s,p,varargin)
% out = double_F(x,y,s,p)
%
% Returns the split domain for the ode given in the function F.
%
% Input "x" and "y" are provided by the ode solver.Note that s.rarray
% should be [1,2,...,k] and s.larray should be [k+1,k+2,...,2k]. See
% STABLAB documentation for more inforamtion about the structure s.

if size(varargin,2) > 0
    out = [(s.R/s.I)*s.F((s.R/s.I)*x,y(s.rarray,:),s,p,varargin);(s.L/s.I)*s.F((s.L/s.I)*x,y(s.larray,:),s,p,varargin)];
else
    out = [(s.R/s.I)*s.F((s.R/s.I)*x,y(s.rarray,:),s,p);(s.L/s.I)*s.F((s.L/s.I)*x,y(s.larray,:),s,p)];
end

