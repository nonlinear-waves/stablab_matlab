function [fv,df] = f(u,p)

[r,dr] = react(u);

fv = p.vstar*r-p.l*u;
df = p.vstar*dr-p.l;
