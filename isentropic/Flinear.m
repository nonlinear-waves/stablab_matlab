function out = Flinear(v,p)

out=v*(v-1+p.a*(v^(-p.gamma)-1));

out = (v-1+p.a*(v^(-p.gamma)-1)) + v*(1-p.a*p.gamma*v^(-p.gamma-1));

