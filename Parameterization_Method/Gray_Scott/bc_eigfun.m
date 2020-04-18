function out = bc_eigfun(ya,yb,lambda,p)





mu1 = sqrt(lambda+p.alpha);
mu2 = sqrt(lambda+1/p.gamma);

z1 = [-mu1;1;0;0];
z2 = [0;0;-mu2;1];

z3 = [mu1;1;0;0];
z4 = [0;0;mu2;1];


out = [[z1,z2].'*ya(1:4);[z3,z4].'*yb(5:8);norm(yb(1:4))^2-1;yb(1:4)-ya(5:8)];



