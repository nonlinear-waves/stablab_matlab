function out = F(x,y,s,p)

   out= [p.c/p.be*(1-y(2,:)-p.be*y(1,:));...
     p.be/p.c*y(2,:)*exp(-1/y(1,:))];

