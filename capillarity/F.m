function out = F(x,y,s,p)

out =  [y(2,:);(y(1,:)-s.UL(1)+p.a*((y(1,:)^(-p.gamma)...
    -s.UL(1)^(-p.gamma))) - y(2,:)./y(1,:))/p.d];
