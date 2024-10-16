function out = guess(x,s)

q = s.scale;

RHS1 = (s.UL+s.UR)/2 + (s.UR-s.UL)*tanh(q*x)/2;
RHS2 =  q*(s.UR-s.UL)*sech(q*x)^2/2;


LHS1 = (s.UL+s.UR)/2 + (s.UR-s.UL)*tanh(-q*x)/2;
LHS2 = q*(s.UR-s.UL)*sech(-q*x)^2/2;

out = [RHS1(1);RHS2(1);LHS1(1);LHS2(1)];