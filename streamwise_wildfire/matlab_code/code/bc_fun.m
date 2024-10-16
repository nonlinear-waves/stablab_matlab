function out = bc_fun(ya,yb,speed,p,s)

AM = profile_jacobian(yb(s.larray),speed,-1,p);
s.LM = orth(projection1(AM,-1,0));

AP = profile_jacobian(yb(s.rarray),speed,1,p);
s.LP = real(orth(projection1(AP,1,0)));

out = [
            ya(s.rarray)-ya(s.larray);              %  matching conditions
            s.LM.' * (yb(s.larray) - s.UL);         % projection at - infinity
            s.LP.'*(yb(s.rarray)-s.UR);             % projection at + infinity
            ya(1)-s.phase;    
         ];





