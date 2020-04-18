function [frontn, frontp] = ep0integrated(hL,hR,D,be,L,maxG)
 
    s=1/3*(hL^2+hL*hR+hR^2);
    K1=-1/3*hL*hR*(hL+hR);
  
    initial=[-3*K1/s,maxG];
    solp = ep0integrated_solve(hL,hR,D,be,L,initial);

    initial=[-3*K1/s,maxG];
    soln = ep0integrated_solve(hL,hR,D,be,-L,initial);
    
  function yp = eval_frontp(xi)
        yp = deval(solp, xi);     
  end

  function yn = eval_frontn(xi)
      yn = deval(soln, xi);
  end

  frontn = @eval_frontn;

  frontp = @eval_frontp;  
  
end
