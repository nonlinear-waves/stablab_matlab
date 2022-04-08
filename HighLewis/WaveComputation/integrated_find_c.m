function [c, front, sol] = integrated_find_c(be)
  c_low = 0;
  c_high = 1;
  initial_scale = 1e-3;

  function v = evector(c)
      v = [-c^2 be*(c^2 +be*exp(-be))];
  end

  function v = evalue(c)
      v = be*exp(-be)/c;
  end

  precision = 1e-6;

  
    fixpt = [1/be,0];
    
    test_c = (c_low + c_high)/2;
    
    initial = fixpt + initial_scale * evector(test_c);

    
    


  while c_high - c_low > precision
    test_c = (c_low + c_high)/2;
    initial = fixpt + initial_scale * evector(test_c);
    sol = integrated_solve(test_c,be,initial);
    

    
    terminus=sol.y(2,end);
    
    
    
    
    
    if terminus < 1
        c_high = test_c;
    else
        c_low = test_c;
    end
    disp(sprintf('%f ... %f', c_low, c_high));
  end

   

  function y = eval_front(xi)
  if xi >= 0
      y = deval(sol, xi)'; % if xi > endpoint this fails. handle somehow?
    else
      y = fixpt - exp(xi*evalue(c)) * initial_scale * evector(c);
    end
  end

  c=test_c; 


      

      
      

figure(1);

  plot(sol.y(1,:),sol.y(2,:));
  front = @eval_front;
  
  hold on;
  
  test_c=.65;
  sol2 = integrated_solve(test_c,be,initial);
  plot(sol2.y(1,:),sol2.y(2,:),':');
  %box off;
  
  hold on;
  
  test_c=.5;
  sol3 = integrated_solve(test_c,be,initial);
  plot(sol3.y(1,:),sol3.y(2,:),'--');
  xlabel('u','FontSize',24)
  ylabel('y','FontSize',24)
  axis([-0.1,1.1,-.1,1.2])
  
  figure(2);
  
  plot(sol.x(1:end),sol.y(1,1:end));
  hold on;
  plot(sol.x(1:end),sol.y(2,1:end),'--');
  %axis([0,sol.x(end),-.1,1.1])
  xlabel('\xi','FontSize',24);
  %box off;
  
  
  
  
    
  
  
  
end
