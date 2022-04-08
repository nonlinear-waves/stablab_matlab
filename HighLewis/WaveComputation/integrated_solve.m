% depends on integrated_ode.m, integrated_ode.m
function sol = integrated_solve(c,be,initial)
  options = odeset('RelTol',1e-13,'AbsTol',1e-13);
  final = 100;
  x_values = [0,final];

  
    ode = integrated_ode(c,be);
  
  sol = ode45(ode, x_values, initial, options);
  test=0;
  while test==0
    if final > 560
      disp('giving up');
      break;
      % if c is too high decline is very slow
      % assume we'd converge by now if that's not the case
    end
    if isnan(sol.y(2,end)) | isinf(sol.y(2,end))
        test=1;
        indf=find(isnan(sol.y(2,:)) | isinf(sol.y(2,:)), 1,'first')-1;
        final=sol.x(indf);
        sol = ode45(ode, [0,final], initial, options);
    elseif sol.y(1,end)<10^(-2)
        test=1;
    else
        final = final * 2;
        disp(['extending to ', num2str(final)]);
        sol = odextend(sol, [], final);
        if isnan(sol.y(2,end)) | isinf(sol.y(2,end))
            test=1;
            indf=find(isnan(sol.y(2,:)) | isinf(sol.y(2,:)), 1,'first')-1;
            final=sol.x(indf);
            sol = ode45(ode, [0,final], initial, options);
        end;
    end;
  end;
end

