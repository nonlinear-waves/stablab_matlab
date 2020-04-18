function out = continuation_guess(x,s_old,s_new)
% out = continuation_guess(x,s_old,s_new)
%
% Ouput gives initial guess for boundary value solver at x where s_old is
% the standard stablab structure s for the previously solved boundary value
% solution and s_new is that for the solution being solved. If v is the
% solution corresponding to s_new and y is the solution corresponding to
% s_old, continuation_guess yields as output v=a*y+b done componenet wise
% to allow phase conditions to be specified and so v matches its end
% states.
%
% The end state and phase condition can not match.

y=deval(x,s_old.sol);
out=zeros(length(y),1);

% coefficients for the guess for the new function v \approx a*y+b done
% componentswise. Positive infinity corresponds to the first column of the
% coefficeint matrices, and negative infinity corresponds to the second
% column of coefficient matrices.
a=zeros(length(s_old.rarray),2);
b=zeros(length(s_old.rarray),2);

% find scaling coefficients
for j=1:length(s_old.rarray)
    
    % determine if the phase condition should be specified for the jth
    % component
    specify_phase=0;
    for k=1:length(s_new.order)
        if j==s_new.order(k)
            specify_phase=1;
            phase_index=j;
        end
    end
    
    % determine coefficients based on type
    vminus=s_new.UL(j);
    vplus=s_new.UR(j);
    yminus=s_old.UL(j);
    yplus=s_old.UR(j);
    
    % Note that it is important that the end state and phase condition not
    % be the same.
    if specify_phase==1
        % case where the phase condition is specified
        vnot=s_new.phase(phase_index);
        ynot=s_old.phase(phase_index);
        vec = [yplus 1 0 0; ynot 1 0 0; 0 0 yminus 1; 0 0 ynot 1]\[vplus; vnot; vminus; vnot];
        a(j,1)=vec(1);
        b(j,1)=vec(2);
        a(j,2)=vec(3);
        b(j,2)=vec(4);
    else
        % case where the phase condition is not specified
        if yplus==yminus
            a(j,1)=1;
            b(j,1)=0;
            a(j,2)=1;
            b(j,2)=0;
        else
            vec=[yplus 1; yminus 1]\[vplus; vminus];
            a(j,1)=vec(1);
            b(j,1)=vec(2);
            a(j,2)=vec(1);
            b(j,2)=vec(2);
        end
    end
end

% make the affine transformation, v=a*y+b
out(s_old.rarray) = a(:,1).*y(s_old.rarray)+b(:,1);
out(s_old.larray)=a(:,2).*y(s_old.larray)+b(:,2);


