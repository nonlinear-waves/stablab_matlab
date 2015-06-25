function out = full_gas_Hugoniot(p,s,var)
% solves the Hugoniot curve 

if strcmp(s.H_ind_var,'tau')
    H = @(p,tau,S)(s.H_fun(p,tau,S));
else
    H = @(p,tau,S)s.H_fun(p,S,tau);
end

left = s.left_H;

Hleft = sign(H(p,var,left));
right = s.right_H;
Hright = sign(H(p,var,right));

if sign(Hleft) == sign(Hright)
    x = linspace(left,right);
    y = zeros(length(x),1);
    for j = 1:length(x)
        y(j) = H(p,var,x(j));
    end
    figure 
    hold on
    plot(x,y,'.-k')
    plot([left,right],[0,0],'-g');
    error('Both end points have the same sign');
end

while (right-left)/left > 1e-14
    
    mid = 0.5*(left+right);
    Hmid = sign(H(p,var,mid));
    if Hleft == Hmid
        left = mid;
    else
        right = mid;
    end
    
end

out = 0.5*(right+left);



