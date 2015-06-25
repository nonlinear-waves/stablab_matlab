function preimage = semibox(left_loc,right_loc,top_loc,left_pnts,top_pnts,right_pnts,ksteps,lambda_steps,zero_dist)
% out = semisquare(box,left_pnts,top_pnts,right_pnts,ksteps,lambda_steps)
%
% Returns a half-contour on which to compute the Evans function. 
%
% Input: a = left coordinate of half box, b = right coordinate of half box, c = top
% coordiante of half box


% Determine if lambda_steps are specified and default to zero if not.
if nargin == 6
    lambda_steps = 0;
end


% Specify the number of points
pr = (right_pnts-1)*ksteps+right_pnts+((right_pnts-1)*ksteps+right_pnts-1)*lambda_steps;
pt = (top_pnts-1)*ksteps+top_pnts+ksteps+1+((top_pnts-1)*ksteps+top_pnts+ksteps)*lambda_steps;
pl =(left_pnts-1)*ksteps+left_pnts+ksteps+1+((left_pnts-1)*ksteps+left_pnts+ksteps)*lambda_steps;

% Construct and combine the parts of the contour.
rr = right_loc+1i*linspace(0,top_loc,pr);
tt = top_loc*1i+linspace(right_loc,left_loc,pt);
ll = left_loc+1i*linspace(top_loc,zero_dist,pl);
preimage=[ rr,tt(2:end) , ll(2:end) ];





