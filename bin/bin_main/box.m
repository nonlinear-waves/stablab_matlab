function preimage = box(left,right,bottom,top,left_pnts,right_pnts,bottom_pnts,top_pnts,ksteps,lambda_steps)
% out = semisquare(box,left_pnts,top_pnts,right_pnts,ksteps,lambda_steps)
%
% Returns a half-contour on which to compute the Evans function. 
%
% Input: a = left coordinate of half box, b = right coordinate of half box, c = top
% coordiante of half box


% Determine if lambda_steps are specified and default to zero if not.
if nargin < 10
    lambda_steps = 0;
end


% Specify the number of points
pr = (right_pnts-1)*ksteps+right_pnts+((right_pnts-1)*ksteps+right_pnts-1)*lambda_steps;
pt = (top_pnts-1)*ksteps+top_pnts+ksteps+1+((top_pnts-1)*ksteps+top_pnts+ksteps)*lambda_steps;
pl =(left_pnts-1)*ksteps+left_pnts+ksteps+1+((left_pnts-1)*ksteps+left_pnts+ksteps)*lambda_steps;
pb = (bottom_pnts-1)*ksteps+bottom_pnts+ksteps+1+((bottom_pnts-1)*ksteps+bottom_pnts+ksteps)*lambda_steps;

% Construct and combine the parts of the contour.
rr = right+1i*linspace(bottom,top,pr);
tt = top*1i+linspace(right,left,pt);
ll = left+1i*linspace(top,bottom,pl);
bb = bottom*1i+linspace(left,right,pb);
preimage=[ rr,tt(2:end) , ll(2:end), bb(2:end) ];





