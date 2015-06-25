function [preimage]=semicirc2(circpnts,imagpnts,innerpnts,ksteps,r,spread,inner_radius,lambda_steps)
% [preimage]=semicirc2(circpnts,imagpnts,innerpnts,ksteps,R,spread,inner_radius,lambda_steps)
%
% Returns a quarter of a semicircle with a small inner circle skirting the
% origin.
%
% Input "circpnts" is the number of points on the circle part, "imagpnts"
% is the number of points on the imaginary axis, "ksteps" is the number of
% kato steps taken between points, "R" is the
% radius of the semicircle, "spread" is a constant that spreads the points
% on the imaginary axis so that they are more dense near the origin, and
% "inner_radius" is how close along the imaginary axis the contour comes to the
% origin. If not specified, "lambda_points" defaults to zero. If specified,
% "lambda_points" is the number of contour points between Kato steps. That
% is, the analytic basis is not computed on the additional lambda_steps
% points unless achieving relative error tolerances of the Evans function ouput
% requires it. If that is the case, the analytic Kato basis is only
% computed on the region where additional Evans function evaluations are
% needed.

% If lambda_points are not specified, default to zero.
if nargin == 7
    lambda_steps = 0;
end

% Specify points
p1=(circpnts-1)*ksteps+circpnts + (((circpnts-1)*ksteps+circpnts)-1)*lambda_steps ;
p2=(imagpnts-1)*ksteps+imagpnts+ksteps+1+(((imagpnts-1)*ksteps+imagpnts+ksteps+1)-1)*lambda_steps;
p3=innerpnts*ksteps+innerpnts+1 + ((innerpnts*ksteps+innerpnts+1)-1)*lambda_steps;

% construct and combine contour parts
theta = linspace(0,pi/2,p1);
ln=linspace(r^(1/spread),inner_radius^(1/spread),p2).^spread;
theta2=linspace(pi/2,0,p3);
preimage=[r*exp(1i*theta) ln(2:end)*1i inner_radius*exp(1i*theta2(2:end))];
    
    
    
    
    
    
    
    