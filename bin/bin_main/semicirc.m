function [preimage]=semicirc(circpnts,imagpnts,ksteps,R,spread,zerodist,lambda_steps)
% [preimage]=semicirc(circpnts,imagpnts,ksteps,R,spread,zerodist,lambda_steps)
%
% Returns a quarter of a semicircle
%
% Input "circpnts" is the number of points on the circle part, "imagpnts"
% is the number of points on the imaginary axis, "ksteps" is the number of
% kato steps taken between contour points, "R" is the
% radius of the semicircle, "spread" is a constant that spreads the points
% on the imaginary axis so that they are more dense near the origin, and
% "zerodist" is how close along the imaginary axis the contour comes to the
% origin. If not specified, "lambda_points" defaults to zero. If specified,
% "lambda_points" is the number of contour points between Kato steps. That
% is, the analytic basis is not computed on the additional lambda_steps
% points unless achieving relative error tolerances of the Evans function ouput
% requires it. If that is the case, the analytic Kato basis is only
% computed on the region where additional Evans function evaluations are
% needed.

% Determine if lambda_steps are specified and default to zero if not.
if nargin == 6
    lambda_steps = 0;
end

% Specify the number of points
p1=(circpnts-1)*ksteps+circpnts+((circpnts-1)*ksteps+circpnts-1)*lambda_steps;
p2=(imagpnts-1)*ksteps+imagpnts+ksteps+1+((imagpnts-1)*ksteps+imagpnts+ksteps)*lambda_steps;

% Construct and combine the parts of the contour.
theta = linspace(0,pi/2,p1);
ln=linspace(R^(1/spread),zerodist^(1/spread),p2).^spread;
preimage=[R*exp(1i*theta) ln(2:end)*1i];

    
    
    
    
    
    