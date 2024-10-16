function D = evans_function(lambda,s,p,options)


% mu_L
[~,df] = f(s.UL(1),p);
a = 1;
b = p.c-Wz(-1,p);
c = @(lambda)df-lambda;
mu_L = @(lambda)(-b+sqrt(b^2-4*a*c(lambda)))/2;

% solve left evans system
muI2 = mu_L(lambda)*eye(2);
ode_evans = @(x,y) (A(x,lambda,s,p)-muI2)*y;
evans_left = ode15s(ode_evans,[s.L,0],[1;mu_L(lambda)],options);

% mu_R
[~,df] = f(s.UR(1),p);
a = 1;
b = p.c-Wz(1,p);
c = @(lambda)df-lambda;
mu_R = @(lambda)(-b-sqrt(b^2-4*a*c(lambda)))/2;

% solve right evans system
muI2 = mu_R(lambda)*eye(2);
ode_evans = @(x,y) (A(x,lambda,s,p)-muI2)*y;
evans_right = ode15s(ode_evans,[s.R,0],[1;mu_R(lambda)],options);

% evan determinant
D = det([deval(evans_left,0),deval(evans_right,0)]);

