function [s,e,m,c] = emcset(s,shock_type,eLR,Evan_type,func,compound_func)
% function [e,m,c] = emcset(shock_type,eL,eR,Evan_type)
%
% Sets the values of the STABLAB structures e, m, and c to 
% default values. Takes as input a string, shock_type, which is either
% "front" or "periodic". The input eL and eR are respectively the 
% dimension of the left and right eigenspaces of the Evans matrix.
% The input Evan_type is an optional string. If not specified, Evan_type
% will be assigned the most advantageous polar coordinate method.
% Evan_type has the following options when shock_type = 'front':
%
% reg_reg_polar
% reg_adj_polar
% adj_reg_polar
% reg_adj_compound
% adj_reg_compound
%
% when shock_type = 'periodic', the choices are:
%
% regular_periodic
% balanced_periodic
% balanced_polar_scaled_periodic
% balanced_polar_periodic
% balanced_scaled_periodic

eL = eLR(1);
eR = eLR(2);

if nargin < 4
    Evan_type = 'default';
end

if nargin < 5
    func = @A;
end

if nargin < 6
    compound_func = @Ak;
end

if strcmp(shock_type,'front')
    [e,m,c] = initialize_front(s,eL,eR,Evan_type,func,compound_func);
    s.A = func;
    s.Ak = compound_func;
elseif strcmp(shock_type,'periodic')
    [s,e,m,c] = initialize_periodic(s,eL,eR,Evan_type);
elseif strcmp(shock_type,'lopatinski')
    [e,m,c] = initialize_lopatinski(func,s,shock_type);
elseif strcmp(shock_type,'lopatinski2')
    [e,m,c] = initialize_lopatinski(func,s,shock_type);
else
    error('user must specify which type of traveling wave is being studied');
end

%%%%%%%%%%%%%%%%%%%%%%%
function [e,m,c] = initialize_lopatinski(func,s,shock_type)

e.evans = shock_type;

c.LA = func;
c.RA = func;

c.stats = 'off';
c.refine = 'off';
c.tol = 0.2;
c.ksteps = 2^5;
c.lambda_steps = 0;
c.basisL = @analytic_basis;
c.basisR = @analytic_basis;
c.evans = @evans;

c.epsl= 0.000000;
c.epsr = 0.000000;
c.Lproj = @projection5;
c.Rproj = @projection5;

m=[];

c.stats = 'off';
c.refine = 'off';
c.tol = 0.2;
c.ksteps = 2^5;
c.lambda_steps = 0;
c.basisL = @analytic_basis;
c.basisR = @analytic_basis;
c.evans = @evans;

c.epsl= 0.000;
c.epsr = 0.000;
c.Lproj = @projection5;
c.Rproj = @projection5;

%dependent structure variables
e.Li = [s.L 0];
e.Ri = [s.R 0];
c.L = s.L;
c.R = s.R;

%%%%%%%%%%%%%%%%%%%%%%%
function [e,m,c] = initialize_front(s,kL,kR,Evan_type,func,compound_func)

m.n = kL+kR;

if strcmp(Evan_type,'default')
    if kL > m.n/2
        e.evans = 'adj_reg_polar';
    elseif kL < m.n/2
        e.evans='reg_adj_polar';
    else
        e.evans='reg_reg_polar';
    end
else
    e.evans = Evan_type;
end

if strcmp(e.evans,'reg_adj_polar')
    c.LA = func;
    e.LA = c.LA;
    c.RA = @Aadj;
    e.RA = c.RA;
    e.kl = kL;
    r = m.n-kR;
elseif strcmp(e.evans,'reg_reg_polar')
    c.LA = func;
    e.LA = c.LA;
    c.RA = func;
    e.RA = c.RA;
    e.kl = kL;
    e.kr = kR;
elseif strcmp(e.evans,'adj_reg_polar')
    c.LA = @Aadj;
    e.LA = @Aadj;
    c.RA = func;
    e.RA = func;
    e.kl = m.n-kL;
    e.kr = kR;
elseif strcmp(e.evans,'reg_adj_compound')
    c.LA = func;
    e.LA = compound_func;
    c.RA = @Aadj;
    e.RA = @Akadj;
    e.kl = kL;
    e.kr = kR;
elseif strcmp(e.evans,'adj_reg_compound')
    c.LA = @Aadj;
    e.LA = @Akadj;
    c.RA = func;
    e.RA = compound_func;
    e.kl = kL;
    e.kr = kR;
elseif strcmp(e.evans,'reg_reg_cheby')
    c.LA = func;
    e.LA = c.LA;
    c.RA = func;
    e.RA = c.RA;
    e.kl = kL;
    e.kr = kR;
    e.NL = 60;
    e.NR = 60;
end

c.stats = 'off';
c.refine = 'off';
c.tol = 0.2;
c.ksteps = 2^5;
c.lambda_steps = 0;
c.basisL = @analytic_basis;
c.basisR = @analytic_basis;
c.evans = @evans;

c.epsl= 0.000000;
c.epsr = 0.000000;
c.Lproj = @projection5;
c.Rproj = @projection5;

m.damping = 0;
m.method = @drury;
m.options = odeset('RelTol',1e-6,'AbsTol',1e-8,'Refine',1,'Stats','off');
m.ode_fun = @ode15s;

%dependent structure variables
e.Li = [s.L 0];
e.Ri = [s.R 0];
c.L = s.L;
c.R = s.R;

%%%%%%%%%%%%%%%%%%%%%%%%
function [s,e,m,c] = initialize_periodic(s,eL,eR,Evan_type)

%
% find center of the wave
%

xdom = linspace(s.sol.x(1),s.sol.x(end),1000);
yran = zeros(length(xdom),1);
for j = 1:length(xdom)
    temp = deval(s.sol,xdom(j));
    yran(j) = temp(1,1);
end
% plot(xdom,yran,'-k')
% return
maxval = yran(1);
xind = 1;
for j = 1:length(yran)
    if yran(j) > maxval
        maxval = yran(j);
        xind = j;
    end
end
s.center = xdom(xind);

% 
% set default structure values
%

n = eL+eR;

c.ksteps = 2^18;
c.lambda_steps = 0;
c.refine = 'off';
c.tol = 0.2;
c.evans = @evans;

e.A = @Aper;
e.Li = [-s.X/2,0];
e.Ri = [s.X/2,0];
e.kl=n;
e.kr=n;

e.dim_eig_L = eL;
e.dim_eig_R = eR;
s.A = e.A;

m.options = odeset('AbsTol',10^(-13), 'RelTol',10^(-13));
m.A = e.A;
m.method = @drury;
m.damping = 0;
m.n = 2*n;

if strcmp(Evan_type,'regular_periodic')
    e.evans='regular_periodic';
    e.Li = [0,s.X];
    e.kl = n;
elseif strcmp(Evan_type,'balanced_scaled_periodic')
    e.evans ='balanced_scaled_periodic' ;
elseif strcmp(Evan_type,'balanced_periodic')
    e.evans = 'balanced_periodic';
elseif strcmp(Evan_type,'balanced_polar_periodic')
    e.evans = 'balanced_polar_periodic';
elseif strcmp(Evan_type,'balanced_polar_scaled_periodic')
    e.evans ='balanced_polar_scaled_periodic';
elseif strcmp(Evan_type,'default')
    e.evans = 'balanced_polar_periodic';
end

m.ode_fun = @ode15s;








