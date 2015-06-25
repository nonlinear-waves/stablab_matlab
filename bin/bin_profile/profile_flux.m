function [p,s] = profile_flux(p,s,s_old)
% [p,s] = profile_flux(p,s)
%
% Solves the profile for the flux formulation. 
% If s.s_old exists, then this
% is used for the initial guess. Otherwise, a tanh
% solution is used for an initial guess. Left and right 
% numerical infinity are expanded as needed to assure
% the end state error is within s.tol, though s.L = s.R
% in this program, which may not be necessary. Uneeded 
% mesh points of the solution are removed to speed up
% interpolation and continuation routines.
%
% The user must include in the input structures the following:
%
% s.phase - a vector of phase conditions
% s.order - a vector indicating the order in which the phase 
%               conditions should be applied. 
% s.F - a function handle to the profile, e.g. s.F = @F,  F(x,y,s,p) = ...
% s.UL, UR - end states of the n-r equations that need to be solved in
%                   flux formulation. UL at - infinity, UR at + infinity
% s.n = n-r in flux formulation (number of profile equations to integrate)
%
% Optional input:
% s.tol - profile endstate maxim absolute error, (defaults to 1e-4)
% s.R_max - maximum allowed interval length on right (defaults to 1000)
% s.L_max - maximum allowed interval length on left (defaults to 1000)

%------------------------------------------------------------
% End states
%------------------------------------------------------------

% end state tolerance
if ~isfield(s,'tol')
    s.tol = 1e-4;
end

% maximum value of R allowed
if ~isfield(s,'R_max')
    s.R_max = 10000;
end

% maximum value of L allowed
if ~isfield(s,'L_max')
    s.L_max = 10000;
end

% check if previous guess is provided
if nargin < 3
    s_old = 'none';
end

% numerical infinity
s.I = 1;
% profile solved on right half domain
s.side=1; 
% array for right hand side
s.rarray=1:1:s.n;
% array for left hand side
s.larray=s.n+1:1:2*s.n; 

% bvp solver projections
AM = s.Flinear(s.UL,p);
s.LM = orth(projection1(AM,-1,0).');
AP = s.Flinear(s.UR,p);
s.LP = orth(projection1(AP,1,0).');

%     Eig_neg  = eig(AM)
%     Eig_pos = eig(AP)
%    
%     LM = s.LM
%     LP = s.LP
    

s.n_phs = s.n-size(s.LM,2)-size(s.LP,2);
if s.n_phs < 1
    disp('Eigenvalues at negative infinity: ');
    disp(eig(AM));
    disp('Eigenvalues at positive infinity: ');
    disp(eig(AP));
    error('profile_flux.m does not solve undercompressive profiles');
end

% bvp tolerances
if ~isfield(s,'bvp_options')
    s.bvp_options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000);
end
% positve numerical infinty

% -----------------------------------------------------------
% solve the profile initially
% -----------------------------------------------------------

[p,s] = profile(p,s,s_old);

% -----------------------------------------------------------
% take out extra mesh points
% -----------------------------------------------------------

% stride = how many points to take out of solution to 
% minimize points in final solution.
stride = 3; 
s.stride = stride;
s_old = s;
mesh = length(s_old.sol.x);
mesh_old = mesh+1;
while mesh < mesh_old
    [p,s] = profile(p,s,s_old);
    s_old = s;
    mesh_old = mesh;
    mesh = length(s_old.sol.x);
end
s.stride = stride;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================================
% **************************************************
%                       profile
% **************************************************
% ================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,s] = profile(p,s,s_old)

%--------------------------------------------------------------------------
% provide initial guess
%--------------------------------------------------------------------------

if (isa(s_old,'struct'))
    
    if isfield(s_old,'solver')
        switch s_old.solver
            case {'bvp4c','bvp5c','bvp5c'}
                solver_type = 'bvp';
            otherwise
                solver_type = 'ode';
                s.stride = 3;
        end
    
        if strcmp(solver_type,'ode')
            pre_guess = @(x)(ode_to_bvp_guess(x,s_old,s));      
        elseif strcmp(solver_type,'bvp')
            pre_guess = @(x)(continuation_guess(x,s_old,s));
        end
    else
        pre_guess = @(x)(continuation_guess(x,s_old,s));
    end
        
        stride = s_old.stride;
        count = 1;
        for j = 1:length(s_old.sol.x)
            if mod(j-1,stride) == 0
                x_dom(count) = s_old.sol.x(j);
                count = count + 1;
            end
        end
        if ~(mod(length(s_old.sol.x)-1,stride) == 0)
           x_dom(count) = s_old.sol.x(end); 
        end
        
        s.I = s_old.I;
        s.L = s_old.L;
        s.R = s_old.R;
  
else
               
        s.I = 1;
        if ~(isfield(s,'R'))
           s.R = 5; 
        end
        s.L = -s.R;

        pre_guess = @(x)(guess(x,s));
        x_dom = linspace(0,1,30);
    
end



%--------------------------------------------------------------------------
% convergence to endstates tolerance
%--------------------------------------------------------------------------

err = s.tol + 1;
while err > s.tol
    pre_bc = @(x,y)(bc(x,y,s));
    pre_ode = @(x,y)(double_F(x,y,s,p));

    solinit = bvpinit(x_dom,pre_guess);
    s.sol = bvp5c(pre_ode,pre_bc,solinit,s.bvp_options);
    
    err1 = max(abs(s.sol.y(s.rarray,end)-s.UR));
    err2 = max(abs(s.sol.y(s.larray,end)-s.UL));
    err  = max(err1,err2);
    
    if isfield(s,'stats')
        if strcmp(s.stats,'on')
            fprintf('Profile boundary error: %4.4g\n',err);
        end
    end
    
%     plot(s.sol.x,s.sol.y,'.-');
%     plot(s.sol.y(1,:),s.sol.y(2,:),'-r')
%     drawnow
    
    if err > s.tol
       s_old = s; 
    end
    if err1 > s.tol
       s.R = 1.1*s.R; 
       s.L = -s.R;
    end
    if err2 > s.tol
       s.L = 1.1*s.L; 
       s.R = -s.L;
    end
    if abs(s.L) > s.L_max
       error(['Could not meet specified tolerance in profile solver without',...
           'exceeding the maximum allowed value of negative infinity.']);
    end
    if abs(s.R) > s.R_max
       error(['Could not meet specified tolerance in profile solver without',...
           'exceeding the maximum allowed value of positive infinity.']);
    end
    if err > s.tol
        pre_guess = @(x)(continuation_guess(x,s_old,s));
        x_dom = s_old.sol.x;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================================
% **************************************************
%                              guess
% **************************************************
% ================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = guess(x,s)
% guess using tanh solution
    
a = 0.5*(s.UL+s.UR);
c = 0.5*(s.UL-s.UR);
out = [a-c*tanh((s.R/s.I)*x);a-c*tanh((s.L/s.I)*x)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================================
% **************************************************
%                                 bc
% **************************************************
% ================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = bc(ya,yb,s)
% out = bc(ya,yb,s,p)
% 
% boundary conditions. We split the problem in half and reflect onto
% the right side

out = [
            ya(s.rarray)-ya(s.larray);               %  matching conditions
            s.LM.' * (yb(s.larray) - s.UL);         % projection at - infinity
            s.LP.'*(yb(s.rarray)-s.UR);             % projection at + infinity
            ya(s.order(1:s.n_phs))-s.phase(s.order(1:s.n_phs));    
         ];
  
