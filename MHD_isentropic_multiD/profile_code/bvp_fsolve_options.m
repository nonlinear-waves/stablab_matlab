function options = bvp_fsolve_options(varargin)

% Default options



if ischar(varargin{1})==1

    % 
    % Start standard options
    %
    
    options.algorithm = optimset('Display','off', 'Jacobian','off', ...
        'Algorithm','Levenberg-Marquardt','TolFun',1e-10);
    
    options.N = 2^8; 
    

    
    
    %
    % End standard options
    %
    
    cnt = 0;

else
    options = varargin{1};
    cnt = 1;
end

  
if cnt == length(varargin)
   return 
end
    
    
while cnt < length(varargin)
    options = set_options(options,varargin{cnt+1},varargin{cnt+2});
    cnt = cnt + 2;
end
    
% -------------------------------------------------------------------------
% set options 
% -------------------------------------------------------------------------

function options = set_options(options,property,choice)

switch property

    case 'algorithm_stats'
        
        if strcmp(choice,'on')
            options.algorithm = optimset(options.algorithm,'Display','iter');
        else
            options.algorithm = optimset(options.algorithm,'Display','off');
        end

end










