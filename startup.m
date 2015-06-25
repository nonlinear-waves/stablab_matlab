function startup(varargin)

% startup.m
%
% Set the search paths for STABLAB_2.0
%
% OPTIONS
%
% 'intlab',''   > start INTLAB 

disp('synclient HorizEdgeScroll=0 HorizTwoFingerScroll=0');

% rng('shuffle');

eval(['addpath ',pwd,'/bin/PolyfitnTools/PolyfitnTools'])
eval(['addpath ',pwd,'/bin'])
eval(['addpath ',pwd,'/bin/bin_profile'])
eval(['addpath ',pwd,'/bin/bin_root_finding'])
eval(['addpath ',pwd,'/bin/bin_main']);

% default values
global start_matlabpool_at_startup;
start_matlabpool_at_startup = 1;

% evaluate variable input arguments
n = length(varargin);
k = 0;
while k < n
    k = k+1;
    com = varargin{k};
    execute_commands(com,varargin(k+1));
    k = k+1;
end

if start_matlabpool_at_startup == 1
    try 
        try
            delete(gcp('nocreate'));
        catch me
            matlabpool close;
        end
    catch me
    end

    try
        try
            parpool
        catch me
            matlabpool open
        end
    catch me
    end
end
clear start_matlabpool_at_startup;

%%%%%%%%%%%%%%%%%%%
function execute_commands(com,vec)

global start_matlabpool_at_startup;

switch com
    % run start up files for Intlab
    case 'intlab'
        
        cur_dir = cd;
        str = cd;
        ind = strfind(str,'Dropbox/stablab20');
%         if isempty(ind)
%                 ind = strfind(str,'Dropbox/stablab20');
%                 if isempty(ind)
%                     error('The current directory is not conatined in Dropbox/STABLAB/');
%                 end
%         end
        cd([str(1:ind-1),'Dropbox/stablab20/bin/Intlab_V6']);
        ls
        startintlab;
        cd(cur_dir);
        
    % don't try opening the Matlab pool
    case 'start matlabpool'

        if strcmp(vec,'off')
            start_matlabpool_at_startup = 0;
        end
        
end




