function [outFunction] = optimize(inFunction)
%OPTIMIZE Optimizes a function passed in by saving parameters.
%   It takes in a function and returns a function.  When the returned
%   function is called, it records the parameters and saves the output to a
%   file.  Thus, optimized functions will perform faster by loading
%   previously saved values.

    outFunction = @ (varargin) optimizeHelper(inFunction, varargin);
end


function solution = optimizeHelper(varargin)
    warning('off','all');
    inFunction = varargin{1};
    parameters = varargin{2};
    parameterString = '';
    stopIndex = 0;
    %parameters
    %parameters(1:1)
    
    % Get the file name.
    for i=1:size(parameters,2)
        temp = parameters{i};

        if isstring(temp) || ischar(temp)
            parameterString = strcat(parameterString,temp);
            if strcmp(temp,'stop')
                stopIndex = i;
            end
        else
            if isa(temp, 'double') && size(temp,1) == 1 && size(temp,2) == 1
                parameterString = strcat(parameterString , num2str(temp, '%.2f'));

            end
        end
 
        parameterString = strcat(parameterString , '__');
    end
    
    fileStr = strcat(pwd,'/generated/', func2str(inFunction), '__', ...
        parameterString);
    fileStr = strcat(strrep(fileStr, '.','p'),'.mat');
    
    % If the file already exists
    if exist(fileStr, 'file') == 2
        fprintf(strcat('Loading file from:  ', strrep(fileStr,pwd,'')));
        ld = load(fileStr);
        solution = ld.solution;
        
    % If it doesn't exist
    else
        fprintf(' Generating and Saving');
        if stopIndex == 0
            solution = inFunction(parameters{:});
        else
            solution = inFunction(parameters{1:stopIndex-1});
        end
        save(fileStr, 'solution');
    end
    
    warning('on','all');
end