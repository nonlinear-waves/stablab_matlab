function [out, preimage2] = contour(c,s,p,m,e,pre_preimage)
% [ out, preimage2 ] = contour2(c,s,p,m,e,preimage)
%
% Returns the Evans function output for the given input. The structures
% c, s, p, m, and e are described in the STABLAB documentation. The 
% input pre_preimage contains the contour points from which the Evans 
% function will be evaluated. The sctructure c should contain a field 
% lambda_steps and ksteps. The positive integer, c.ksteps indicates how
% many Kato steps will be taken between points on which the Evans function
% is initially evaluated. The positive integer c.lambda_steps indicates how
% many additional points are specified between Kato steps in the contour,
% pre_preimage, for optional evaluation of the Evans function if needed to 
% obtain desired relative error, if specified. 
%
% Example: Suppose we want to evaluate the Evans function on a contour with
% 3 points. In order to get accurate results we determine that we need to
% take 2 Kato steps bewteen each point. Then our preimage contour will have
% entries [ 1 2 3 4 5 6 7] with entries 1, 4, and 7 corresponding to the 3
% points we want to Evaluate the Evans function on and entries 2 and 3 intermediate
% points on our contour between 1 and 4, and points 5 and 6 intermediate points on
% the preimage contour between 4 and 7. If we want our Evans
% fucntion output to vary in relative error bewteen consecutive points by
% less than some tolerance, c.tol, then we set c.refine = 'on'. Now the Evans
% fucntion will be evaluated also on entries 2,3,5, and 6 if needed. Perhaps we
% can achieve relative error only evaluating the Evans function on the
% additional point 3. Then the Evans function is not evaluated on points
% 2, 5, and 6. Suppose now that even after evaluating the Evans function on
% all the preiamge points we don't meet relative tolerance requirements.
% Perhaps the region where we don't meet tolerance is only in one small
% region. We don't want to slow the computation way down by computing the
% Kato basis numerically on additional points on all of the contour. If we
% specify c.lambda_steps = 1, for example, then between each Kato step, we
% can compute the Evans function on an extra point if needed without
% originally computing the Evans function on that point. So now our
% preimage contour has 13 entries with the original points we compute on
% residing in entries 1 7, and 13. The entries on which the Kato basis are
% computed are 1,3,5,7,9,11,13. The entries on which the Kato basis can be
% computed if needed and the Evans fucntion evaluated are 2,4,6,8,10,12.
%
% If desired, one can set c.check = 'on'. Then the Evans fucntion is first
% evaluated on the last entry of the preimage contour to determine if the
% output and the conjugate of the output are within specified relative
% error, c.tol. This can save time computing the Evans fucntion on half
% a contour approaching the origin if the Evans function is going to fail
% at the origin. 
%
% If c.stats = 'on', then a waitbar showing computation 
% progress is displayed. If c.stats = 'print', then computation progress is
% printed to the command window instead of to a waitbar (use this option
% for parallel computing)

if isfield(c,'stats')
    if strcmp(c.stats,'on')
        cstats = 1;
    else
        cstats = 0;
    end
    if strcmp(c.stats,'print')
        disp('Finding the Kato basis');
    end
else
    cstats = 0;
end

% Find the subset on which Kato basis is initially evaluated
pre_index = 1:(c.lambda_steps+1):length(pre_preimage);
preimage = pre_preimage(pre_index);

% Find the subset on which Evans fucntion is initially evaluated
[lbasis,lproj] = c.basisL(c.Lproj,c.L,preimage,s,p,c.LA,1,c.epsl);
[rbasis, rproj] = c.basisR(c.Rproj,c.R,preimage,s,p,c.RA,-1,c.epsr);
index =  1:(c.ksteps+1):length(preimage);
lbasis2 = lbasis(:,:,index);
rbasis2 = rbasis(:,:,index);
preimage2 = preimage(index);

out = zeros(1,length(preimage2));

% Makes sure the contour can be sucessfully computed close enough to the
% origin to satisfy tolerance before computing everything.
if isfield(c,'check')
    if strcmp(c.check,'on')
        try 
            out(1) = c.evans(lbasis2(:,:,1),rbasis2(:,:,1),preimage2(:,1),s,p,m,e);
            near_origin = c.evans(lbasis2(:,:,end),rbasis2(:,:,end),preimage2(:,end),s,p,m,e);
            out(end) = near_origin/out(1);
                        if abs(conj(out(end))-out(end))/abs(out(end)) > c.tol
                            out(1)
                            out(end)
                            abs(conj(out(end))-out(end))/abs(out(end))
                            error('The Evans function does not satisfy tolerance at the endpoint of the contour');
                        end
        catch me
              error('The Evans function failed to compute at the end point.');
        end
    end
end

% Compute the Evans fucntion on the intial contour
if cstats == 1
    h = waitbar(0,'Initial. ');
    time1 = tic;
    for j=1:length(index) 
            out(j) = c.evans(lbasis2(:,:,j),rbasis2(:,:,j),preimage2(:,j),s,p,m,e);
            curr_time = toc(time1);
            time_left = curr_time*(length(index)/j)-curr_time;
            waitbar(j/length(index),h,['Initial. Est time left: ',num2str(time_left), 'sec']);
    end
    close(h);
else
    if isfield(c,'stats')
        if strcmp(c.stats,'print')
            disp('Computing the Evans function on the first set of points')
        end
    end
    if isfield(c,'debug')
        for j=1:length(index) 
                out(j) = c.evans(lbasis2(:,:,j),rbasis2(:,:,j),preimage2(j),s,p,m,e);
        end    
    else
        parfor j=1:length(index) 
                out(j) = c.evans(lbasis2(:,:,j),rbasis2(:,:,j),preimage2(j),s,p,m,e);
        end
    end
end

% check if relative error tolerance has been specified
if isfield(c,'refine')
    if ~(strcmp(c.refine,'on'))
        return
    end
else
    return
end
 
% Refine the mesh on which the Evans function is computed until requested tolerance is
% achieved using the Kato steps as needed.

 rel_error = relative_error(out);
 if isfield(c,'stats')
     if strcmp(c.stats,'print')
         fprintf('\nRelative Error: %4.4g\n',rel_error);
     end
 end
 if rel_error > c.tol 
      c.fail = 'off';
      [preimage2,out,index] = refine_contour(index,preimage,lbasis,rbasis,preimage2,out,s,p,m,c,e,cstats);
 end
 
 if isfield(c,'best_refine')
     if c.lambda_steps == 0
        return; 
     end
 end
 
 % use lambda_steps to further refine the mesh as needed
 rel_error = relative_error(out);
 if isfield(c,'stats')
     if strcmp(c.stats,'print')
         fprintf('\nRelative Error: %4.4g\n',rel_error);
     end
 end
 if rel_error > c.tol
     % test if lambda_steps have been specified
    if c.lambda_steps == 0
        atemp = out/out(1);
        atemp = [atemp(1:end-1),fliplr(conj(atemp))];
%         hold on
%         plot(atemp,'.-k')
%         for jk = 1:length(atemp)-1
%             if abs((atemp(jk+1)-atemp(jk))/atemp(jk))>c.tol
%                 plot(real(atemp(jk)),imag(atemp(jk)),'.r')
%             end
%         end
        error('Not enough lambda_steps specified to obtain desired raltive error');
    end
    
    final_out_index = 0;
    final_out = out;
    final_preimage2 = preimage2;
     
    % compute the Evans function on additional points in the problem areas
    % using lambda_steps
     for j=1:length(out)-1
        if abs(out(j+1)-out(j))/min(abs(out(j)),abs(out(j+1))) > c.tol  

            % Find the needed points, basis, etc. between Kato steps where
            % tolerance is too large
            temp_preimage = ...
                pre_preimage((index(j)-1)*c.lambda_steps+index(j):1:...
                (index(j+1)-1)*c.lambda_steps+index(j+1));
      
            temp_lbasis = c.basisL(c.Lproj,c.L,temp_preimage,s,p,c.LA,1,c.epsl,lbasis(:,:,index(j)),lproj(:,:,index(j)));
            temp_rbasis = c.basisR(c.Rproj,c.R,temp_preimage,s,p,c.RA,-1,c.epsr,rbasis(:,:,index(j)),rproj(:,:,index(j)));
            temp_index = [1,length(temp_preimage)];
            temp_preimage2 = [temp_preimage(1),temp_preimage(end)];
            temp_out = [out(j) out(j+1)];
            
            % compute the Evans function
            if isfield(c,'best_refine')
                if strcmp(c.best_refine,'on')
                   c.fail = 'off';
                else
                    c.fail = 'on';
                end
            else
                c.fail = 'on';
            end
            [temp_preimage2, temp_out, temp_index] = ...
                refine_contour(temp_index,temp_preimage,temp_lbasis,temp_rbasis,temp_preimage2,temp_out,s,p,m,c,e,cstats);
            
            % merge the new computations with old ones
            final_out = [final_out(1:j+final_out_index),temp_out(2:end-1),final_out(j+final_out_index+1:end)];
            final_preimage2 = [final_preimage2(1:j+final_out_index),temp_preimage2(2:end-1),final_preimage2(j+final_out_index+1:end)];
            final_out_index = final_out_index + length(temp_preimage2)-2;
            
        end
     end
 else
     return;
 end

 % assign merged data to ouput variables
 out = final_out;
 preimage2 = final_preimage2;
 
 % ------------------------------------------------------------------------
 % function [preimage2,out,index] = refine_contour(index,preimage,lbasis,rbasis,preimage2,out,s,p,m,c,e)
 % ------------------------------------------------------------------------
 
 function [preimage2,out,index] = refine_contour(index,preimage,lbasis,rbasis,preimage2,out,s,p,m,c,e,cstats)
 % Compute the Evans fucntion on additional points as needed to achieve
 % requested relative error.
 
 break_while = 'false';
 rel_error = c.tol+1;
 while rel_error > c.tol
     k=1;
     
     % find new points that need to be computed
    clear refine_index
    for j=1:length(out)-1
        if abs(out(j+1)-out(j))/min(abs(out(j)),abs(out(j+1))) > c.tol  
            if index(j+1)-index(j) > 1
                refine_index(k) = round(0.5*(index(j+1)+index(j)));
            else
                if strcmp(c.fail,'on')
                    error('Not enough contour points to meet specified tolerance in refine_contour');
                end
                break_while = 'true';
                break
            end    
            k = k+1;
        end
    end
    
    if strcmp(break_while,'true')
        break;
    end
       
    % gather the needed entries
    refine_lbasis = lbasis(:,:,refine_index);
    refine_rbasis = rbasis(:,:,refine_index);
    refine_preimage = preimage(refine_index); 
    
    % compute the Evans function on the needed entries
    clear refine_out
    if cstats == 1
        h = waitbar(0,['Refine.e rel error = ',num2str(rel_error),' > ',num2str(c.tol)]);
        time1 = tic;
        for j=1:length(refine_index)
            refine_out(j) = c.evans(refine_lbasis(:,:,j),refine_rbasis(:,:,j),refine_preimage(:,j),s,p,m,e);
            curr_time = toc(time1);
            time_left = curr_time*(length(refine_index)/j)-curr_time;
            waitbar(j/length(refine_index),h,['Refine. Rel error = ',num2str(rel_error),'. Est time left: ',num2str(time_left),' sec']);
        end
        close(h);
    else
        parfor j=1:length(refine_index)
            refine_out(j) = c.evans(refine_lbasis(:,:,j),refine_rbasis(:,:,j),refine_preimage(:,j),s,p,m,e);
        end
    end
    
   % merge the two vectors to obtain the refined mesh
   [preimage2,out,index]=merge(preimage2,out,index,refine_preimage,refine_out,refine_index);
   rel_error=relative_error(out);
   if isfield(c,'stats')
        if strcmp(c.stats,'print')
            fprintf('\nRelative Error: %4.4g\n',rel_error);
        end
   end
   
%     figure
%     plot(out,'.-k')
%     drawnow;
   
 end
 
 % ------------------------------------------------------------------------
 % 
 % ------------------------------------------------------------------------
 
 
 
 
 
 
 
 
    