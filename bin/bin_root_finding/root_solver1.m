function rts = root_solver1(box,tol,p,s,e,m,c)
% out = root_solver1(box,p,s,e,m,c)
%
% Subdivides boxes until boxes containing roots are within tolerance.
% The input box = [a b c d] is a box with lower left coordinate (a,b)
% and upper right coordinate (c,d). The input tol is a bound on the
% relative tolerance of the location of the roots. The inputs p, s, e, m,
% and c are the standard STABLAB structures.
%
% Alternatively, if c.moments = 'on', then the program subdivides boxes
% until there is only one root in the box and then computes the moment to
% determine the location of the root. When there is only one root inside a
% box, computing the moments generally provides a good approximation.

% Future improvements: Make it so the program re-uses the contours already
% computed. Also, make it so it doesn't store contours if not necessary.

rts = [];

box_s{1}{1} = cont_box(box,c); % box
box_s{1}{2} = box; % box coordinates
num_boxes = 1; 
location = 1;
num_rts = 0;
box_s{1}{3} = 0; % box number of parent

while location <= num_boxes
    
    % Don't evaluate a box if we already know there are no roots in it
    parent_id = box_s{location}{3};
    if parent_id > 0        
        if box_s{parent_id}{4} == box_s{parent_id}{5}
           location = location + 1;
           continue 
        end
    end
    
    % plot the box that we are currently computing
    if strcmp(c.pic_stats,'on')
        hold on
        plot(get_corners(box_s{location}{2}),'.-m')
        drawnow;
    end
    
    % Evans function call
    [D,domain] = c.root_fun(c,s,p,m,e,box_s{location}{1});
        
    % winding number processing
    wnd = winding_number(D);
    box_s{location}{4} = wnd; % winding number of this box
    box_s{location}{5} = 0; % number of roots found in children of this box
    if parent_id > 0
       box_s{parent_id}{5} = box_s{parent_id}{5} + wnd;
    end
    
    % find roots or subdivide box
    if wnd > 0
 
        % box coordinates
        aa = box_s{location}{2}(1);
        bb = box_s{location}{2}(2);
        cc = box_s{location}{2}(3);
        dd = box_s{location}{2}(4);

        % use the method of moments to find root when wnd = 1
        if isfield(c,'moments')
            if strcmp(c.moments,'on')
                if wnd == 1
                    do_solve = 1;
                    if isfield(c,'moments_tol')
                       if strcmp(c.moments_tol,'on')
                           norm([aa,bb]-[cc,dd])
                           if norm([aa,bb]-[cc,dd]) > tol
                              do_solve = 0; 
                           end
                       end
                    end
                    
                    if do_solve == 1
                        num_rts = num_rts+1;
                        rts(num_rts) = moments_roots(domain,D);
                        if strcmp(c.pic_stats,'on')
                             hold on
                             plot(real(rts),imag(rts),'.k','MarkerSize',18);
                             drawnow;
                        end
                        location = location + 1;
                        continue
                    end
                    
                end
            end
        end
        
        % use the relative size of diagonal of box to find root when wnd = 1
        if norm([aa,bb]-[cc,dd])/min(norm([aa,bb]),norm([cc,dd])) < tol
            num_rts = num_rts+1;
            rts(num_rts) = 0.5*(aa+cc)+0.5*(bb+dd)*1i;
            if strcmp(c.pic_stats,'on')
                 hold on
                 plot(real(rts),imag(rts),'.k','MarkerSize',18);
                 drawnow;
            end
            location = location + 1;
            continue
        end
        
        % if wnd > 1, divide boxes
        num_boxes = num_boxes + 1;
        temp_box = [aa,bb,0.5*(aa+cc),0.5*(bb+dd)];
        box_s{num_boxes}{1} = cont_box(temp_box,c);
        box_s{num_boxes}{2} = temp_box;
        box_s{num_boxes}{3} = location;
        
        num_boxes = num_boxes + 1;
        temp_box = [0.5*(aa+cc),bb,cc,0.5*(bb+dd)];
        box_s{num_boxes}{1} = cont_box(temp_box,c);
        box_s{num_boxes}{2} = temp_box;
        box_s{num_boxes}{3} = location;

        num_boxes = num_boxes + 1;
        temp_box = [0.5*(aa+cc),0.5*(bb+dd),cc,dd];
        box_s{num_boxes}{1} = cont_box(temp_box,c);
        box_s{num_boxes}{2} = temp_box;
        box_s{num_boxes}{3} = location;

        num_boxes = num_boxes + 1;
        temp_box = [aa,0.5*(bb+dd),0.5*(aa+cc),dd];
        box_s{num_boxes}{1} = cont_box(temp_box,c);
        box_s{num_boxes}{2} = temp_box;
        box_s{num_boxes}{3} = location;
  
    end
    location = location + 1;
end


%------------------------------------------------------------
function out = cont_box(box,c)
% make contour out of box coordinates

line1 = linspace(box(1),box(3),10+9*c.ksteps+(10+9*c.ksteps-1)*c.lambda_steps)+box(2)*1i;
line2 = linspace(box(2),box(4),10+9*c.ksteps+(10+9*c.ksteps-1)*c.lambda_steps)*1i+box(3);
line3 = fliplr(linspace(box(1),box(3),10+9*c.ksteps+(10+9*c.ksteps-1)*c.lambda_steps))+box(4)*1i;
line4 = fliplr(linspace(box(2),box(4),10+9*c.ksteps+(10+9*c.ksteps-1)*c.lambda_steps))*1i+box(1);

out = [line1 line2(2:end) line3(2:end) line4(2:end)];
%-------------------------------------------------------------
function out = get_corners(cont)
% prep box coordinates for plotting

out = [cont(1)+1i*cont(2),cont(3)+1i*cont(2),cont(3)+1i*cont(4),cont(1)+1i*cont(4),cont(1)+1i*cont(2)];








