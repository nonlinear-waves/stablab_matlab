clc; clear all; beep off;  close all; drawnow;

file_name = mfilename;

cd('../core_code');
curr_dir = cd;

%
% parameters
%

p.gamma = 5/3;
u1_p_vals = linspace(0.9,0.4,1000);


%
% preimage contour
%

pts = 1000;
dom_a = 1e-5;
dom_b = 1;
spread = 4;
preimage = linspace(dom_b^(1/spread),dom_a^(1/spread),pts).^spread;

%
% find transition curve
%

tol = 1e-3;

cnt = 0;
ind = 1;
for j = 1:length(u1_p_vals)

    disp(j)
    p.u1_p = u1_p_vals(j);
    
    if j > 2
        h_guess = 2*data{ind-1}.p.h1-data{ind-2}.p.h1
        if h_guess <= 1
            h_guess = 1.0001;
        end
        h1_low = max(1+0.5*abs(h_guess-1),h_guess-0.1);
        h1_up = h_guess + 0.1;
    else
        h1_up = 3;
        h1_low = 2;
    end

    try
        p.h1 = h1_up;
        [w,prew] = lop_for_transition(p,preimage);
        w_up = w(end);

        p.h1 = h1_low;
        [w,prew] = lop_for_transition(p,preimage);
        w_low = w(end);

        if sign(w_up) == sign(w_low)
           error('signs are not opposite'); 
        end

        h1_mid = 0.5*(h1_up+h1_low);

        while norm(h1_up-h1_low)/norm(h1_low) > tol
            p.h1 = h1_mid;
            [w,prew] = lop_for_transition(p,preimage);
            w_mid = w(end);
            if sign(w_mid) == sign(w_low)
                h1_low = h1_mid;
            else 
                h1_up = h1_mid;
            end
            h1_mid = 0.5*(h1_low+h1_up);

        end

        p.h1 = 0.5*(h1_low+h1_up);

        hold on;
        plot(1-p.u1_p,p.h1-1,'.k','MarkerSize',18);
        drawnow;

        cnt = cnt + 1;
        data{cnt}.p = p;
        
        cd('../data');
        save(file_name,'data');
        cd('../core_code');
        ind = ind + 1;
        
    catch me
        disp(me)
        cd(curr_dir);
    end
end


cd(curr_dir);






