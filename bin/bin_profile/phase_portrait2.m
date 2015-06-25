function phase_portrait2(ode,jacobian,restpnt,box,time,plot_more,num1,num2,p,s)
% phase_portrait(ode_in,jacobian,restpnt,box_in,time,plot_more,num1,num2)
%
% Plots the phase portrait for the ode given by the function handle ode.
% The jacobian is a function handle to the Jacobian of the system. The input restpnt
% contains the rest points of the sytem. The input time is the length of
% time the ode solver computes. If more than the connections between the
% rest point is desired, plot_more can be set to 'on' and num1 and num2
% specified to break the bounding box into an evenly spaced mesh through
% which solutions in phase space pass. The input box is used in the command
% axis(box), so should be input of that type.

pre_init=@(x,y)(init2(x,y,box,ode,s,p));

minbox=min(box(2)-box(1),box(4)-box(3));

width=3;
hold on
axis(box);
options = odeset('AbsTol',10^(-8), 'RelTol',10^(-8));
[vs,hs]=size(restpnt);
for j=1:vs
    linearization = jacobian([restpnt(j,1),restpnt(j,2)],p);
    [eigenvec,eigenval]=eigs(linearization);
    
    skip = 0;
    if real(eigenval(1,1))>0
        tspan=[0,time];
    elseif real(eigenval(1,1))<0
        tspan=[0,-time];
    else
        skip = 1;
    end
    
    if skip == 0;
        ynot=real(10^(-5)*eigenvec(:,1)+restpnt(j,:)');
        plot_this(ynot,tspan,minbox,width,pre_init,box)

        ynot=real(-10^(-5)*eigenvec(:,1)+restpnt(j,:)');
        plot_this(ynot,tspan,minbox,width,pre_init,box)
    end
    
    skip = 0;
    if real(eigenval(2,2))>0
        tspan=[0,time];
    elseif real(eigenval(2,2))<0
        tspan=[0,-time];
    else
        skip = 1;
    end
    
    if skip == 0
        ynot=real(10^(-5)*eigenvec(:,2)+restpnt(j,:)');
        plot_this(ynot,tspan,minbox,width,pre_init,box)

        ynot=real(-10^(-5)*eigenvec(:,2)+restpnt(j,:)');
        plot_this(ynot,tspan,minbox,width,pre_init,box)
    end
       
end

for j=1:vs
    plot(restpnt(j,1),restpnt(j,2),'.k','MarkerSize',36);
end

if strcmp(plot_more,'on')
    options=odeset('AbsTol',10^(-8), 'RelTol',10^(-8));
    for j=1:num1
        for k=1:num2
            x = box(1)+(1/(num1+1))*(box(2)-box(1))*j;
            y = box(3)+(1/(num2+1))*(box(4)-box(3))*k;
            ynot=[x,y];   
%             plot(ynot(1),ynot(2),'.g','MarkerSize',18);
            sol=ode45(pre_init,[0,time],ynot,options);
            plot(sol.y(1,:),sol.y(2,:),'-k','LineWidth',1);
%             plot(sol.y(1,end),sol.y(2,end),'.b');
            sol=ode45(pre_init,[0,-time],ynot,options);
            plot(sol.y(1,:),sol.y(2,:),'-k','LineWidth',1);
%             plot(sol.y(1,end),sol.y(2,end),'.r');
            
            if isfield(s,'arrow')
                if strcmp(s.arrow,'on')
                    vector_field(x,y,p,pre_init,box)
                end
            end
    
        end
    end
end

h = gca;
set(h,'FontSize',18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_this(ynot,tspan,minbox,width,ode,box)
   
    options=odeset('AbsTol',10^(-8), 'RelTol',10^(-8));
    %sol=ode45(@init,tspan,ynot,options,box,ode);
    sol = ode45(ode,tspan,ynot,options);
    plot(sol.y(1,:),sol.y(2,:),'-k','LineWidth',width);

    M=length(sol.x);
    max_dist=10*norm([sol.y(1,1),sol.y(2,1)]-[sol.y(1,M),sol.y(2,M)]);
    index=M;
    for k=2:M
        dist1=norm([sol.y(1,k),sol.y(2,k)]-[sol.y(1,M),sol.y(2,M)]);
        dist2=norm([sol.y(1,k),sol.y(2,k)]-[sol.y(1,1),sol.y(2,1)]);
        temp = (dist1+1)^2+(dist2+1)^2;
        if temp < max_dist
            max_dist=temp;
            index=k;            
        end
    end

%      if tspan(end)<tspan(1)
%          arrow([sol.y(1,index),sol.y(2,index)],[sol.y(1,index-1),sol.y(2,index-1)]);
%      else
%          arrow([sol.y(1,index-1),sol.y(2,index-1)],[sol.y(1,index),sol.y(2,index)]);
%      end


     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vector_field(x,y,p,ode,box)

            vec = ode(0,[x,y]);
            
            
            scl = 0.1;
            

            arrow([x,y],[x,y]+scl*vec.');

     




