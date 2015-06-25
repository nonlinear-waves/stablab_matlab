function phase_portrait_3D(p,s,rest_points,ode,box,jacobian,num,time)

minbox=1;
width=1;
[sx,sy]=size(rest_points);

for k = 1:sx
    rp = rest_points(k,:);
    linearization = jacobian([rp(1);rp(2);rp(3)],p);
    [eigenvec,eigenval] = eigs(linearization);
    
    for j=1:3    
        if real(eigenval(j,j))>0
             tspan=[0,time];
        elseif real(eigenval(j,j))<0
            tspan=[0,-time];
        else
            tspan=[0,time];
        end
        hold on
        col=[0,0,0];
        
        plot3(rp(1),rp(2),rp(3),'.-k','MarkerSize',18);
        ynot=10^(-5)*real(eigenvec(:,j))+rp.';
        plot_this(ynot,tspan,minbox,width,ode,box,p,s,col)

        ynot=-10^(-5)*real(eigenvec(:,j))+rp.';
        plot_this(ynot,tspan,minbox,width,ode,box,p,s,col)
    end
    
end

if num > 0
    plot_grid(tspan,ode,box,p,s,num,width)
end

% xlabel('u')
% ylabel('u_x')
% zlabel('u_{xx}')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_grid(tspan,ode,box,p,s,num,width)

    
    for j=1:num
        for k=1:num
            for s=1:num
                ynot=[box(1,1)+(box(1,2)-box(1,1))*(j/(num+1)),...
                    box(2,1)+(box(2,2)-box(2,1))*(k/(num+1)),...
                    box(3,1)+(box(3,2)-box(3,1))*(s/(num+1))];
                options=odeset('AbsTol',10^(-6), 'RelTol',10^(-8));
                sol=ode45(@init,tspan,ynot,options,s,p,ode,box);
                plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:),'-g','LineWidth',width);
                sol=ode45(@init,-tspan,ynot,options,s,p,ode,box);
                plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:),'-g','LineWidth',width);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_this(ynot,tspan,minbox,width,ode,box,p,s,col)
    
    options=odeset('AbsTol',10^(-6), 'RelTol',10^(-8));
    sol=ode45(@init,tspan,ynot,options,s,p,ode,box);
    plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:),'-k','LineWidth',width,'Color',col);    
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=init(x,y,s,p,ode,box)

out=ode(x,y,s,p);

if y(1) > box(1,2)
    out=0;
elseif y(1) < box(1,1)
    out=0;
end
        
if y(2) > box(2,2)
    out=0;
elseif y(2) < box(2,1)
    out=0;
end

if y(3) > box(3,2)
    out=0;
elseif y(3) < box(3,1)
    out=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     

     
     
     
     
     