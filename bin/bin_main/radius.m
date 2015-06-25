function R = radius(c,s,p,m,e)
%R = radius(c,s,p,m,e)
%
% Find radius of a semicircle enclosing possible unstable eigenvalues of
% D(lambda) for rNS. Fit with $C_1 e ^{C_2\sqrt{\lambda}}$

c.lambda_steps = 0;
c.stats = 'off';
c.refine='off';
err=1;
mult=1;
prev_err = 10^(-12);
while err > c.tol
    R=2^mult;

    mult=mult+1;
    if R > c.max_R
        error('Needed radius for convergence exceeds tolerance');
    end
    
    line=linspace(R,min(2*R,R+50),2+c.ksteps);
    D=contour(c,s,p,m,e,line);
    D=real(D);
    
    % Curve fitting with C1*exp(C2*sqrt(lambda)).
    C2 = (log(D(end))-log(D(1)))/(sqrt(line(end))-sqrt(line(1)));
    C1 = D(end)*exp(-C2*sqrt(line(end)));
    
    theta = linspace(0,pi/2,2+c.ksteps);
    curve = R*exp(1i*theta);
    
    DRi  = contour(c,s,p,m,e,curve);
    err1 = abs(DRi(end)-C1*exp(C2*sqrt(curve(end))))/abs(DRi(end));
    err2 = abs(conj(DRi(end))-C1*exp(C2*sqrt(curve(end))))/abs(DRi(end));
    
    err = min(err1,err2);

%     fprintf('Radius solve: R = %4.4g, error = %4.4g, error convergence = %4.4g\n',R,err,(err/prev_err));
    if err/prev_err > 0.8
       prev_err = err;
       err = 1;
    else
       prev_err = err;
    end
   
    % Compute and compare the fit on more points of the semicircle
    if err <= c.tol
        theta=linspace(0,pi/2,7*(c.ksteps)+8);
        curve=R*exp(1i*theta);
        DR=contour(c,s,p,m,e,curve);
        index=1:(c.ksteps+1):length(curve);
        curve2=curve(index);
        half_err=0;
        for j=1:length(DR)
            err1=abs(DR(j)-C1*exp(C2*sqrt(curve2(j))))/abs(DR(j));
            err2=abs(conj(DR(j))-C1*exp(C2*sqrt(curve2(j))))/abs(DR(j));
            half_err=max(half_err,min(err1,err2));
        end
        
        if half_err > c.tol
            err=1;
        end
    end
end
























% % % function R = radius(c,s,p,m,e)
% % % %R = radius(c,s,p,m,e)
% % % %
% % % % Find radius of a semicircle enclosing possible unstable eigenvalues of
% % % % D(lambda) for rNS. Fit with $C_1 e ^{C_2\sqrt{\lambda}}$
% % % 
% % % stats = c.stats;
% % % c.stats = 'off';
% % % c.refine='off';
% % % err=1e10;
% % % mult=1;
% % % prev_err = 10^(-12);
% % % while err > c.Rtol
% % % R=2^mult;
% % % 
% % %     mult=mult+1;
% % %     if R > c.max_R
% % %         error('Needed radius for convergence exceeds tolerance');
% % %     end
% % %     
% % %     line=linspace(R,min(2*R,R+50),2+c.ksteps);
% % %     D=contour(c,s,p,m,e,line);
% % %     D=real(D);
% % %     D = D/D(1);
% % %     
% % %     % Curve fitting with C1*exp(C2*sqrt(lambda)).
% % %     C2 = (log(D(end))-log(D(1)))/(sqrt(line(end))-sqrt(line(1)));
% % %     C1 = D(end)*exp(-C2*sqrt(line(end)));
% % %    
% % %     theta = linspace(0,pi/2,2+c.ksteps);
% % %     curve = R*exp(1i*theta);
% % %     
% % %     DRi  = contour(c,s,p,m,e,curve);
% % %     DRi = DRi/DRi(1);
% % %     err1 = abs(DRi(end)-C1*exp(C2*sqrt(curve(end))))/abs(DRi(end));
% % %     err2 = abs(conj(DRi(end))-C1*exp(C2*sqrt(curve(end))))/abs(DRi(end));
% % %     
% % %     err = min(err1,err2);
% % % 
% % %     if strcmp(stats,'on')
% % %         fprintf('Radius solve: R = %4.4g, error = %4.4g, error convergence = %4.4g, C1 = %4.4g, C2 = %4.4g\n',R,err,(err/prev_err),C1,C2);
% % %     end
% % %     if err/prev_err > 1e10
% % %        prev_err = err;
% % %        err = 1;
% % %     else
% % %        prev_err = err;
% % %     end
% % %    
% % %     % Compute and compare the fit on more points of the semicircle
% % %     if err <= c.Rtol
% % %         theta=linspace(0,pi/2,7*(c.ksteps)+8);
% % %         curve=R*exp(1i*theta);
% % %         DR=contour(c,s,p,m,e,curve);
% % %         DR = DR/DR(1);
% % %         index=1:(c.ksteps+1):length(curve);
% % %         curve2=curve(index);
% % %         half_err=0;
% % %         for j=1:length(DR)
% % %             err1=abs(DR(j)-C1*exp(C2*sqrt(curve2(j))))/abs(DR(j));
% % %             err2=abs(conj(DR(j))-C1*exp(C2*sqrt(curve2(j))))/abs(DR(j));
% % %             half_err=max(half_err,min(err1,err2));
% % %         end
% % %         
% % %         
% % %         guess = C1*exp(C2*sqrt(curve));
% % % %         guess2 = conj(C1*exp(C2*sqrt(curve)));
% % %         clf;
% % %         hold on;
% % %         plot(guess,'.-k')
% % %         plot(DR,'o-r');
% % % %         plot(guess2,'.-b');
% % % 
% % %     
% % %         if half_err > c.Rtol
% % %             err=1;
% % %         end
% % %     end
% % % end





















