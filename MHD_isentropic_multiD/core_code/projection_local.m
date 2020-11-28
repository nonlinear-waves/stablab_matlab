function [P,Q1] = projection_local(A,x,lambda,s,p,posneg,eps)
% [P,Q1] = projection2(matrix,posneg,eps)
%
% Returns a projector P and spanning set Q1 of the invariant subspace
% associated with the given matrix and specified subspace.
%
% Input "matrix" is the matrix from which the eigenprojection comes,
% "posneg" is 1,-1, or 0 if the unstable, stable, or center space is
% sought. The input eps gives a bound on how small the eigenvalues sought
% can be, which is desirable when a zero mode should be avoided.

% Uses Schur decomposition to get a basis for the generalized eigenspace

% determine if lambda is such that the correct subspaces are obvious

lambda0 = 0.5;
step_size = 0.1;

if real(lambda) > lambda0
    
    matrix = A(x,lambda,s,p);



    [U,T] = schur(matrix,'complex');
    E = ordeig(T);
    k = length(find(posneg*real(E)>eps));
    US = ordschur(U,T,posneg*real(E)>eps);
    Q1 = US(:,1:k);

    [U,T] = schur(-matrix,'complex');
    E = ordeig(T);
    k = length(find(posneg*real(E)>-eps));
    US = ordschur(U,T,posneg*real(E)>-eps);
    Q2 = US(:,1:k);

    R = [Q1 Q2];
    L = inv(R);

    P = zeros(size(matrix));
    for k=1:size(Q1,2)
        P = P + R(:,k)*L(k,:);
    end
    
    return
    
end

% 
% left hand side
%

if posneg == 1

    lam = lambda0+imag(lambda);
    mat = A(x,lam,s,p);



    [R,D] = eig(mat);
    ind = find(real(diag(D))==max(real(diag(D))));
    z = R(:,ind);

    for h = 0:step_size:1
        mat = A(x,h*lambda+(1-h)*lam,s,p);
        [R,D] = eig(mat);
        y = abs(z'*R);
        ind = find(y==max(y));
        z = R(:,ind);

    end

    L = inv(R);
    P = zeros(size(R));

    if posneg == 1
        index = ind;
    elseif posneg == -1
        index = [];
        for k = 1:5
            if ~(k==ind)
                index = [index,k] ;
            end
        end
    end

    for j=index
        P = P + R(:,j)*L(j,:);
    end

    Q1 = P*R(:,index);

    return

end

%
% right hand side
%

steps = 10; 

if posneg == -1
   
   lam = lambda0+imag(lambda);
    mat = A(x,lam,s,p);

    [R,D] = eig(mat);
    diagD = diag(D);
    
    [unused,id] = sort(real(diagD));
    
    ind = id(1:3);
     
    z = R(:,ind);

    for h = linspace(0,1,steps)
        mat = A(x,h*lambda+(1-h)*lam,s,p);
        [R,D] = eig(mat);
        
        
        y1 = abs(z(:,1)'*R);
        ind1 = find(y1==max(y1));
        
%         ind10=ind1
        
%         dis = 2;
%         ind1 = 0;
%         for k = 1:size(R,2)
%             dist = min(norm(z(:,1)-R(:,k)),norm(z(:,1)+R(:,k)));
%             if  dist < dis
%                ind1 = k; 
%                dis = dist;
%             end
%         end
        

        y2 = abs(z(:,2)'*R);
        ind2 = find(y2==max(y2));
        y3 = abs(z(:,3)'*R);
        ind3 = find(y3==max(y3));
        z = R(:,[ind1,ind2,ind3]);

    end


    L = inv(R);
    P = zeros(size(R));


    index = [ind1,ind2,ind3];

    for j=index
        P = P + R(:,j)*L(j,:);
    end

    Q1 = P*R(:,index);
    
    return;
    
end








