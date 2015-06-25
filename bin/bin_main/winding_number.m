function out=winding_number(w)
% out=winding_number(w)
%
% Returns the winding number of a contour w
%
% Input "w" should be a closed contour not passing through zero.
%
% If f is analytic and nonzero at each point of a simple closed positively
% oriented contour C and is meromorphic inside C, then
% winding_number(w)=N0-Np where w=f'(C)./f(C) and N0 and Np are respectively
% the number of zeros and poles of f inscide C (multiplicity included).
%
% The change in the argument between any two points of w should be less than
% Pi for winding_number(w) to be accurate. 

%Computes the winding number of the contour
out=0;
for k=1:length(w)-1    
    if imag(w(k+1))==0 && real(w(k+1))<0
        kp=pi*sign(imag(w(k)));
    else
        kp=imag(log(w(k+1)/norm(w(k+1))));
    end
    if imag(w(k))==0 && real(w(k))<0
        kc=pi*sign(imag(w(k+1)));
    else 
       kc=imag(log(w(k)/norm(w(k))));
    end
    opt1=kp-kc;
    opt2=-(2*pi-abs(opt1))*sign(opt1);
    if min(abs(opt1),abs(opt2))==abs(opt1)
        out=out+opt1;
    else 
        out=out+opt2;
    end
end
out=round(out/(2*pi));
