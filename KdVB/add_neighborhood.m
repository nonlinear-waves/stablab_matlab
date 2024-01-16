function out = add_neighborhood(cf,r)

% add error

out = cf;
if maxi(sup(abs(imag(cf))))>0
    out(1) = out(1)+iv(-r,r)+1i*iv(-r,r);
else
    out(1) = out(1) + iv(-r,r);
end