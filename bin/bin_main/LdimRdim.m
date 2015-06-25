function out = LdimRdim(Amat,s,p)

out = zeros(2,1);

egs = eig(Amat(s.L,1,s,p));

indx = find(real(egs)>0);

out(1) = length(indx);

indx = find(real(egs)<0);

out(2) = length(indx);
