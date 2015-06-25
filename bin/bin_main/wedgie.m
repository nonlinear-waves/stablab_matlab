function out = wedgie(M)
% out = wedgie(M)
%
% Take wedgie product of columns of M.

[n,k] = size(M);
temp = sortrows(combnk(1:n,k));
out = zeros(size(temp,1),1);
for j=1:size(temp,1)
    out(j) = det(M(temp(j,:),:));
end