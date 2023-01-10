function out = eval_mani(cf,theta1,theta2)

M = size(cf,1);
N = size(cf,2);

out = zeros(4,1);
for j = 1:M
    for k = 1:N
        out = out + squeeze(cf(M-j+1,N-k+1,:))*theta1^(M-j)*theta2^(N-k);
    end
end