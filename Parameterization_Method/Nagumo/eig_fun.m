function out = eig_fun(xi)


out = zeros(2,length(xi));

for j = 1:length(xi)
   out(:,j) = local_eig(xi(j));
end

function out = local_eig(x)

if x <= 0
   out = [exp(2.*x).*(2-2.*tanh(x)-sech(x).^2); 
          2.*exp(2.*x).*(2-2.*tanh(x)-sech(x).^2) + ...
          exp(2.*x).*(-2.*sech(x).^2+2.*sech(x).^2.*tanh(x))];
elseif x > 0
   out = [exp(-2.*x).*(2+2.*tanh(x)-sech(x).^2);
   -2.*exp(-2.*x).*(2+2.*tanh(x)-sech(x).^2) + ...
          exp(-2.*x).*(2.*sech(x).^2+2.*sech(x).^2.*tanh(x))];
end




