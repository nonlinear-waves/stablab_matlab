function out = get_eigenvalues(p,choices,left,right,num_nodes)

if nargin < 5
    num_nodes = 10;
end

s.options = odeset('RelTol',1e-8,'AbsTol',1e-8);

lambda_vals = linspace(left,min(right,sqrt(p.g/p.L0)),num_nodes);

D = zeros(size(lambda_vals));
for j = 1:length(lambda_vals)

    lambda = lambda_vals(j);

    D(j) = local_evans(lambda,s,p);

end

index = find(D(2:end).*D(1:end-1)<0);

out = zeros(size(choices));

for cnt = 1:length(choices)
    choice = choices(cnt);
    a = lambda_vals(index(end-choice+1));
    fa = D(index(end-choice+1));
    b = lambda_vals(index(end-choice+1)+1);
    while (b-a) > 1e-8
        c = 0.5*(a+b);
        fc = local_evans(c,s,p);
        if sign(fc)==sign(fa)
            a = c;
        else
            b = c;
        end
    end
    out(cnt) = c;
end
