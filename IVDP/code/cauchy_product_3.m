function R = cauchy_product_3(P1, P2, P3, n1, n2)

    sum = 0;
    for k1 = 0:n1
        for k2 = 0:n2 
            for j1 = 0:k1
                for j2 = 0:k2
                    sum = sum + P1(n1-k1+1, n2-k2+1)* ...
                        P2(k1-j1+1, k2-j2+1) ...
                        * P3(j1+1, j2+1);
                end
            end
        end
    end
        
    %removing higher and lower order terms
     sum = sum - P1(n1+1, n2+1)* P2(1, 1)* P3(1, 1)- P1(1, 1) ...
         * P2(n1+1, n2+1)* P3(1, 1)-...
          P1(1, 1)* P2(1, 1)* P3(n1+1, n2+1);

     R = sum;

        
end
     
