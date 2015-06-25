 function [mergepre,mergepost,mergeindex]=merge(pre1,post1,index1,pre2,post2,index2)
 % [mergepre,mergepost,mergeindex]=merge(pre1,post1,index1,pre2,post2,index2)
 %
 % Takes as input preimages pre1 and pre2 and corresponding function output
 % output post1 and post2 (which may contain output for more than one
 % parameter as in solving the periodic Evans function), and index1 and
 % index2 describing how the two sets of data should be merged together.
 % Returns the two sets of data merged together.
        
    % determine the number of data sets in output contained in variable
    % post
    [n,sy]=size(post1);
    
    % initialize variables
    mergepre=zeros(n,length(pre1)+length(pre2));
    mergepost=zeros(n,size(mergepre,2));
    mergeindex=zeros(1,size(mergepre,2));
    k=1;
    r=1;
    for j=1:length(mergepre)
        
        % attach the last entries if merging has occurred up to the last
        % value of one variable
        if k > length(index1)
            mergepre=[mergepre(1:j-1) pre2(r:end)];
            mergepost(1:n,:)=[mergepost(:,1:j-1) post2(:,r:end)];
            mergeindex=[mergeindex(1:j-1) index2(r:end)];
            break
        end
        if r > length(index2)
            mergepre=[mergepre(1:j-1) pre1(k:end)];            
            mergepost(1:n,:)=[mergepost(1:n,1:j-1) post1(1:n,k:end)];
            mergeindex=[mergeindex(1:j-1) index1(k:end)];
            break
        end
        
        % merge entries where the variables have intersecting values
        if index1(k) < index2(r)
            mergepre(j)=pre1(k);
            mergepost(1:n,j)=post1(1:n,k);
            mergeindex(j)=index1(k);
            k=k+1;
        else
           mergepre(j)=pre2(r);
           mergepost(1:n,j)=post2(1:n,r);
           mergeindex(j)=index2(r);
           r=r+1;
        end
    end