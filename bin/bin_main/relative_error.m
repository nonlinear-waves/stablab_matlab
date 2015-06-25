function out=relative_error(x)
% out=relative_error(x)
%
% Returns max(|x(j+1)-x(j)|/max(|x(j)|,|x(j+1)|)
%
% Input "x" is a vector whose relative error is sought 

out1=(x(2:end)-x(1:end-1))./x(1:end-1);
out2=(x(2:end)-x(1:end-1))./x(2:end);
out1=max(abs(out1));
out2=max(abs(out2));
out=max(out1,out2);
