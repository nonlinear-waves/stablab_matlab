clc; clear all; close all; beep off; curr_dir = cd;
% 
% batch_nums = [58,55,40,41,6,1,60,61,62];


batch_nums = [61,62];

d = {};
cd('../data');
cnt = 0;
for j = 1:length(batch_nums)
    ld = load(['batch',num2str(batch_nums(j)),'.mat']);
    try
        data = ld.data;
    catch me
       data = ld.d; 
    end
    for k = 1:length(data)
       cnt = cnt + 1;
       d{cnt} = data{k};
    end
end
cd(curr_dir);

figure;
hold on;
for k=1:length(d)
    
    
    p = d{k}.p;
    if 1-p.u1_p <= 0.39
        plot(1-p.u1_p,p.h1-1,'.k','MarkerSize',18);
    end

end

h = xlabel('u_1^--u_1^+');
set(h,'FontSize',18);
h = ylabel('h_1-H_*');
set(h,'FontSize',18);
h = gca;
set(h,'FontSize',18);
h = title('Inviscid stability boundary');
set(h,'FontSize',18);
