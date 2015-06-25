function plot_evans(w,style,type)
% plot the Evans function

if nargin < 2
    style = 'normalized';
end
if nargin < 3
    type = 'evans';
end

figure;
hold on;

if strcmp(style,'normalized')
    plot(w/w(1),'.-k','LineWidth',1,'MarkerSize',18);
else
    plot(w,'.-k','LineWidth',1,'MarkerSize',18);
end

plot(0,0,'+r','MarkerSize',10,'LineWidth',2);

h = gca;
set(h,'FontSize',22);
set(h,'LineWidth',2);

h = xlabel('Re(\lambda)');
set(h,'FontSize',22);

h = ylabel('Im(\lambda)');
set(h,'FontSize',22);

if strcmp(type,'evans')
    h = title('Evans Function');
elseif strcmp(type,'lopatinski')
    h = title('Lopatinski determinant');
end
set(h,'FontSize',22);


