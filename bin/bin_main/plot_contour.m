function plot_contour(w)
% plot the Evans function


figure;
hold on;

plot(w,'.-k','LineWidth',1,'MarkerSize',18);


plot(0,0,'+r','MarkerSize',10,'LineWidth',2);

h = gca;
set(h,'FontSize',22);
set(h,'LineWidth',2);

h = xlabel('Re(\lambda)');
set(h,'FontSize',22);

h = ylabel('Im(\lambda)');
set(h,'FontSize',22);


h = title('contour');
set(h,'FontSize',22);


