function plot_profile(s)


hold on;
dom = linspace(s.L,s.R,4000);
y = zeros(2,length(dom));
for j = 1:length(dom)
    y(:,j) = soln(dom(j),s);
end
plot(dom,y,'-k','LineWidth',2);
h = xlabel('x');
set(h,'FontSize',18);
h = ylabel("u,u'");
set(h,'FontSize',18);
h = gca;
set(h,'FontSize',18);