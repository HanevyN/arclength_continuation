x = -.25:0.1:1.25;


f = @(x) x- x.^2;
figure(1), clf
plot(x, f(x), 'k.', 'Markersize', 10)
hold on
yline(0,'b', 'LineWidth',1)
xlabel("$\omega$", 'Interpreter','latex', 'FontSize',20)
ylabel("$-\alpha_i$", 'Interpreter','latex', 'FontSize',20)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

plot(x(3:4), f(x(3:4)), 'b', x(13:14), f(x(13:14)) ,'b', 'LineWidth',2 )
plot(0, f(0), 'ksquare', 'MarkerSize',10)
plot( 1, f(1), 'ksquare','Markersize',10)

