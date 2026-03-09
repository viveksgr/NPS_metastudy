function plot_scatter_linear(Px,Py,colorIdx)

figure()
scatter(Px,Py,36,colorIdx,'filled', 'MarkerEdgeColor', 'k')
hold on
P = [ones(size(Px)) Px];
b = P\Py;
plot(Px,P*b)
cr = corrcoef(Px,Py);
[~,p]=r2t(cr(2),length(Px));
str = sprintf('r= %.3f; p = %.3f',cr(2),p);
title(str)

