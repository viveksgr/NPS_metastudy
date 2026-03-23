function plot_scatter_linear(Px,Py,colorIdx)

% figure()
scatter(Px,Py,18,colorIdx,'filled', 'MarkerEdgeColor', 'k')
hold on
P = [ones(size(Px)) Px];
b = P\Py;
plot(Px,P*b)
[cr,p] = corrcoef(Px,Py);

str = sprintf('r= %.3f; p = %.3f',cr(2),p(2));
title(str)

