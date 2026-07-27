% MAKE_DISP_FIG  Predisposition scatter (the first analysis): per study,
% each subject's off-arm loading for class A vs their real on-arm loading.
%   x = P(A | other-arm state)   "latent" A-ness while NOT in A
%   y = P(A | own-arm state)     actual A-ness while in A
% Positive slope => off-arm loading predicts real loading = predisposition.
% Matches corr(pBA,pAA): study 5 r=+0.19, study 7 r=+0.57, study 9 r=+0.22.
setup;
S = load('datamat.mat'); cl_mat = S.cl_mat;
GRP = {'Cognitive','Conditioning','Control','Mindfulness', ...
       'Placebo','Placebo_C','Remifentanil','Social'};

R = run_mixture_withinsubj('COORDS',{'C1','C2','C3','NPS','SIIPS'});
nS = numel(R);

figure('Position',[50 50 330*nS 320]);
for i = 1:nS
    x = R(i).pBA; y = R(i).pAA;              % off-arm (latent) vs on-arm (real), class A
    col = cl_mat(strcmp(GRP,R(i).A),:);
    [r,p] = corr(x,y); n = numel(x);
    subplot(1,nS,i); hold on;
    scatter(x, y, 44, col, 'filled', 'MarkerFaceAlpha',0.7, 'MarkerEdgeColor','k');
    b = polyfit(x,y,1); xx = linspace(min(x),max(x),50);
    star=''; if p<0.001,star='***';elseif p<0.01,star='**';elseif p<0.05,star='*';end
    plot(xx, polyval(b,xx), '-', 'Color', col*0.7, 'LineWidth', 2.5);
    xlabel(sprintf('off-arm loading  P(%s | other arm)', R(i).A),'Interpreter','none');
    ylabel(sprintf('on-arm loading  P(%s | %s)', R(i).A, R(i).A),'Interpreter','none');
    title(sprintf('study %s:  %s\nr=%+.2f %s  (p=%.3f, n=%d)', R(i).study, R(i).A, r, star, p, n), ...
        'Interpreter','none','FontSize',9);
    axis square; box off; grid on; hold off;
end
sgtitle('Predisposition: does off-arm loading predict real on-arm loading?','FontWeight','bold');
exportgraphics(gcf, 'outputs/fig_disposition.png', 'Resolution', 200);
disp('predisposition scatter written');
