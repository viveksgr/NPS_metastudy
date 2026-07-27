% MAKE_MIX_FIGS  Within-subject mixture figures (5D space).
%   fig_mix_paired.png - paired on-target shift: per subject, loading on the
%       arm they are IN vs the OTHER arm. Downward slope = responsiveness.
%   fig_mix_trait.png  - trait stability: per-subject correlation between the
%       two arms' full profiles, one dot per subject, by study.
setup;
S = load('datamat.mat'); cl_mat = S.cl_mat;
GRP = {'Cognitive','Conditioning','Control','Mindfulness', ...
       'Placebo','Placebo_C','Remifentanil','Social'};

R = run_mixture_withinsubj('COORDS',{'C1','C2','C3','NPS','SIIPS'});
nS = numel(R);
% colour each study by its B-arm intervention
bcol = zeros(nS,3);
for i=1:nS, bcol(i,:) = cl_mat(strcmp(GRP,R(i).B),:); end

% ---------------- Figure 1: paired on-target shift ----------------
figure('Position',[50 50 300*nS 300]);
for i = 1:nS
    subplot(1,nS,i); hold on;
    d = R(i).diagv; o = R(i).offv; n = R(i).n;
    for k = 1:n
        col = [0.7 0.7 0.7];
        if d(k) < o(k), col = [0.93 0.6 0.6]; end   % subjects bucking the trend
        plot([1 2], [d(k) o(k)], '-', 'Color', col, 'LineWidth', 0.6);
    end
    plot([1 2], [mean(d) mean(o)], '-', 'Color', bcol(i,:), 'LineWidth', 3);
    errorbar([1 2], [mean(d) mean(o)], [sem(d) sem(o)], 'Color', bcol(i,:), ...
        'LineWidth', 1.5, 'CapSize', 10, 'linestyle','none');
    plot(1, mean(d), 'o', 'MarkerFaceColor', bcol(i,:), 'MarkerEdgeColor','k','MarkerSize',9);
    plot(2, mean(o), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', bcol(i,:), ...
        'LineWidth',1.5, 'MarkerSize',9);
    set(gca,'XTick',[1 2],'XTickLabel',{'current arm','other arm'}, ...
        'XTickLabelRotation',20,'TickLabelInterpreter','none');
    xlim([0.7 2.3]); ylabel('mixture loading');
    title(sprintf('study %s: %s / %s\nacc=%.2f  p=%.1e', R(i).study, R(i).A, R(i).B, ...
        R(i).accR, R(i).pR), 'Interpreter','none','FontSize',9);
    box off; hold off;
end
sgtitle('Paired on-target shift (5D): loading follows the current arm','FontWeight','bold');
exportgraphics(gcf, 'outputs/fig_mix_paired.png', 'Resolution', 200);

% ---------------- Figure 2: trait stability ----------------
figure('Position',[50 50 520 420]); hold on;
yline(0,'k:','LineWidth',1);
for i = 1:nS
    y = R(i).trr; y = y(~isnan(y)); n = numel(y);
    x = i + (rand(n,1)-0.5)*0.34;
    scatter(x, y, 30, bcol(i,:), 'filled', 'MarkerFaceAlpha',0.65, 'MarkerEdgeColor','k');
    plot([i-0.30 i+0.30], [mean(y) mean(y)], '-', 'Color', bcol(i,:), 'LineWidth', 4);
    text(i, 1.20, sprintf('%s / %s\nr=%+.2f', R(i).A, R(i).B, R(i).traitr), ...
        'HorizontalAlignment','center', 'FontWeight','bold', 'Color', bcol(i,:)*0.75, 'FontSize',9);
end
lbls = arrayfun(@(r) sprintf('study %s', r.study), R, 'UniformOutput', false);
set(gca,'XTick',1:nS,'XTickLabel',lbls,'TickLabelInterpreter','none','FontSize',10);
xlim([0.4 nS+0.6]); ylim([-1.05 1.42]);
ylabel('within-subject profile correlation (arm A vs arm B)');
title({'Trait stability: does a subject keep their mixture across arms?', ...
    'high = trait-dominated (intervention nudges);  low/neg = intervention overwrites'}, ...
    'FontSize',9,'FontWeight','normal');
box off; hold off;
exportgraphics(gcf, 'outputs/fig_mix_trait.png', 'Resolution', 200);

disp('mixture figures written');

function s = sem(x), x=x(~isnan(x)); s = std(x)/sqrt(numel(x)); end
