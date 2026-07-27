function plotTest2Acc(R, cl_mat, ttl, varargin)
% plotTest2Acc  Bar plots of Test-2 decoding accuracy from the full space,
%               the readout axis alone, and the null space alone.
%
%   plotTest2Acc(R, cl_mat, ttl)                     % overall: 3 bars
%   plotTest2Acc(R, cl_mat, ttl, 'byGroup', true)    % per intervention
%
%   R      - struct from run_test2 (needs acc, perm, boot, accG, pG, bootG)
%   cl_mat - Nx3 RGB, one row per intervention (GRP order), used in byGroup mode.
%
%   Error bars are bootstrap-over-studies 95% intervals. Asterisks mark
%   significance against the PERMUTATION FLOOR (group labels shuffled). The
%   dashed line is that floor, i.e. the empirical chance level for this
%   pipeline - unlike the Test-1 isotropy line it is a genuine floor that
%   bars are expected to sit above.
%
%   NB bar heights are NOT dimension-normalized: full uses k dims, null k-1,
%   readout 1. Comparing null vs readout by height therefore conflates
%   subspace identity with dimensionality - see the rand-1d / rand-(k-1)d
%   controls in run_test2 for the dimension-matched comparison.
%
%   Styling follows common_functions/plotFixedEffects.m.

    p = inputParser;
    addParameter(p, 'byGroup', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    byGroup = p.Results.byGroup;

    sets   = R.sets;
    labels = {'full (z)','readout (a)','null (z^\perp)'};
    GRP    = R.GRP;
    nS     = numel(sets);

    if ~byGroup
        mu = zeros(nS,1); lo = zeros(nS,1); hi = zeros(nS,1); pv = zeros(nS,1);
        for i = 1:nS
            mu(i) = R.acc.(sets{i});
            ci    = prctile(R.boot.(sets{i}), [2.5 97.5]);
            lo(i) = mu(i)-ci(1); hi(i) = ci(2)-mu(i);
            pv(i) = (1+sum(R.perm.(sets{i}) >= mu(i)))/(numel(R.perm.(sets{i}))+1);
        end
        flo = mean(R.perm.full);

        figure('Position',[0.5 0.5 430 300]); hold on;
        cols = [0.35 0.35 0.35; 0.80 0.45 0.15; 0.20 0.40 0.70];
        for i = 1:nS
            bar(i, mu(i), 'FaceColor', cols(i,:), 'EdgeColor','k');
        end
        errorbar(1:nS, mu, lo, hi, 'k','linestyle','none','LineWidth',1.2,'CapSize',8);
        hF = plot([0.4 nS+0.6], [flo flo], 'k--', 'LineWidth', 1.3);

        yr = max(mu+hi);
        for i = 1:nS
            s = pStars(pv(i));
            if isempty(s), continue; end
            text(i, mu(i)+hi(i)+0.03*yr, s, 'Color','r','HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom','FontWeight','bold');
        end

        set(gca,'XTick',1:nS,'XTickLabel',labels,'XTickLabelRotation',45, ...
                'TickLabelInterpreter','tex');
        xlim([0.4 nS+0.6]); ylim([0 yr*1.20]);
        ylabel('balanced accuracy');
        legend(hF,'permutation floor','Location','northoutside','Box','off');
        title(ttl,'Interpreter','none'); box off; hold off;

    else
        n = numel(GRP);
        keep = find(any(~isnan(R.accG),2));
        nk = numel(keep);

        figure('Position',[0.5 0.5 680 350]); hold on;
        w = 0.26;
        offs = ([1 2 3]-2)*w;                     % three bars per intervention
        shade = [0.30 0.60 1.00];                 % full=lightest -> null=saturated
        h = gobjects(nS,1);
        for q = 1:nk
            i = keep(q); col = cl_mat(i,:);
            for s = 1:nS
                fc = 1-shade(s)*(1-col);
                hb = bar(q+offs(s), R.accG(i,s), w, 'FaceColor', fc, 'EdgeColor','k');
                if q==1, h(s) = hb; end
            end
        end

        yrAll = [];
        for s = 1:nS
            mu = R.accG(keep,s);
            ci = prctile(R.bootG.(sets{s})(:,keep), [2.5 97.5], 1)';
            lo = mu-ci(:,1); hi = ci(:,2)-mu;
            lo(isnan(lo))=0; hi(isnan(hi))=0;
            errorbar((1:nk)'+offs(s), mu, lo, hi, 'k','linestyle','none', ...
                     'LineWidth',1.0,'CapSize',4);
            yrAll = [yrAll; mu+hi]; %#ok<AGROW>
        end
        yr = max(yrAll);
        for q = 1:nk
            i = keep(q);
            for s = 1:nS
                st = pStars(R.pG(i,s));
                if isempty(st), continue; end
                mu = R.accG(i,s);
                ci = prctile(R.bootG.(sets{s})(:,i), 97.5);
                if isnan(ci), ci = mu; end
                text(q+offs(s), ci+0.02*yr, st, 'Color','r', ...
                     'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                     'FontWeight','bold','FontSize',8);
            end
        end

        flo = mean(R.perm.full);
        plot([0.4 nk+0.6],[flo flo],'k--','LineWidth',1.2);

        set(gca,'XTick',1:nk,'XTickLabel',GRP(keep),'XTickLabelRotation',45, ...
                'TickLabelInterpreter','none');
        xlim([0.4 nk+0.6]); ylim([0 yr*1.22]);
        ylabel('balanced accuracy (per-group recall)');
        legend(h, labels, 'Location','northoutside','Orientation','horizontal', ...
               'Box','off','Interpreter','tex');
        title(ttl,'Interpreter','none'); box off; hold off;
    end
end

function s = pStars(p)
    if     isnan(p),  s = '';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = '';
    end
end
