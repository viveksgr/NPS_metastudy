function plotTest1AB(R, cl_mat, ttl, varargin)
% plotTest1AB  Bar plots of the Test-1 separation split: on-axis (a) vs
%              null-space (b) squared separation, with SE across pairs.
%
%   plotTest1AB(R, cl_mat, ttl)                      % overall: 2 bars
%   plotTest1AB(R, cl_mat, ttl, 'byGroup', true)     % per intervention
%
%   R      - struct returned by run_test1 (needs pairs, permA/permB, pAg/pBg)
%   cl_mat - Nx3 RGB, one row per intervention (GRP order), used in byGroup mode.
%
%   Asterisks mark significance against the PERMUTATION NOISE FLOOR (labels
%   shuffled), i.e. "is this component larger than sampling noise produces".
%   NB this figure therefore speaks to whether a and b are individually real,
%   NOT to whether the null space is enriched relative to the readout axis -
%   that comparison needs the isotropy/per-dimension reference (see run_test1
%   benchmark 2), which is deliberately not drawn here.
%
%   Styling follows common_functions/plotFixedEffects.m.

    p = inputParser;
    addParameter(p, 'byGroup', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'GRP', {});
    parse(p, varargin{:});
    byGroup = p.Results.byGroup;
    GRP = p.Results.GRP;
    if isempty(GRP)
        if isfield(R,'GRP'), GRP = R.GRP;
        else, GRP = {'Cognitive','Conditioning','Control','Mindfulness', ...
                     'Placebo','Placebo_C','Remifentanil','Social'};
        end
    end

    P = R.pairs;
    aAll = [P.a]'; bAll = [P.b]';

    if ~byGroup
        % ---------- overall: mean a vs mean b ----------
        mu = [mean(aAll), mean(bAll)];
        se = [sem(aAll),  sem(bAll)];
        pv = [getfielddef(R,'pA',NaN), getfielddef(R,'pB',NaN)];

        figure('Position',[0.5 0.5 400 280]); hold on;
        cols = [0.45 0.45 0.45; 0.20 0.40 0.70];
        for i = 1:2
            bar(i, mu(i), 'FaceColor', cols(i,:), 'EdgeColor','k');
        end
        errorbar(1:2, mu, se, 'k', 'linestyle','none', 'LineWidth',1.2, 'CapSize',8);

        yr = max(mu+se);
        for i = 1:2
            s = pStars(pv(i));
            if isempty(s), continue; end
            text(i, mu(i)+se(i)+0.04*yr, s, 'Color','r', ...
                 'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                 'FontWeight','bold');
        end

        set(gca,'XTick',1:2,'XTickLabel',{'on-axis (a)','null-space (b)'}, ...
                'XTickLabelRotation',45,'TickLabelInterpreter','none');
        xlim([0.4 2.6]); ylim([0 yr*1.20]);
        ylabel('squared separation');
        title(ttl,'Interpreter','none'); box off; hold off;

    else
        % ---------- per intervention ----------
        n = numel(GRP);
        mA = nan(n,1); sA = nan(n,1); mB = nan(n,1); sB = nan(n,1); cnt = zeros(n,1);
        for i = 1:n
            m = strcmp({P.g1}, GRP{i}) | strcmp({P.g2}, GRP{i});
            if ~any(m), continue; end
            cnt(i) = sum(m);
            mA(i) = mean(aAll(m)); sA(i) = sem(aAll(m));
            mB(i) = mean(bAll(m)); sB(i) = sem(bAll(m));
        end
        idx = find(cnt > 0);
        pAg = getfielddef(R,'pAg',nan(n,1));
        pBg = getfielddef(R,'pBg',nan(n,1));

        figure('Position',[0.5 0.5 620 340]); hold on;
        w = 0.38;
        hA = gobjects(1); hB = gobjects(1);
        for q = 1:numel(idx)
            i = idx(q); col = cl_mat(i,:);
            h1 = bar(q-w/2, mA(i), w, 'FaceColor', 1-0.45*(1-col), 'EdgeColor','k');
            h2 = bar(q+w/2, mB(i), w, 'FaceColor', col, 'EdgeColor','k');
            if q==1, hA = h1; hB = h2; end
        end
        errorbar((1:numel(idx))-w/2, mA(idx), sA(idx), 'k','linestyle','none', ...
                 'LineWidth',1.1,'CapSize',5);
        errorbar((1:numel(idx))+w/2, mB(idx), sB(idx), 'k','linestyle','none', ...
                 'LineWidth',1.1,'CapSize',5);

        yr = max([mA(idx)+sA(idx); mB(idx)+sB(idx)]);
        for q = 1:numel(idx)
            i = idx(q);
            put_star(q-w/2, mA(i), sA(i), pAg(i), yr);
            put_star(q+w/2, mB(i), sB(i), pBg(i), yr);
        end

        set(gca,'XTick',1:numel(idx),'XTickLabel',GRP(idx), ...
                'XTickLabelRotation',45,'TickLabelInterpreter','none');
        xlim([0.4 numel(idx)+0.6]); ylim([0 yr*1.22]);
        ylabel('squared separation');
        legend([hA hB], {'on-axis (a)','null-space (b)'}, ...
               'Location','northoutside','Orientation','horizontal','Box','off');
        title(ttl,'Interpreter','none'); box off; hold off;
    end
end

function put_star(x, mu, se, pv, yr)
    s = pStars(pv);
    if isempty(s), return; end
    if isnan(se), se = 0; end
    text(x, mu+se+0.03*yr, s, 'Color','r', 'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', 'FontWeight','bold', 'FontSize', 9);
end

function s = pStars(p)
    if     isnan(p),  s = '';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = '';
    end
end

function v = getfielddef(S, f, d)
    if isfield(S,f), v = S.(f); else, v = d; end
end

function s = sem(x)
    x = x(~isnan(x));
    if numel(x) < 2, s = NaN; else, s = std(x)/sqrt(numel(x)); end
end
