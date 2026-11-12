function stylePublicationFigure(fig, mode)
%STYLEPUBLICATIONFIGURE Enforce high-contrast publication typography & axes styling.
%   STYLEPUBLICATIONFIGURE(fig, 'init') configures a newly created figure
%   so MATLAB R2025a+ dark desktop themes never bleed light-gray text onto
%   exported figures.
%   STYLEPUBLICATIONFIGURE(fig, 'finalize') walks all TiledChartLayouts,
%   Axes, ColorBars, Legends, and ConstantLines in fig and enforces crisp
%   WCAG AAA contrast, legible font sizes, and clean gridlines.

    if nargin < 2 || isempty(mode)
        mode = 'finalize';
    end

    ink = [0.06, 0.09, 0.16];       % #0f172a deep slate-black (17.8:1 on white)
    panelBg2D = [1.00, 1.00, 1.00]; % Crisp white 2D plot canvas
    panelBg3D = [0.965, 0.975, 0.99]; % Subtle cool-white 3D box walls
    gridCol = [0.38, 0.44, 0.52];
    borderCol = [0.64, 0.70, 0.78];

    % Always force light figure theme in MATLAB R2025a+ so defaults are dark-on-light
    if isprop(fig, 'Theme')
        try
            theme(fig, 'light');
        catch
            try %#ok<TRYNC>
                fig.Theme = 'light';
            end
        end
    end
    set(fig, 'Color', [1.0, 1.0, 1.0], 'InvertHardcopy', 'off');

    if strcmpi(mode, 'init')
        return;
    end

    % 1. Style TiledChartLayout super-titles & spacing
    tls = findall(fig, 'Type', 'tiledlayout');
    for k = 1:numel(tls)
        tl = tls(k);
        tl.TileSpacing = 'loose';
        tl.Padding = 'compact';
        if ~isempty(tl.Title) && strlength(string(tl.Title.String)) > 0
            tl.Title.Color = ink;
            tl.Title.FontSize = max(tl.Title.FontSize, 14);
            tl.Title.FontWeight = 'bold';
        end
        if ~isempty(tl.Subtitle) && strlength(string(tl.Subtitle.String)) > 0
            tl.Subtitle.Color = [0.20, 0.25, 0.34];
            tl.Subtitle.FontSize = max(tl.Subtitle.FontSize, 11);
        end
    end

    % 2. Style all Cartesian Axes
    axList = findall(fig, 'Type', 'axes');
    for k = 1:numel(axList)
        ax = axList(k);
        [azView, elView] = view(ax);
        is3D = ~(abs(azView) < 1e-6 && abs(elView - 90) < 1e-6);

        if is3D
            ax.Color = panelBg3D;
        else
            ax.Color = panelBg2D;
        end

        ax.FontSize = max(ax.FontSize, 10);
        ax.LineWidth = max(ax.LineWidth, 1.05);
        ax.GridColor = gridCol;
        ax.GridAlpha = 0.28;
        ax.MinorGridColor = gridCol;
        ax.MinorGridAlpha = 0.16;
        ax.XColor = ink;
        ax.ZColor = ink;

        if numel(ax.YAxis) == 2
            % Dual yyaxis: darken each side if needed for >= 5.5:1 contrast on white
            for s = 1:2
                cSide = darkenForWhiteContrast(ax.YAxis(s).Color);
                ax.YAxis(s).Color = cSide;
                ax.YAxis(s).FontSize = max(ax.YAxis(s).FontSize, 10);
                if ~isempty(ax.YAxis(s).Label)
                    ax.YAxis(s).Label.Color = cSide;
                    ax.YAxis(s).Label.FontSize = max(ax.YAxis(s).Label.FontSize, 10.5);
                    ax.YAxis(s).Label.FontWeight = 'bold';
                end
            end
        else
            ax.YColor = ink;
            if ~isempty(ax.YLabel)
                ax.YLabel.Color = ink;
                ax.YLabel.FontSize = max(ax.YLabel.FontSize, 10.5);
                ax.YLabel.FontWeight = 'bold';
            end
        end

        if ~isempty(ax.XLabel)
            ax.XLabel.Color = ink;
            ax.XLabel.FontSize = max(ax.XLabel.FontSize, 10.5);
            ax.XLabel.FontWeight = 'bold';
        end
        if ~isempty(ax.ZLabel)
            ax.ZLabel.Color = ink;
            ax.ZLabel.FontSize = max(ax.ZLabel.FontSize, 10.5);
            ax.ZLabel.FontWeight = 'bold';
        end
        if ~isempty(ax.Title)
            ax.Title.Color = ink;
            ax.Title.FontSize = max(ax.Title.FontSize, 11.5);
            ax.Title.FontWeight = 'bold';
        end

        % Style xline / yline ConstantLines depending on whether ax has heatmap/contourf
        hasFilledField = ~isempty(findall(ax, 'Type', 'image')) || ...
                         ~isempty(findall(ax, 'Type', 'contour', 'Fill', 'on'));
        cLines = findall(ax, 'Type', 'constantline');
        for m = 1:numel(cLines)
            cl = cLines(m);
            cl.FontSize = max(cl.FontSize, 9.0);
            cl.FontWeight = 'bold';
            if hasFilledField
                % Keep high-visibility light or dark color chosen by caller
                if mean(cl.Color) > 0.75
                    cl.Color = [1.0, 1.0, 1.0];
                else
                    cl.Color = ink;
                end
            else
                cl.Color = darkenForWhiteContrast(cl.Color);
            end
        end
    end

    % 3. Style all ColorBars
    cbs = findall(fig, 'Type', 'colorbar');
    for k = 1:numel(cbs)
        cb = cbs(k);
        cb.Color = ink;
        cb.FontSize = max(cb.FontSize, 9.5);
        cb.LineWidth = 1.0;
        if ~isempty(cb.Label)
            cb.Label.Color = ink;
            cb.Label.FontSize = max(cb.Label.FontSize, 10.5);
            cb.Label.FontWeight = 'bold';
        end
    end

    % 4. Style all Legends
    lgs = findall(fig, 'Type', 'legend');
    for k = 1:numel(lgs)
        lg = lgs(k);
        lg.TextColor = ink;
        lg.Color = [1.0, 1.0, 1.0];
        lg.EdgeColor = borderCol;
        lg.LineWidth = 0.9;
        lg.FontSize = max(lg.FontSize, 8.5);
    end
end

function cOut = darkenForWhiteContrast(cIn)
%DARKENFORWHITECONTRAST Ensure RGB color has strong contrast against white.
    if ~isnumeric(cIn) || numel(cIn) ~= 3
        cOut = [0.08, 0.12, 0.20];
        return;
    end
    lum = 0.2126 * cIn(1) + 0.7152 * cIn(2) + 0.0722 * cIn(3);
    if lum > 0.42
        scale = 0.36 / max(lum, 1e-6);
        cOut = max(min(cIn .* scale, 1.0), 0.0);
    elseif lum < 0.03
        cOut = [0.12, 0.16, 0.24];
    else
        cOut = cIn;
    end
end
