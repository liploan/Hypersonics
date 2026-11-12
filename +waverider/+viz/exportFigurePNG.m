function exportFigurePNG(fig, pngPath, dpi)
%EXPORTFIGUREPNG Export a MATLAB figure to a true RGBA transparent-background PNG.
%   EXPORTFIGUREPNG(fig, pngPath, dpi) finalizes high-contrast publication
%   typography on fig and writes a 32-bit RGBA PNG at resolution dpi (default
%   200) whose outer figure background is genuinely transparent (Alpha = 0).
%
%   Method:
%     1. Renders identical pixel grids against pure white (fig.Color = [1 1 1])
%        and pure black (fig.Color = [0 0 0]) via MATLAB's OpenGL/painters
%        rasterizer (`print -dpng`), recovering exact per-pixel anti-aliased
%        opacity:
%            alpha = 1 - mean(I_white - I_black, 3)
%        and premultiplied foreground color C_pre = I_black.
%     2. Adds a tight 2.5-pixel white contrast halo strictly behind dark
%        outer-margin text/tick strokes (where neighboring pixels are
%        transparent). On light backgrounds the white halo blends invisibly;
%        on dark backgrounds (e.g., GitHub Dark Mode #0d1117 or dark slides)
%        every title, axis label, tick number, and colorbar label remains
%        crisp and legible while preserving a transparent background.

    if nargin < 3 || isempty(dpi)
        dpi = 200;
    end

    waverider.viz.stylePublicationFigure(fig, 'finalize');
    drawnow;

    origColor = get(fig, 'Color');
    origInvert = get(fig, 'InvertHardcopy');
    origPPM = get(fig, 'PaperPositionMode');

    tmpWhite = [tempname, '_w.png'];
    tmpBlack = [tempname, '_b.png'];
    cleanupObj = onCleanup(@() deleteTempFiles(tmpWhite, tmpBlack));

    set(fig, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
    resFlag = sprintf('-r%d', round(dpi));

    % Pass 1: Render over pure white canvas
    set(fig, 'Color', [1.0, 1.0, 1.0]);
    drawnow;
    print(fig, tmpWhite, '-dpng', resFlag);

    % Pass 2: Render over pure black canvas (identical pixel dimensions)
    set(fig, 'Color', [0.0, 0.0, 0.0]);
    drawnow;
    print(fig, tmpBlack, '-dpng', resFlag);

    % Restore figure properties
    set(fig, 'Color', origColor, 'InvertHardcopy', origInvert, 'PaperPositionMode', origPPM);

    Iw = double(imread(tmpWhite)) ./ 255.0;
    Ib = double(imread(tmpBlack)) ./ 255.0;

    % Exact Porter-Duff alpha recovery
    alpha0 = 1.0 - mean(Iw - Ib, 3);
    alpha0 = min(max(alpha0, 0.0), 1.0);
    % Snap numerical quantization noise (< 1/255) to exact 0 or 1
    alpha0(alpha0 < 0.008) = 0.0;
    alpha0(alpha0 > 0.992) = 1.0;

    % Identify dark ink strokes located in the transparent outer margin
    % (i.e., opaque dark pixels within ~10 px of transparent canvas pixels)
    transMask = double(alpha0 < 0.05);
    kNear = ones(21, 21) ./ 441.0;
    nearTransparent = conv2(transMask, kNear, 'same') > 0.02;

    % Unpremultiplied luminance of foreground pixels
    fgR = Ib(:, :, 1) ./ max(alpha0, 1e-4);
    fgG = Ib(:, :, 2) ./ max(alpha0, 1e-4);
    fgB = Ib(:, :, 3) ./ max(alpha0, 1e-4);
    fgLum = 0.2126 .* fgR + 0.7152 .* fgG + 0.0722 .* fgB;

    outerDarkStroke = double((alpha0 > 0.15) & (fgLum < 0.55) & nearTransparent);

    % Build smooth 3-px radial contrast halo behind outer dark text/ticks
    rHalo = max(2, round(dpi / 60)); % ~3 px at 200 DPI
    [xx, yy] = meshgrid(-rHalo:rHalo, -rHalo:rHalo);
    rr = sqrt(xx.^2 + yy.^2);
    kHalo = max(0.0, 1.0 - rr ./ (rHalo + 0.75));
    kHalo = kHalo ./ max(sum(kHalo(:)), 1e-8);

    % Multi-pass soft dilation for smooth anti-aliased white backing
    haloField = conv2(outerDarkStroke, kHalo, 'same');
    haloField = conv2(min(haloField .* 3.2, 1.0), kHalo, 'same');
    alphaHalo = min(max(haloField .* 1.35, 0.0), 0.92);
    alphaHalo(alphaHalo < 0.02) = 0.0;

    % Composite foreground (C_pre = Ib, alpha0) OVER white halo (RGB=1, alphaHalo)
    alphaOut = alpha0 + (1.0 - alpha0) .* alphaHalo;
    bgContrib = (1.0 - alpha0) .* alphaHalo;

    R_out = (Ib(:, :, 1) + bgContrib) ./ max(alphaOut, 1e-6);
    G_out = (Ib(:, :, 2) + bgContrib) ./ max(alphaOut, 1e-6);
    B_out = (Ib(:, :, 3) + bgContrib) ./ max(alphaOut, 1e-6);
    RGB_out = min(max(cat(3, R_out, G_out, B_out), 0.0), 1.0);

    % Zero RGB where completely transparent for clean PNG compression
    zeroMask = (alphaOut == 0);
    for c = 1:3
        ch = RGB_out(:, :, c);
        ch(zeroMask) = 0.0;
        RGB_out(:, :, c) = ch;
    end

    % Auto-crop outer transparent margins with uniform padding
    [rows, cols] = find(alphaOut > 0.01);
    if ~isempty(rows)
        pad = round(dpi * 0.08); % ~16 px at 200 DPI
        r1 = max(1, min(rows) - pad);
        r2 = min(size(alphaOut, 1), max(rows) + pad);
        c1 = max(1, min(cols) - pad);
        c2 = min(size(alphaOut, 2), max(cols) + pad);
        RGB_out = RGB_out(r1:r2, c1:c2, :);
        alphaOut = alphaOut(r1:r2, c1:c2);
    end

    outDir = fileparts(pngPath);
    if ~isempty(outDir) && ~isfolder(outDir)
        mkdir(outDir);
    end
    imwrite(uint8(round(RGB_out .* 255.0)), pngPath, 'png', ...
        'Alpha', uint8(round(alphaOut .* 255.0)));
end

function deleteTempFiles(f1, f2)
    if isfile(f1), delete(f1); end
    if isfile(f2), delete(f2); end
end
