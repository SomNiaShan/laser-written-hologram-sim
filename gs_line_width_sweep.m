%% Line-width sweep for the scanned-line model (designs ① and ③)
% Runs gs_starry_night_great_wave_hologram.m, then redraws holograms ① and ③
% as continuous lines along x (round ends, as in gs_great_wave_scanned_lines.m)
% for each width in `widths`, and prints transmission, efficiency and the
% correlation with the target. The figure shows the +1-order far field of ③
% for every width. Widths must not exceed the 5 um site pitch.
gs_starry_night_great_wave_hologram;
close all
widths = [3 3.5 4 4.5 5];
designs = [1 3];
res = zeros(numel(widths), 4, numel(designs));   % transm, eff, corr, corr5x5
keepFar = cell(1, numel(widths));
for w = 1:numel(widths)
    r = widths(w)/2;
    inBand = abs(pixelY) <= r;
    disk = hypot(pixelX, pixelY) <= r;
    tiles = {disk, disk | (inBand & pixelX >= 0), inBand, disk | (inBand & pixelX <= 0)};
    T = cellfun(@(t) pixelPhase*double(t)*pixelPhase.', tiles, 'UniformOutput', false);
    fill = cellfun(@(t) nnz(t)/samplesPerSite^2, tiles);
    for j = 1:numel(designs)
        k = designs(j);
        d = dotSites{k};
        L = circshift(d, [0 1]); R = circshift(d, [0 -1]);
        m = {d & ~L & ~R, d & ~L & R, d & L & R, d & L & ~R};
        F = siteTransform.*fftc(double(openSites{k}));
        for t = 1:4, F = F - T{t}.*fftc(double(m{t})); end
        I = abs(F).^2;
        c1 = corrcoef(I(pictureBox), target(pictureBox));
        S = conv2(I, smoothing, 'same');
        c2 = corrcoef(S(pictureBox), target(pictureBox));
        res(w, :, j) = [(nnz(openSites{k}) - sum(fill.*cellfun(@nnz, m)))/sites, ...
            sum(I(pictureBox))/(N^2*sites), c1(1, 2), c2(1, 2)];
        if k == 3, keepFar{w} = I; end
    end
end
for j = 1:numel(designs)
    fprintf('\nDesign %d (dots: eff %.3f%%, corr5x5 %.3f)\n', designs(j), 100*efficiency(designs(j)), smoothCorrelation(designs(j)));
    fprintf('  width  transm    eff     corr  corr5x5\n');
    for w = 1:numel(widths)
        fprintf('  %4.1f  %5.1f%%  %6.3f%%  %6.3f  %6.3f\n', widths(w), 100*res(w,1,j), 100*res(w,2,j), res(w,3,j), res(w,4,j));
    end
end
fig = figure('Position', [40 40 1500 520]);
lay = tiledlayout(fig, 1, numel(widths), 'TileSpacing', 'compact', 'Padding', 'compact');
title(lay, '③ scanned lines: far field (+1 order) for different line widths');
for w = 1:numel(widths)
    ax = nexttile(lay); I = keepFar{w};
    imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), I(viewRows, viewCols), [0, 1.3*mean(I(bright))]);
    axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, blueLight);
    xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
    title(ax, {sprintf('width %.1f µm', widths(w)), sprintf('eff %.2f%%, corr %.3f', 100*res(w,2,2), res(w,4,2))});
end
function F = fftc(x)
F = fftshift(fft2(ifftshift(x)));
end
