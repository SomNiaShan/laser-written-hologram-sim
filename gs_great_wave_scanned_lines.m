%% The Great Wave hologram written by continuous scanning instead of single dots
% Runs gs_starry_night_great_wave_hologram.m (same holograms, same far-field
% target), then rewrites every hologram as the laser would draw it when it
% scans along x and switches on and off: each run of neighbouring dot sites
% in a row becomes one continuous line, a stadium of width lineWidthUm with
% round ends (the line swept by a round spot), with no gaps between sites.
% The far field stays exact: every dot site is one of four tile types
% (single dot, left end, middle, right end of a line), so
%   F = C(f).*DFT(open sites) - sum over tiles T(f).*DFT(sites of that type).
% Needs lineWidthUm <= sitePitchUm so that neighbouring rows do not overlap.

gs_starry_night_great_wave_hologram;   % dot model: holograms, target, metrics
dotFarFields = farFields;
dotStats = [transmitted; efficiency; brightness; correlation; smoothCorrelation];
close all

lineWidthUm = dotDiameterUm;          % width of the scanned line (spot size)
if lineWidthUm > sitePitchUm
    error('Lines from neighbouring rows overlap (width > pitch).');
end

%% Tiles: the part of one site covered by each kind of line piece
r = lineWidthUm/2;
inBand = abs(pixelY) <= r;             % the straight part of a line along x
disk = hypot(pixelX, pixelY) <= r;     % a round end centred on the site
tileNames = {'single', 'left end', 'middle', 'right end'};
tiles = {disk, disk | (inBand & pixelX >= 0), inBand, disk | (inBand & pixelX <= 0)};
tileTransforms = cellfun(@(t) pixelPhase*double(t)*pixelPhase.', tiles, 'UniformOutput', false);
tileFill = cellfun(@(t) nnz(t)/samplesPerSite^2, tiles);

%% Far fields of the scanned versions of ①, ②, ③
lineFarFields = cell(1, 3);
lineMasks = cell(1, 3);
[runCount, singleRuns, medianRunSites, lineTransmitted] = deal(zeros(1, 3));
for k = 1:3
    m = splitIntoLinePieces(dotSites{k});
    lineMasks{k} = m;
    F = siteTransform.*fft2c(double(openSites{k}));
    for t = 1:4
        F = F - tileTransforms{t}.*fft2c(double(m{t}));
    end
    lineFarFields{k} = F;
    lineTransmitted(k) = (nnz(openSites{k}) - sum(tileFill.*cellfun(@nnz, m)))/sites;
    runs = runLengths(dotSites{k});
    runCount(k) = numel(runs);
    singleRuns(k) = nnz(runs == 1);
    medianRunSites(k) = median(runs);
end

%% Metrics, same definitions as the dot model
[lineEff, lineBright, lineCorr, lineSmooth, lineDark] = deal(zeros(1, 3));
for k = 1:3
    intensity = abs(lineFarFields{k}).^2;
    smoothIntensity = conv2(intensity, smoothing, 'same');
    lineEff(k) = sum(intensity(pictureBox))/(N^2*sites);
    lineBright(k) = mean(intensity(pictureBox))/meanPictureFull;
    c = corrcoef(intensity(pictureBox), target(pictureBox));      lineCorr(k) = c(1, 2);
    c = corrcoef(smoothIntensity(pictureBox), target(pictureBox)); lineSmooth(k) = c(1, 2);
    lineDark(k) = mean(intensity(darkSky))/mean(intensity(bright));
end
fprintf('\nScanned lines: width %g µm, along x, round ends\n', lineWidthUm);
fprintf('%-58s %8s %8s %8s %8s %8s\n', '', 'transm', 'eff', 'bright', 'corr', 'corr5x5');
for k = 1:3
    fprintf('%-58s %7.1f%% %7.3f%% %8.3f %8.3f %8.3f\n', [labels{k} ' — dots'], 100*dotStats(1, k), ...
        100*dotStats(2, k), dotStats(3, k), dotStats(4, k), dotStats(5, k));
    fprintf('%-58s %7.1f%% %7.3f%% %8.3f %8.3f %8.3f\n', [labels{k} ' — lines'], 100*lineTransmitted(k), ...
        100*lineEff(k), lineBright(k), lineCorr(k), lineSmooth(k));
    fprintf('    %d lines, %d of them single sites, median length %g sites (%g µm)\n', ...
        runCount(k), singleRuns(k), medianRunSites(k), medianRunSites(k)*sitePitchUm);
end

%% Figure: close-ups and far fields, dots versus lines for ③
fig3 = figure('Name', 'Dots versus scanned lines', 'Position', [60 30 1400 900]);
layout = tiledlayout(fig3, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('③ written as %g µm dots or as continuous %g µm lines scanned along x', ...
    dotDiameterUm, lineWidthUm));
ax = nexttile(layout);
image(ax, zoomAxisUm(zoomCols), zoomAxisUm(zoomRows), dotImage(openSites{3}(zoomRows, zoomCols), ...
    dotSites{3}(zoomRows, zoomCols), dotPixels, displayColours));
axis(ax, 'xy'); axis(ax, 'image'); xlabel(ax, 'x (µm)'); ylabel(ax, 'y (µm)');
title(ax, {'Dots (as simulated so far)', sprintf('%d × %d sites', zoomSites, zoomSites)});
ax = nexttile(layout);
image(ax, zoomAxisUm(zoomCols), zoomAxisUm(zoomRows), lineImage(openSites{3}(zoomRows, zoomCols), ...
    cellfun(@(m) m(zoomRows, zoomCols), lineMasks{3}, 'UniformOutput', false), tiles, displayColours));
axis(ax, 'xy'); axis(ax, 'image'); xlabel(ax, 'x (µm)'); ylabel(ax, 'y (µm)');
title(ax, {'Scanned lines', 'same sites, joined along x with round ends'});
ax = nexttile(layout);
runs = runLengths(dotSites{3});
histogram(ax, runs, 0.5:1:min(max(runs), 40) + 0.5);
xlabel(ax, 'line length (sites of 5 µm)'); ylabel(ax, 'number of lines');
title(ax, {sprintf('%d lines in ③', numel(runs)), sprintf('median %g sites, %d single sites', ...
    median(runs), nnz(runs == 1))});
farTitles = {'③ dots', '③ scanned lines', '① scanned lines (reference)'};
farData = {dotFarFields{3}, lineFarFields{3}, lineFarFields{1}};
farStats = [dotStats(2, 3) dotStats(5, 3); lineEff(3) lineSmooth(3); lineEff(1) lineSmooth(1)];
for k = 1:3
    ax = nexttile(layout);
    intensity = abs(farData{k}).^2;
    imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), intensity(viewRows, viewCols), ...
        [0, 1.3*mean(intensity(bright))]);
    axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, blueLight);
    xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
    title(ax, {farTitles{k}, sprintf('efficiency %.2f%%, correlation %.3f after 5 × 5 averaging', ...
        100*farStats(k, 1), farStats(k, 2))});
end

%% Local functions
function m = splitIntoLinePieces(d)
% Sites of each tile type for runs of true along x (columns).
leftOn = circshift(d, [0 1]);    % the site to the left also has a dot
rightOn = circshift(d, [0 -1]);  % the site to the right also has a dot
m = {d & ~leftOn & ~rightOn, d & ~leftOn & rightOn, d & leftOn & rightOn, d & leftOn & ~rightOn};
end

function runs = runLengths(d)
% Lengths, in sites, of every run of true along x.
padded = [false(size(d, 1), 1), d, false(size(d, 1), 1)]';
edges = diff(padded(:));
runs = find(edges == -1) - find(edges == 1);
end

function F = fft2c(field)
F = fftshift(fft2(ifftshift(field)));
end

function rgb = dotImage(open, dots, dotPixels, colours)
samples = size(dotPixels, 1);
rgb = kron(double(open), ones(samples)).*reshape(colours(1, :), 1, 1, 3) + ...
    kron(double(dots), double(dotPixels)).*reshape(colours(3, :) - colours(1, :), 1, 1, 3);
end

function rgb = lineImage(open, masks, tiles, colours)
samples = size(tiles{1}, 1);
covered = zeros(size(open)*samples);
for t = 1:numel(tiles)
    covered = covered + kron(double(masks{t}), double(tiles{t}));
end
rgb = kron(double(open), ones(samples)).*reshape(colours(1, :), 1, 1, 3) + ...
    covered.*reshape(colours(3, :) - colours(1, :), 1, 1, 3);
end
