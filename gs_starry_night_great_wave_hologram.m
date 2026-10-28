%% Laser-written Starry Night whose blue-light far field shows The Great Wave
% Same sample and model as gs_starry_night_laser_hologram.m:
%   sky (blue in the picture) -> not processed, transparent;
%   black parts               -> processed to block light completely;
%   yellow parts              -> carry the 0/1 hologram: each "0" site gets one
%                                yellow-type dot that blocks blue light.
% Every site is rasterized with samplesPerSite^2 pixels, and the far field
% of the full raster is computed exactly with the convolution theorem,
%   F = C(f).*DFT(open sites) - D(f).*DFT(dot sites),
% within the first diffraction zone of the grid. Dots must not overlap.
%
% The far-field target is a picture instead of text: the brightness of
% great_wave.png with its dark background set to zero, placed in the +1
% order above the zero order. Its 180-degree-rotated twin appears below the
% zero order. GS reproduces the grey levels on average, with speckle grain on
% top. The target is divided by |D| at the picture position so that the dot's
% own diffraction envelope does not shade the picture. Scalar model with
% uniform coherent illumination of the written square. Uses base MATLAB only.

%% Parameters
folder = fileparts(mfilename('fullpath'));
farFieldPictureFile = fullfile(folder, 'great_wave.png');   % what the far field should show
nearFieldPictureFile = fullfile(folder, 'starry_night_v2.png');
pictureColours = [1 55 175; 0 0 0; 253 239 3]/255;  % sky, black, yellow: classifies the near-field picture
displayColours = [1 55 178; 0 0 0; 255 240 0]/255;  % near-field preview: sky, black, bright-yellow dots
sitesPerSide = 1000;        % laser sites per side: 5 mm at a 5 um pitch
sitePitchUm = 5;            % site grid pitch
dotDiameterUm = 5;          % one laser dot
samplesPerSite = 20;        % raster pixels per site side; each dot is a pixelated disk
wavelengthNm = 405;         % blue light, fully blocked by both dot types
padFactor = 2;              % FFT size / sites per side
carrierPeriodSites = 10/3;  % fringe period; sets the picture's distance from the zero order
pictureSizeFraction = 0.38; % far-field picture side / full far-field window
gsIterations = 200;
rng(1);                     % reproducible random initial phase

%% Site grid and near-field picture classes
if dotDiameterUm > sitePitchUm
    error('Dots overlap (diameter > pitch); the convolution model needs d <= p.');
end
N = padFactor*sitesPerSide;
siteIdx = N/2 - sitesPerSide/2 + (1:sitesPerSide);
[~, Y] = meshgrid((0:N-1) - N/2);  % site units; row index increases with y
written = false(N);                % the illuminated square
written(siteIdx, siteIdx) = true;
nearFieldRgb = imread(nearFieldPictureFile);
siteClass = zeros(N);              % 1 sky, 2 black, 3 yellow; 0 outside the square
siteClass(siteIdx, siteIdx) = fitClassesToGrid(nearestColour(nearFieldRgb, pictureColours), ...
    size(pictureColours, 1), sitesPerSide);
yellow = siteClass == 3;
notBlack = written & siteClass ~= 2;

%% One rasterized site and its Fourier transform at the far-field samples
pixelOffsetsUm = ((0:samplesPerSite-1) + 0.5)/samplesPerSite*sitePitchUm - sitePitchUm/2;
[pixelX, pixelY] = meshgrid(pixelOffsetsUm);
dotPixels = hypot(pixelX, pixelY) <= dotDiameterUm/2;  % the pixelated dot
dotFill = nnz(dotPixels)/samplesPerSite^2;              % dot area / site area
f = ((0:N-1) - N/2)/(N*sitePitchUm);                    % cycles per micrometre
pixelPhase = exp(-2i*pi*f(:)*pixelOffsetsUm)/samplesPerSite;  % per unit site width
siteTransform = pixelPhase*ones(samplesPerSite)*pixelPhase.';  % C(f), rows fy, columns fx
dotTransform = pixelPhase*double(dotPixels)*pixelPhase.';     % D(f)

%% Far-field target
[targetIntensity, placed] = farFieldPicture(farFieldPictureFile, N, pictureSizeFraction);  % centred
carrierShift = round(N/carrierPeriodSites);  % far-field rows from the zero order to the picture centre
target = circshift(targetIntensity, [carrierShift, 0]);  % the picture in the +1 order
pictureBox = circshift(placed, [carrierShift, 0]);
[rows, cols] = find(pictureBox);
if min(rows) <= N/2 + 1 || max(rows) == N
    warning('The picture reaches the zero order or wraps round; adjust carrierPeriodSites or pictureSizeFraction.');
end
gsTarget = sqrt(targetIntensity)./circshift(abs(dotTransform), [-carrierShift, 0]);  % undo |D| at the picture
initialPhase = 2*pi*rand(N);

%% GS phases, then 0/1 carrier encoding into dots
% carrierShift/N cycles per site moves the +1 order by exactly carrierShift rows.
carrier = 2*pi*carrierShift/N*Y;
phaseFull = gerchbergSaxton(gsTarget, double(written), initialPhase, gsIterations);
phaseYellow = gerchbergSaxton(gsTarget, double(yellow), initialPhase, gsIterations);
dotsFull = written & cos(carrier + phaseFull) <= 0;  % "0" sites of the binary hologram
labels = {'① Hologram over the whole square (reference)', ...
    '② Starry Night: hologram ① kept in the yellow areas only', ...
    '③ Starry Night: hologram designed for the yellow areas'};
openSites = {written, notBlack, notBlack};
dotSites = {dotsFull, dotsFull & yellow, yellow & cos(carrier + phaseYellow) <= 0};
hologramSites = {written, yellow, yellow};

%% Far field of the rasterized dots
farFields = cell(1, 3);
for k = 1:3
    farFields{k} = siteTransform.*fft2c(double(openSites{k})) - dotTransform.*fft2c(double(dotSites{k}));
end

%% Metrics
% Powers are relative to the light incident on the written square. With the
% field per unit site area, the power in a set of far-field samples is
% sum(abs(F).^2)/N^2 in units of one site's incident power. Correlations
% compare intensity with the target inside the picture box, raw and after a
% 5 x 5 sample average that removes most speckle grain (about 2 samples).
sites = sitesPerSide^2;
smoothing = ones(5)/25;
bright = target > 0.75;                    % the brightest parts of the picture
darkSky = pictureBox & target < 0.02;      % the picture's dark background
meanPictureFull = mean(abs(farFields{1}(pictureBox)).^2);
[transmitted, coverage, efficiency, brightness, correlation, smoothCorrelation, darkBackground] = ...
    deal(zeros(1, 3));
for k = 1:3
    intensity = abs(farFields{k}).^2;
    smoothIntensity = conv2(intensity, smoothing, 'same');
    transmitted(k) = (nnz(openSites{k}) - dotFill*nnz(dotSites{k}))/sites;
    coverage(k) = dotFill*nnz(dotSites{k})/nnz(hologramSites{k});
    efficiency(k) = sum(intensity(pictureBox))/(N^2*sites);
    brightness(k) = mean(intensity(pictureBox))/meanPictureFull;
    c = corrcoef(intensity(pictureBox), target(pictureBox));
    correlation(k) = c(1, 2);
    c = corrcoef(smoothIntensity(pictureBox), target(pictureBox));
    smoothCorrelation(k) = c(1, 2);
    darkBackground(k) = mean(intensity(darkSky))/mean(intensity(bright));
end
sizeMm = sitesPerSide*sitePitchUm/1000;
fprintf(['\nLaser-written sample: %.2f mm square, %d × %d sites, %g µm dots on a %g µm pitch ' ...
    '(%d × %d pixels per site), %d nm light; far-field picture of %d × %d samples\n' ...
    'Area: sky %.1f%%, black %.1f%%, yellow %.1f%%\n'], sizeMm, sitesPerSide, sitesPerSide, ...
    dotDiameterUm, sitePitchUm, samplesPerSite, samplesPerSite, wavelengthNm, max(rows) - min(rows) + 1, ...
    max(cols) - min(cols) + 1, 100*nnz(siteClass == 1)/sites, 100*nnz(siteClass == 2)/sites, ...
    100*nnz(yellow)/sites);
for k = 1:3
    fprintf(['%s\n    transmitted %5.1f%%  dot coverage %4.1f%%  efficiency %.3f%%  ' ...
        'relative brightness %.3f  correlation %.3f (%.3f after 5 × 5 averaging)\n'], labels{k}, ...
        100*transmitted(k), 100*coverage(k), 100*efficiency(k), brightness(k), correlation(k), ...
        smoothCorrelation(k));
end

%% Figures
siteAxisMm = ((0:sitesPerSide-1) - sitesPerSide/2)*sitePitchUm/1000;
thetaDeg = asind(f*wavelengthNm/1000);
viewMargin = round(0.02*N);
viewRows = max(min(rows) - viewMargin, 1):min(max(rows) + viewMargin, N);
viewCols = max(min(cols) - viewMargin, 1):min(max(cols) + viewMargin, N);
level = linspace(0, 1, 256)';
blueLight = [0.12*level, 0.35*level, level.^0.6];  % black to saturated blue, no white

% Plain 2-by-3 grid: in R2026a, row-spanning tiles let tile titles overlap the layout title.
fig = figure('Name', 'Starry Night sample with The Great Wave in the far field', 'Position', [40 30 1400 900]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
layout.OuterPosition = [0 0.15 1 0.85];  % room for the explanation below
title(layout, sprintf(['Laser-written %g × %g mm sample: Starry Night in the near field, ' ...
    'The Great Wave in the %d nm far field'], sizeMm, sizeMm, wavelengthNm));
subtitle(layout, {['Top: the sample as written (blue = unprocessed and clear, black = opaque, ' ...
    'yellow = sites carrying a yellow dot)'], ['Bottom: far-field intensity in the +1 diffraction ' ...
    'order, where the picture forms (the full far field is in Figure 2)']});
for k = 1:3
    ax = nexttile(layout, k);
    image(ax, siteAxisMm, siteAxisMm, siteColours(openSites{k}(siteIdx, siteIdx), ...
        dotSites{k}(siteIdx, siteIdx), displayColours));
    axis(ax, 'xy'); axis(ax, 'image');
    xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
    title(ax, {labels{k}, sprintf('transmission %.1f%%, dot coverage %.0f%%', ...
        100*transmitted(k), 100*coverage(k))});

    ax = nexttile(layout, 3 + k);
    intensity = abs(farFields{k}).^2;
    imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), intensity(viewRows, viewCols), ...
        [0, 1.3*mean(intensity(bright))]);
    axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, blueLight);
    xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
    title(ax, {sprintf('efficiency %.2f%%, relative brightness %.3f', 100*efficiency(k), brightness(k)), ...
        sprintf('correlation with the target %.3f (%.3f after 5 × 5 averaging)', correlation(k), ...
        smoothCorrelation(k))});
end
addCaption(fig, 0.14, { ...
    sprintf(['How it is made. The laser writes %g µm dots on a %g µm grid. Blue parts of the painting stay ' ...
    'unprocessed and clear, black parts are written opaque, and the yellow parts hold a binary (0/1) ' ...
    'hologram: every site where the hologram is 0 receives one yellow dot, which blocks %d nm light.'], ...
    dotDiameterUm, sitePitchUm, wavelengthNm), ...
    ['How the picture forms. The Gerchberg–Saxton algorithm finds a phase φ whose far-field (Fraunhofer) ' ...
    'pattern is The Great Wave; φ is written into a 0/1 carrier grating, cos(2πy/P + φ) > 0, so the ' ...
    'picture appears in the +1 diffraction order.'], ...
    ['Reading the panels. ② is the original idea: hologram ① covered by the painting. The picture survives ' ...
    'because every part of a Fourier hologram contributes to the whole image, but it turns hazy because ' ...
    'the yellow areas are fragmented. ③ designs the hologram for the yellow areas from the start and ' ...
    'recovers most of the detail.'], ...
    ['Numbers. Transmission is the share of the incident light that passes the sample, and dot coverage the ' ...
    'share of the hologram area under yellow dots. Efficiency is the share of the incident light that lands ' ...
    'in the picture. Correlation compares the far-field intensity with the target; averaging over 5 × 5 ' ...
    'samples removes most of the speckle grain.']});

% Close-up on the most even boundary between yellow dots and clear sky.
zoomSites = 40;
skyCount = conv2(ones(zoomSites, 1), ones(1, zoomSites), double(siteClass(siteIdx, siteIdx) == 1), 'valid');
yellowCount = conv2(ones(zoomSites, 1), ones(1, zoomSites), double(siteClass(siteIdx, siteIdx) == 3), 'valid');
[colGrid, rowGrid] = meshgrid((1:size(skyCount, 2)) - (size(skyCount, 2) + 1)/2, ...
    (1:size(skyCount, 1)) - (size(skyCount, 1) + 1)/2);
score = min(skyCount, yellowCount) - 0.01*hypot(rowGrid, colGrid);  % ties go to the centre
[~, best] = max(score(:));
[zoomRow, zoomCol] = ind2sub(size(score), best);
zoomRows = siteIdx(zoomRow + (0:zoomSites-1));
zoomCols = siteIdx(zoomCol + (0:zoomSites-1));
zoomAxisUm = @(idx) ([idx(1), idx(end)] - (N/2 + 1))*sitePitchUm;

fig2 = figure('Name', 'The Great Wave: target, full far field and sample', 'Position', [80 30 1100 900]);
layout = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
layout.OuterPosition = [0 0.12 1 0.88];  % room for the explanation below
ax = nexttile(layout);
imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), target(viewRows, viewCols), [0 1]);
axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, blueLight);
xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
title(ax, {'Far-field target', 'brightness of The Great Wave, dark background set to zero'});
ax = nexttile(layout);
intensity = abs(farFields{3}).^2;
imagesc(ax, thetaDeg, thetaDeg, intensity, [0, 1.3*mean(intensity(bright))]);
axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, blueLight);
xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
title(ax, {'③ Full far field', 'the picture, the zero order and the twin'});
labelX = thetaDeg(max(cols)) + 0.1;
pictureCentreDeg = thetaDeg(round((min(rows) + max(rows))/2));
text(ax, labelX, pictureCentreDeg, {'+1 order:', 'the picture'}, 'Color', 'w', 'FontSize', 9);
text(ax, labelX, 0, {'zero order', '(saturated)'}, 'Color', 'w', 'FontSize', 9);
text(ax, labelX, -pictureCentreDeg, {'−1 order:', 'rotated twin'}, 'Color', 'w', 'FontSize', 9);
ax = nexttile(layout);
pictureSizeMm = size(nearFieldRgb, [2 1])/max(size(nearFieldRgb, [1 2]))*sizeMm;
image(ax, pictureSizeMm(1)*[-0.5 0.5], pictureSizeMm(2)*[-0.5 0.5], flipud(nearFieldRgb));
axis(ax, 'xy'); axis(ax, 'image');
xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
title(ax, {'Near-field design: Starry Night', 'blue = unprocessed, black = opaque, yellow = hologram dots'});
ax = nexttile(layout);
image(ax, zoomAxisUm(zoomCols), zoomAxisUm(zoomRows), dotImage(openSites{3}(zoomRows, zoomCols), ...
    dotSites{3}(zoomRows, zoomCols), dotPixels, displayColours));
axis(ax, 'xy'); axis(ax, 'image');
xlabel(ax, 'x (µm)'); ylabel(ax, 'y (µm)');
title(ax, {sprintf('③ Close-up of %d × %d sites', zoomSites, zoomSites), ...
    sprintf('each %g µm dot is drawn with %d × %d pixels, as in the simulation', dotDiameterUm, ...
    samplesPerSite, samplesPerSite)});
addCaption(fig2, 0.11, { ...
    sprintf(['Why there is a twin. A 0/1 transmission is real, so its far-field intensity is point-symmetric ' ...
    'and the −1 order repeats the picture rotated by 180°. The zero order is light the hologram does not ' ...
    'diffract; it spreads into the diffraction pattern of the clear parts of the sample. The dark sky of the ' ...
    'picture keeps a faint speckle background (%.0f%% of the bright level in ③), so the picture''s square ' ...
    'shows faintly.'], 100*darkBackground(3)), ...
    sprintf(['How the far field is computed. Every site is drawn with %d × %d pixels. The far field of that ' ...
    'whole pixel image follows exactly from the convolution theorem, F = C(f)·DFT(clear sites) − ' ...
    'D(f)·DFT(dot sites), where C and D are the Fourier transforms of one clear site and one pixelated dot.'], ...
    samplesPerSite, samplesPerSite)});

%% Local functions
function [phase, farField, errorHistory] = gerchbergSaxton(targetAmplitude, aperture, initialPhase, iterations)
% Alternate between the far-field amplitude target and a phase-only
% hologram whose amplitude equals the aperture (zero outside it).
% Arrays stay in FFT order inside the loop; the outputs are centred.
target = ifftshift(targetAmplitude);
aperture = ifftshift(aperture);
farConstrained = target .* exp(1i*ifftshift(initialPhase));
errorHistory = zeros(iterations, 1);
for k = 1:iterations
    hologram = aperture .* unitPhasor(ifft2(farConstrained));
    farField = fft2(hologram);
    errorHistory(k) = amplitudeError(farField, target);
    farConstrained = target .* unitPhasor(farField);
end
phase = fftshift(angle(hologram));
farField = fftshift(farField);
end

function z = unitPhasor(z)
% exp(1i*angle(z)) without the trigonometric calls; zero stays zero.
z = z ./ max(abs(z), realmin);
end

function e = amplitudeError(farField, targetAmplitude)
% RMS distance between unit-energy far-field and target amplitudes.
amplitude = abs(farField);
e = norm(amplitude(:)/norm(amplitude(:)) - targetAmplitude(:)/norm(targetAmplitude(:)));
end

function F = fft2c(field)
F = fftshift(fft2(ifftshift(field)));
end

function [intensity, placed] = farFieldPicture(file, N, sizeFraction)
% Brightness of a single-hue picture (its brightest channel), dark
% background set to zero and maximum to one, area-averaged onto a centred
% square of sizeFraction*N far-field samples with its aspect ratio kept.
% placed marks that square. Rows are flipped so the picture reads upright
% with axis xy.
value = max(double(imread(file))/255, [], 3);
value = (value - min(value(:)))/(max(value(:)) - min(value(:)));
side = round(sizeFraction*N);
scale = max(size(value))/side;            % picture pixels per far-field sample
k = max(2*round((scale - 1)/2) + 1, 1);   % odd box close to one sample
value = conv2(flipud(value), ones(k)/k^2, 'same');
offsets = ((1:side) - (side + 1)/2)*scale;
[xq, yq] = meshgrid((size(value, 2) + 1)/2 + offsets, (size(value, 1) + 1)/2 + offsets);
square = N/2 - floor(side/2) + (1:side);
intensity = zeros(N);
intensity(square, square) = interp2(value, xq, yq, 'linear', 0);
placed = false(N);
placed(square, square) = true;
end

function classes = nearestColour(rgb, colours)
% Index of the closest palette colour for every picture pixel.
pixels = double(reshape(rgb, [], 3))/255;
distance = zeros(size(pixels, 1), size(colours, 1));
for c = 1:size(colours, 1)
    distance(:, c) = sum((pixels - colours(c, :)).^2, 2);
end
[~, classes] = min(distance, [], 2);
classes = reshape(classes, size(rgb, 1), size(rgb, 2));
end

function siteClass = fitClassesToGrid(classes, classCount, sitesPerSide)
% Majority class of the picture area each site covers. The picture keeps
% its aspect ratio (centred); sites outside it get class 1 (unprocessed).
% Rows are flipped so the picture reads upright with axis xy.
classes = flipud(classes);
scale = max(size(classes))/sitesPerSide;  % picture pixels per site
k = max(2*round((scale - 1)/2) + 1, 1);   % odd box close to one site
offsets = ((1:sitesPerSide) - (sitesPerSide + 1)/2)*scale;
[xq, yq] = meshgrid((size(classes, 2) + 1)/2 + offsets, (size(classes, 1) + 1)/2 + offsets);
fractions = zeros(sitesPerSide, sitesPerSide, classCount);
for c = 1:classCount
    fractions(:, :, c) = interp2(conv2(double(classes == c), ones(k)/k^2, 'same'), xq, yq, 'linear', 0);
end
[~, siteClass] = max(fractions, [], 3);
end

function rgb = siteColours(open, dots, colours)
% One colour per site: clear sites show the sky colour, black sites are
% black, and sites carrying a dot are yellow.
rgb = double(open).*reshape(colours(1, :), 1, 1, 3) + ...
    double(dots).*reshape(colours(3, :) - colours(1, :), 1, 1, 3);
end

function rgb = dotImage(open, dots, dotPixels, colours)
% Sites drawn with the same pixelated dot that the far field uses.
samples = size(dotPixels, 1);
rgb = kron(double(open), ones(samples)).*reshape(colours(1, :), 1, 1, 3) + ...
    kron(double(dots), double(dotPixels)).*reshape(colours(3, :) - colours(1, :), 1, 1, 3);
end

function addCaption(fig, height, paragraphs)
% Explanatory paragraphs along the bottom of a figure.
annotation(fig, 'textbox', [0.015 0.005 0.97 height], 'String', paragraphs, 'EdgeColor', 'none', ...
    'FontSize', 9, 'Interpreter', 'none', 'VerticalAlignment', 'top');
end
