%% Laser-written Starry Night hologram: rasterized round dots, blue-light far field
% The sample is written by laser, dot by dot, on a square grid of sites:
%   sky (blue in the picture) -> not processed, transparent;
%   black parts               -> processed to block light completely;
%   yellow parts              -> carry the 0/1 hologram: each "0" site gets one
%                                yellow-type dot, a disk that blocks blue light,
%                                and "1" sites stay clear.
% Light also passes the clear corners around every dot: a 5 um dot on a 5 um
% grid blocks only about pi/4 of its site. The black parts are taken as
% fully opaque, i.e. written densely enough that neighbouring dots overlap.
%
% Every site is rasterized with samplesPerSite^2 pixels, so each dot is a
% pixelated disk. The whole sample is the site map convolved with that
% one-site raster, so by the convolution theorem the far field of the full
% raster is exactly
%   F = C(f).*DFT(open sites) - D(f).*DFT(dot sites),
% with C and D the Fourier transforms of the rasterized open site and dot.
% Inside the first diffraction zone of the grid, |sin(theta)| < lambda/(2p),
% this equals an FFT of the full raster (20000 x 20000 pixels by default)
% without forming it. Dots must not overlap (d <= p); overlapping dots would
% need their union rasterized explicitly. The site maps are zero-padded by
% padFactor, sampling the far-field intensity at its Nyquist rate.
%
% The hologram is designed as in gs_starry_night_binary_hologram.m: a GS
% phase for the text, written into a 0/1 carrier grating. The +1 order is
% weighted by D(f), which falls off away from the axis, so the GS target is
% divided by |D| at the text position to keep the outer letters as bright as
% the middle ones. Scalar model with uniform coherent illumination of the
% written square. Uses base MATLAB only.

%% Parameters
targetText = 'Vincent van Gogh';  % a cell array of strings gives several lines
pictureFile = fullfile(fileparts(mfilename('fullpath')), 'starry_night_v2.png');
pictureColours = [1 55 175; 0 0 0; 253 239 3]/255;  % sky, black, yellow: classifies the picture
displayColours = [1 55 178; 0 0 0; 255 240 0]/255;  % near-field preview: sky, black, bright-yellow dots
sitesPerSide = 1000;        % laser sites per side: 5 mm at a 5 um pitch
sitePitchUm = 5;            % site grid pitch
dotDiameterUm = 5;          % one laser dot
samplesPerSite = 20;        % raster pixels per site side; each dot is a pixelated disk
wavelengthNm = 405;         % blue light, fully blocked by both dot types
padFactor = 2;              % FFT size / sites per side
carrierPeriodSites = 4;     % binary fringe period; sets the text's distance from the zero order
textWidthFraction = 0.8;    % text width / full far-field window
gsIterations = 200;
rng(1);                     % reproducible random initial phase

%% Site grid and picture classes
if dotDiameterUm > sitePitchUm
    error('Dots overlap (diameter > pitch); the convolution model needs d <= p.');
end
N = padFactor*sitesPerSide;
siteIdx = N/2 - sitesPerSide/2 + (1:sitesPerSide);
[~, Y] = meshgrid((0:N-1) - N/2);  % site units; row index increases with y
written = false(N);                % the illuminated square
written(siteIdx, siteIdx) = true;
pictureRgb = imread(pictureFile);
siteClass = zeros(N);              % 1 sky, 2 black, 3 yellow; 0 outside the square
siteClass(siteIdx, siteIdx) = fitClassesToGrid(nearestColour(pictureRgb, pictureColours), ...
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

%% Target
textAmplitude = renderText(targetText, N, textWidthFraction);  % centred
carrierShift = round(N/carrierPeriodSites);  % far-field rows from the zero order to the text
letters = circshift(textAmplitude > 0, [carrierShift, 0]);     % text in the +1 order
[rows, cols] = find(letters);
if min(rows) <= N/2 + 1
    warning('The text reaches the zero order and its twin; reduce carrierPeriodSites or textWidthFraction.');
end
gsTarget = textAmplitude./circshift(abs(dotTransform), [-carrierShift, 0]);  % undo |D| at the text
initialPhase = 2*pi*rand(N);

%% GS phases, then 0/1 carrier encoding into dots
% carrierShift/N cycles per site moves the +1 order by exactly carrierShift rows.
carrier = 2*pi*carrierShift/N*Y;
phaseFull = gerchbergSaxton(gsTarget, double(written), initialPhase, gsIterations);
phaseYellow = gerchbergSaxton(gsTarget, double(yellow), initialPhase, gsIterations);
dotsFull = written & cos(carrier + phaseFull) <= 0;  % "0" sites of the binary hologram
labels = {'① 方形全孔径：全息图写满整个方形', '② 星夜样品：① 只留在黄色区域', '③ 星夜样品：按黄色区域设计'};
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
% sum(abs(F).^2)/N^2 in units of one site's incident power.
sites = sitesPerSide^2;
margin = round(0.02*N);
signalWindow = false(N);
signalWindow(max(min(rows)-margin, 1):min(max(rows)+margin, N), ...
    max(min(cols)-margin, 1):min(max(cols)+margin, N)) = true;
letterCore = erodeMask(letters, 3);  % keeps edge blur out of speckle statistics
meanLetterFull = mean(abs(farFields{1}(letters)).^2);
[transmitted, coverage, efficiency, brightness, contrast, correlation] = deal(zeros(1, 3));
for k = 1:3
    intensity = abs(farFields{k}).^2;
    transmitted(k) = (nnz(openSites{k}) - dotFill*nnz(dotSites{k}))/sites;
    coverage(k) = dotFill*nnz(dotSites{k})/nnz(hologramSites{k});
    efficiency(k) = sum(intensity(letters))/(N^2*sites);
    brightness(k) = mean(intensity(letters))/meanLetterFull;
    contrast(k) = std(intensity(letterCore))/mean(intensity(letterCore));
    c = corrcoef(intensity(signalWindow), double(letters(signalWindow)));
    correlation(k) = c(1, 2);
end
sizeMm = sitesPerSide*sitePitchUm/1000;
fprintf(['\n激光直写：%.2f mm 见方，%d×%d 个点位，%g µm 圆点、%g µm 间距（每个点位 %d×%d 像素），' ...
    '%d nm 蓝光，远场 "%s"\n天空 %.1f%%，黑色 %.1f%%，黄色 %.1f%% 的面积\n'], sizeMm, sitesPerSide, ...
    sitesPerSide, dotDiameterUm, sitePitchUm, samplesPerSite, samplesPerSite, wavelengthNm, ...
    strjoin(cellstr(targetText), ' / '), 100*nnz(siteClass == 1)/sites, ...
    100*nnz(siteClass == 2)/sites, 100*nnz(yellow)/sites);
for k = 1:3
    fprintf(['%s\n    透过光功率 %5.1f%%  全息区黄点覆盖 %4.1f%%  衍射效率 %.3f%%  ' ...
        '字母相对亮度 %.3f  散斑对比 %.2f  相关系数 %.3f\n'], labels{k}, 100*transmitted(k), ...
        100*coverage(k), 100*efficiency(k), brightness(k), contrast(k), correlation(k));
end

%% Figures
siteAxisMm = ((0:sitesPerSide-1) - sitesPerSide/2)*sitePitchUm/1000;
thetaDeg = asind(f*wavelengthNm/1000);
viewMargin = round(0.05*N);
viewHalfRows = min(max(rows) - (N/2 + 1) + viewMargin, N/2 - 1);  % text and, by symmetry, its twin
viewHalfCols = min(max(abs(cols - (N/2 + 1))) + viewMargin, N/2 - 1);
viewRows = N/2 + 1 + (-viewHalfRows:viewHalfRows);
viewCols = N/2 + 1 + (-viewHalfCols:viewHalfCols);
level = linspace(0, 1, 256)';
blueLight = [0.12*level, 0.35*level, level.^0.6];  % black to saturated blue, no white

% Plain 2-by-3 grid: in R2026a, row-spanning tiles let tile titles overlap the layout title.
fig = figure('Name', '激光直写星夜全息图', 'Position', [40 40 1400 860]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('激光直写 %g×%g mm：%g µm 圆点、%g µm 间距；%d nm 蓝光下远场为 "%s"', ...
    sizeMm, sizeMm, dotDiameterUm, sitePitchUm, wavelengthNm, strjoin(cellstr(targetText), ' ')));
subtitle(layout, {'上：近场（蓝 = 未加工透光，黑 = 黑色加工区，黄 = 写黄色圆点的点位；圆点实际大小见局部放大图）', ...
    '下：蓝光远场强度，上方 +1 级为文字，中心为 0 级，下方 −1 级为旋转 180° 的共轭像；色标上限为字母区平均强度的 1.5 倍，0 级饱和'});
for k = 1:3
    ax = nexttile(layout, k);
    image(ax, siteAxisMm, siteAxisMm, siteColours(openSites{k}(siteIdx, siteIdx), ...
        dotSites{k}(siteIdx, siteIdx), displayColours));
    axis(ax, 'xy'); axis(ax, 'image');
    xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
    title(ax, {labels{k}, sprintf('透过光功率 %.1f%%，全息区黄点覆盖 %.0f%%', ...
        100*transmitted(k), 100*coverage(k))});

    ax = nexttile(layout, 3 + k);
    intensity = abs(farFields{k}).^2;
    imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), intensity(viewRows, viewCols), ...
        [0, 1.5*mean(intensity(letters))]);
    axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, blueLight);
    xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
    title(ax, {sprintf('衍射效率 %.3f%%，相对亮度 %.3f', 100*efficiency(k), brightness(k)), ...
        sprintf('散斑对比 %.2f，与目标相关系数 %.3f', contrast(k), correlation(k))});
end

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

figure('Name', '激光直写：原图与圆点局部', 'Position', [80 120 1000 460]);
layout = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout);
pictureSizeMm = size(pictureRgb, [2 1])/max(size(pictureRgb, [1 2]))*sizeMm;
image(ax, pictureSizeMm(1)*[-0.5 0.5], pictureSizeMm(2)*[-0.5 0.5], flipud(pictureRgb));
axis(ax, 'xy'); axis(ax, 'image');
xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
title(ax, '原图：蓝 = 不加工，黑 = 黑色加工，黄 = 全息图（黄色圆点）');
ax = nexttile(layout);
image(ax, zoomAxisUm(zoomCols), zoomAxisUm(zoomRows), dotImage(openSites{3}(zoomRows, zoomCols), ...
    dotSites{3}(zoomRows, zoomCols), dotPixels, displayColours));
axis(ax, 'xy'); axis(ax, 'image');
xlabel(ax, 'x (µm)'); ylabel(ax, 'y (µm)');
title(ax, sprintf('③ 近场局部 %d×%d 个点位：每个 %g µm 圆点由 %d×%d 像素构成', zoomSites, zoomSites, ...
    dotDiameterUm, samplesPerSite, samplesPerSite));

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

function amplitude = renderText(textString, N, widthFraction)
% Rasterize bold text off-screen and centre it in an N-by-N array whose
% row index increases with y, so it reads upright with axis xy.
fig = figure('Visible', 'off', 'Color', 'k', 'Units', 'pixels', 'Position', [0 0 1000 400]);
ax = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off');
textObject = text(ax, 0.5, 0.5, textString, 'Units', 'normalized', 'Color', 'w', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontUnits', 'normalized', 'FontSize', 0.5, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
extent = textObject.Extent;  % shrink long or multi-line text to fit the canvas
if max(extent(3:4)) > 0.95
    textObject.FontSize = textObject.FontSize*0.95/max(extent(3:4));
end
frame = getframe(fig);
close(fig);
glyph = frame.cdata(:,:,1) > 127;
[rows, cols] = find(glyph);
glyph = double(flipud(glyph(min(rows):max(rows), min(cols):max(cols))));
width = round(widthFraction*N);
height = round(width*size(glyph, 1)/size(glyph, 2));
[xq, yq] = meshgrid(linspace(1, size(glyph, 2), width), linspace(1, size(glyph, 1), height));
amplitude = zeros(N);
amplitude(N/2 - floor(height/2) + (1:height), N/2 - floor(width/2) + (1:width)) = ...
    interp2(glyph, xq, yq, 'linear') >= 0.5;
end

function core = erodeMask(mask, radius)
% Pixels whose whole disk of the given radius lies inside the mask.
[dx, dy] = meshgrid(-radius:radius);
disk = double(hypot(dx, dy) <= radius);
core = conv2(double(mask), disk, 'same') > sum(disk(:)) - 0.5;
end
