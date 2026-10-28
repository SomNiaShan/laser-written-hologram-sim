%% Binary (0/1) amplitude hologram whose far field reads STAR, masked by a star
% Each hologram pixel transmits either 0 or 1. GS first designs a phase
% phi(x,y) whose far field is the text, exactly as in gs_star_hologram.m.
% The phase is then written into a binary carrier grating, the same
% encoding as a fork grating:
%   t = aperture .* (cos(2*pi*y/P + phi) > 0)
% Per pixel, t = 1/2 + (1/pi)*(exp(1i*(2*pi*y/P + phi)) + c.c.) - ..., so
% the +1 order reproduces the GS far field shifted by the carrier frequency
% 1/P (N/P far-field samples), the -1 order is
% its 180-degree-rotated conjugate (twin), the zero order is the aperture's
% diffraction pattern, and the +-3, +-5 ... orders scatter weakly. Any real
% 0/1 transmission has a centrosymmetric far field, so the twin cannot be
% removed, and at most (1/pi)^2 ~ 10% of the incident power reaches the text.
% (Enforcing 0/1 inside the GS loop gave more speckle than this encoding.)
%
% The star mask blocks pixels outside a five-pointed star; the masked
% hologram stays 0/1. Case 3 designs the GS phase for the star aperture
% before encoding. Model assumptions are those of gs_star_hologram.m:
% scalar Fraunhofer diffraction, uniform coherent plane-wave illumination,
% 100% pixel fill factor (pixel sinc envelope and pixel-grid orders ignored),
% hologram zero-padded by padFactor. Uses base MATLAB only.

%% Parameters
targetText = 'STAR';
holoPixels = 512;           % hologram pixels per side
padFactor = 2;              % FFT size / hologram size
pixelPitchUm = 8;           % hologram pixel pitch
wavelengthNm = 532;
carrierPeriodPixels = 8;    % binary fringe period; sets the text's distance from the zero order
textWidthFraction = 0.5;    % text width / full far-field window
starWidthFraction = 0.98;   % star tip-to-tip width / hologram width
gsIterations = 200;
rng(1);                     % reproducible random initial phase

%% Hologram plane, star mask and target
N = padFactor*holoPixels;
holoIdx = N/2 - holoPixels/2 + (1:holoPixels);
[X, Y] = meshgrid((0:N-1) - N/2);  % pixel units; row index increases with y
squareAperture = zeros(N);
squareAperture(holoIdx, holoIdx) = 1;
starOuterRadius = starWidthFraction*holoPixels/2/cosd(18);
starMask = double(starPolygon(X, Y, starOuterRadius)) .* squareAperture;

textAmplitude = renderText(targetText, N, textWidthFraction);  % centred GS target
carrierShift = round(N/carrierPeriodPixels);  % far-field rows from the zero order to the text
letters = circshift(textAmplitude > 0, [carrierShift, 0]);     % text in the +1 order
[rows, cols] = find(letters);
if min(rows) <= N/2 + 1
    warning('The text reaches the zero order and its twin; reduce carrierPeriodPixels or textWidthFraction.');
end
initialPhase = 2*pi*rand(N);

%% GS phases, then binary carrier encoding
% carrierShift/N cycles per pixel moves the +1 order by exactly carrierShift rows.
carrier = 2*pi*carrierShift/N*Y;
phaseFull = gerchbergSaxton(textAmplitude, squareAperture, initialPhase, gsIterations);
phaseStar = gerchbergSaxton(textAmplitude, starMask, initialPhase, gsIterations);
binaryFull = squareAperture .* (cos(carrier + phaseFull) > 0);  % case 1
binaryMasked = starMask .* binaryFull;                          % case 2: case 1 behind the mask
binaryStar = starMask .* (cos(carrier + phaseStar) > 0);        % case 3

%% Metrics
% Powers are relative to the light incident on the square hologram; the
% unnormalized FFT gives sum(abs(F(:)).^2) = N^2*sum(abs(t(:)).^2).
labels = {'① 方形全孔径 GS → 二值编码', '② 五角星掩膜遮挡 ①', '③ 以五角星为孔径的 GS → 二值编码'};
holograms = {binaryFull, binaryMasked, binaryStar};
farFields = cellfun(@fft2c, holograms, 'UniformOutput', false);
incident = sum(squareAperture(:));
margin = round(0.02*N);
signalWindow = false(N);
signalWindow(max(min(rows)-margin, 1):min(max(rows)+margin, N), ...
    max(min(cols)-margin, 1):min(max(cols)+margin, N)) = true;
letterCore = erodeMask(letters, 6);  % keeps edge blur out of speckle statistics
meanLetterFull = mean(abs(farFields{1}(letters)).^2);
[transmitted, efficiency, brightness, contrast, correlation] = deal(zeros(1, 3));
for k = 1:3
    intensity = abs(farFields{k}).^2;
    transmitted(k) = sum(holograms{k}(:))/incident;
    efficiency(k) = sum(intensity(letters))/(N^2*incident);
    brightness(k) = mean(intensity(letters))/meanLetterFull;
    contrast(k) = std(intensity(letterCore))/mean(intensity(letterCore));
    c = corrcoef(intensity(signalWindow), double(letters(signalWindow)));
    correlation(k) = c(1, 2);
end
fprintf('\n0/1 二值全息图，远场 "%s"（%d×%d 像素，载频周期 %d 像素，%d 次 GS 迭代）\n', ...
    targetText, holoPixels, holoPixels, carrierPeriodPixels, gsIterations);
for k = 1:3
    fprintf('%s\n    透过光功率 %5.1f%%  衍射效率 %5.2f%%  字母相对亮度 %.2f  散斑对比 %.2f  相关系数 %.3f\n', ...
        labels{k}, 100*transmitted(k), 100*efficiency(k), brightness(k), contrast(k), correlation(k));
end

%% Figures
holoAxisMm = ((0:holoPixels-1) - holoPixels/2)*pixelPitchUm/1000;
thetaDeg = asind(((0:N-1) - N/2)/(N*pixelPitchUm)*wavelengthNm/1000);
viewMargin = round(0.05*N);
viewHalfRows = min(max(rows) - (N/2 + 1) + viewMargin, N/2 - 1);  % text and, by symmetry, its twin
viewHalfCols = min(max(abs(cols - (N/2 + 1))) + viewMargin, N/2 - 1);
viewRows = N/2 + 1 + (-viewHalfRows:viewHalfRows);
viewCols = N/2 + 1 + (-viewHalfCols:viewHalfCols);

% Plain 2-by-3 grid: in R2026a, row-spanning tiles let tile titles overlap the layout title.
fig = figure('Name', '0/1 二值全息图与五角星掩膜', 'Position', [40 40 1400 860]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('0/1 二值振幅全息图的远场为 "%s"：用五角星掩膜遮挡后依然是 "%s"', targetText, targetText));
subtitle(layout, {'上：全息图透过率（白 = 1 透光，黑 = 0 不透光）', ...
    '下：远场强度，上方 +1 级为文字，中心为 0 级，下方 −1 级为旋转 180° 的共轭像；色标上限为字母区平均强度的 2 倍，0 级饱和'});
for k = 1:3
    ax = nexttile(layout, k);
    imagesc(ax, holoAxisMm, holoAxisMm, holograms{k}(holoIdx, holoIdx), [0 1]);
    axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, gray);
    xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
    title(ax, {labels{k}, sprintf('透过光功率 %.1f%%', 100*transmitted(k))});

    ax = nexttile(layout, 3 + k);
    intensity = abs(farFields{k}).^2;
    imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), intensity(viewRows, viewCols), ...
        [0, 2*mean(intensity(letters))]);
    axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, hot);
    xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
    title(ax, {sprintf('衍射效率 %.2f%%，相对亮度 %.2f', 100*efficiency(k), brightness(k)), ...
        sprintf('散斑对比 %.2f，与目标相关系数 %.3f', contrast(k), correlation(k))});
end

figure('Name', '0/1 二值全息图：目标、掩膜与局部', 'Position', [80 120 1200 360]);
layout = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout);
imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), double(letters(viewRows, viewCols)));
axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, gray);
xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
title(ax, '目标：+1 级处的文字');
ax = nexttile(layout);
imagesc(ax, holoAxisMm, holoAxisMm, starMask(holoIdx, holoIdx));
axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, gray);
xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
title(ax, sprintf('五角星掩膜：透过 %.1f%% 的面积', 100*sum(starMask(:))/incident));
ax = nexttile(layout);
zoomPixels = holoPixels/2 + (-31:32);
imagesc(ax, 1000*holoAxisMm(zoomPixels), 1000*holoAxisMm(zoomPixels), ...
    binaryFull(holoIdx(zoomPixels), holoIdx(zoomPixels)), [0 1]);
axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, gray);
xlabel(ax, 'x (µm)'); ylabel(ax, 'y (µm)');
title(ax, {'① 中心 64×64 像素', sprintf('周期 %d 像素的载频条纹被 GS 相位扭曲', carrierPeriodPixels)});

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

function inside = starPolygon(X, Y, outerRadius)
% Regular five-pointed star, one tip up, bounding box centred at the origin.
vertexAngle = 90 + 36*(0:9);
radius = outerRadius*repmat([1, cosd(72)/cosd(36)], 1, 5);
vx = radius.*cosd(vertexAngle);
vy = radius.*sind(vertexAngle);
vy = vy - (max(vy) + min(vy))/2;
inside = inpolygon(X, Y, vx, vy);
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
