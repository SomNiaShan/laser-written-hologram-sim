%% GS hologram whose far field reads STAR, then masked by a five-pointed star
% A phase-only Fourier hologram is designed with Gerchberg-Saxton (GS) so
% that its far field spells targetText. A star-shaped mask then blocks the
% hologram outside the star. Every hologram pixel contributes to every
% far-field point, so the masked hologram still reconstructs the text, only
% dimmer and with coarser speckle: its far field is the original far field
% convolved with the diffraction pattern of the star aperture. Case 3 reruns
% GS with the star as the hologram-plane support, i.e. designs for the mask.
%
% Scalar Fraunhofer model: uniform coherent plane-wave illumination,
% continuous phase, 100% pixel fill factor. The pixel sinc envelope and the
% higher diffraction orders of the pixel grid are ignored. The hologram is
% zero-padded by padFactor; padFactor = 2 samples the far-field intensity at
% its Nyquist rate, so GS cannot hide speckle between FFT samples.
% Uses base MATLAB only.

%% Parameters
targetText = 'STAR';
holoPixels = 512;          % hologram pixels per side
padFactor = 2;             % FFT size / hologram size
pixelPitchUm = 8;          % hologram pixel pitch
wavelengthNm = 532;
textWidthFraction = 0.6;   % text width / full far-field window
starWidthFraction = 0.98;  % star tip-to-tip width / hologram width
gsIterations = 200;
rng(1);                    % reproducible random initial phase

%% Hologram plane, star mask and target
N = padFactor*holoPixels;
holoIdx = N/2 - holoPixels/2 + (1:holoPixels);
[X, Y] = meshgrid((0:N-1) - N/2);  % pixel units; row index increases with y
squareAperture = zeros(N);
squareAperture(holoIdx, holoIdx) = 1;
starOuterRadius = starWidthFraction*holoPixels/2/cosd(18);
starMask = double(starPolygon(X, Y, starOuterRadius)) .* squareAperture;

targetAmplitude = renderText(targetText, N, textWidthFraction);
letters = targetAmplitude > 0;
initialPhase = 2*pi*rand(N);

%% Case 1: GS with the full square hologram
[phaseFull, farFull, errorFull] = gerchbergSaxton( ...
    targetAmplitude, squareAperture, initialPhase, gsIterations);

%% Case 2: block case 1 outside the star
farMasked = fft2c(starMask .* exp(1i*phaseFull));
errorMasked = amplitudeError(farMasked, targetAmplitude);

%% Case 3: GS with the star as the hologram aperture
[phaseStar, farStar, errorStar] = gerchbergSaxton( ...
    targetAmplitude, starMask, initialPhase, gsIterations);

%% Metrics
% Intensities share one scale: the same incident field, unnormalized FFT.
labels = {'① 方形全孔径 GS', '② 五角星掩膜遮挡 ①', '③ 以五角星为孔径的 GS'};
farFields = {farFull, farMasked, farStar};
apertures = {squareAperture, starMask, starMask};
phases = {phaseFull, phaseFull, phaseStar};
letterCore = erodeMask(letters, 6);  % keeps edge blur out of speckle statistics
meanLetterFull = mean(abs(farFull(letters)).^2);
[transmitted, efficiency, brightness, contrast, correlation] = deal(zeros(1, 3));
for k = 1:3
    intensity = abs(farFields{k}).^2;
    transmitted(k) = sum(apertures{k}(:))/sum(squareAperture(:));
    efficiency(k) = sum(intensity(letters))/sum(intensity(:));
    brightness(k) = mean(intensity(letters))/meanLetterFull;
    contrast(k) = std(intensity(letterCore))/mean(intensity(letterCore));
    c = corrcoef(intensity(:), double(letters(:)));
    correlation(k) = c(1, 2);
end
fprintf('\n远场 "%s"（%d×%d 全息图，%d 次 GS 迭代）\n', targetText, holoPixels, holoPixels, gsIterations);
for k = 1:3
    fprintf('%-14s 透过光功率 %5.1f%%  进入字母 %5.1f%%  字母相对亮度 %.2f  散斑对比 %.2f  相关系数 %.3f\n', ...
        labels{k}, 100*transmitted(k), 100*efficiency(k), brightness(k), contrast(k), correlation(k));
end

%% Figures
holoAxisMm = ((0:holoPixels-1) - holoPixels/2)*pixelPitchUm/1000;
thetaDeg = asind(((0:N-1) - N/2)/(N*pixelPitchUm)*wavelengthNm/1000);
viewRows = N/2 + 1 + (-round(0.27*N):round(0.27*N));
viewCols = N/2 + 1 + (-round(0.36*N):round(0.36*N));

% Plain 2-by-3 grid: in R2026a, row-spanning tiles let tile titles overlap the layout title.
fig = figure('Name', 'GS 全息图与五角星掩膜', 'Position', [40 40 1400 860]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('GS 相位全息图的远场为 "%s"：用五角星掩膜遮挡后依然是 "%s"', targetText, targetText));
subtitle(layout, '上：全息图相位（HSV 色环，黑色为遮挡）　下：远场强度（各图色标上限为字母区平均强度的 2 倍）');
for k = 1:3
    ax = nexttile(layout, k);
    image(ax, holoAxisMm, holoAxisMm, ...
        phaseImage(phases{k}(holoIdx, holoIdx), apertures{k}(holoIdx, holoIdx)));
    axis(ax, 'xy'); axis(ax, 'image');
    xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
    title(ax, {labels{k}, sprintf('透过光功率 %.1f%%', 100*transmitted(k))});

    ax = nexttile(layout, 3 + k);
    intensity = abs(farFields{k}).^2;
    imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), intensity(viewRows, viewCols), ...
        [0, 2*mean(intensity(letters))]);
    axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, hot);
    xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
    title(ax, {sprintf('相对亮度 %.2f，进入字母 %.1f%%', brightness(k), 100*efficiency(k)), ...
        sprintf('散斑对比 %.2f，与目标相关系数 %.3f', contrast(k), correlation(k))});
end

figure('Name', 'GS 目标、掩膜与收敛', 'Position', [80 80 1200 340]);
layout = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout);
imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), targetAmplitude(viewRows, viewCols));
axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, gray);
xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
title(ax, '目标远场');
ax = nexttile(layout);
imagesc(ax, holoAxisMm, holoAxisMm, starMask(holoIdx, holoIdx));
axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, gray);
xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
title(ax, sprintf('五角星掩膜：透过 %.1f%% 的面积', 100*transmitted(2)));
ax = nexttile(layout);
semilogx(ax, 1:gsIterations, errorFull, 1:gsIterations, errorStar, 'LineWidth', 1.5);
yline(ax, errorMasked, '--', '② 遮挡后（未再迭代）', 'LineWidth', 1.2);
grid(ax, 'on');
xlabel(ax, 'GS 迭代次数'); ylabel(ax, '归一化振幅误差');
legend(ax, '① 方形全孔径', '③ 五角星孔径', 'Location', 'southwest');
title(ax, 'GS 收敛');

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

function rgb = phaseImage(phase, aperture)
% Cyclic HSV colours for phase; blocked pixels (aperture = 0) are black.
index = round(mod(phase, 2*pi)/(2*pi)*255) + 1;
rgb = ind2rgb(index, hsv(256)) .* aperture;
end
