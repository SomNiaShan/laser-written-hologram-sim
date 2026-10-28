%% Binary (0/1) hologram reading "Vincent van Gogh", masked by Starry Night
% Same method as gs_star_binary_hologram.m: GS designs a phase phi(x,y) for
% the text, then the phase is written into a 0/1 carrier grating
%   t = aperture .* (cos(2*pi*y/P + phi) > 0),
% whose +1 order shows the text, -1 order its 180-degree-rotated twin and
% zero order the aperture's own diffraction pattern. Here the mask is a
% picture and only its yellow pixels transmit. The yellow regions are many
% separate stars, swirls and dashes, so the mask diffracts into a broad
% halo instead of a compact spot: behind the mask much of the light forms a
% haze around the letters, and the zero order becomes a bright cloud. A
% short carrier period keeps the text clear of that cloud. Case 3 designs
% the GS phase for the yellow aperture before encoding.
%
% Model assumptions are those of gs_star_hologram.m: scalar Fraunhofer
% diffraction, uniform coherent plane-wave illumination, 100% pixel fill
% factor (pixel sinc envelope and pixel-grid orders ignored), hologram
% zero-padded by padFactor. Uses base MATLAB only.

%% Parameters
targetText = 'Vincent van Gogh';  % a cell array of strings gives several lines
maskImageFile = fullfile(fileparts(mfilename('fullpath')), 'starry_night.png');
holoPixels = 512;           % hologram pixels per side
padFactor = 2;              % FFT size / hologram size
pixelPitchUm = 8;           % hologram pixel pitch
wavelengthNm = 532;
carrierPeriodPixels = 4;    % binary fringe period; sets the text's distance from the zero order
textWidthFraction = 0.8;    % text width / full far-field window
gsIterations = 200;
rng(1);                     % reproducible random initial phase

%% Hologram plane, picture mask and target
N = padFactor*holoPixels;
holoIdx = N/2 - holoPixels/2 + (1:holoPixels);
[~, Y] = meshgrid((0:N-1) - N/2);  % pixel units; row index increases with y
squareAperture = zeros(N);
squareAperture(holoIdx, holoIdx) = 1;
maskRgb = imread(maskImageFile);
yellow = maskRgb(:,:,1) > 127 & maskRgb(:,:,2) > 127 & maskRgb(:,:,3) < 128;
pictureMask = zeros(N);
pictureMask(holoIdx, holoIdx) = fitToHologram(yellow, holoPixels);

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
phasePicture = gerchbergSaxton(textAmplitude, pictureMask, initialPhase, gsIterations);
binaryFull = squareAperture .* (cos(carrier + phaseFull) > 0);         % case 1
binaryMasked = pictureMask .* binaryFull;                              % case 2: case 1 behind the mask
binaryPicture = pictureMask .* (cos(carrier + phasePicture) > 0);      % case 3

%% Metrics
% Powers are relative to the light incident on the square hologram; the
% unnormalized FFT gives sum(abs(F(:)).^2) = N^2*sum(abs(t(:)).^2).
labels = {'① 方形全孔径 GS → 二值编码', '② 星夜黄色掩膜遮挡 ①', '③ 以黄色区域为孔径的 GS → 二值编码'};
holograms = {binaryFull, binaryMasked, binaryPicture};
farFields = cellfun(@fft2c, holograms, 'UniformOutput', false);
incident = sum(squareAperture(:));
margin = round(0.02*N);
signalWindow = false(N);
signalWindow(max(min(rows)-margin, 1):min(max(rows)+margin, N), ...
    max(min(cols)-margin, 1):min(max(cols)+margin, N)) = true;
letterCore = erodeMask(letters, 3);  % keeps edge blur out of speckle statistics
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
fprintf('\n0/1 二值全息图，远场 "%s"（%d×%d 像素，载频周期 %d 像素，%d 次 GS 迭代），黄色掩膜透过 %.1f%% 面积\n', ...
    strjoin(cellstr(targetText), ' / '), holoPixels, holoPixels, carrierPeriodPixels, gsIterations, ...
    100*sum(pictureMask(:))/incident);
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
fig = figure('Name', '0/1 二值全息图与星夜掩膜', 'Position', [40 40 1400 860]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('0/1 二值振幅全息图的远场为 "%s"：只让星夜图的黄色区域透光后依然可读', ...
    strjoin(cellstr(targetText), ' ')));
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

figure('Name', '0/1 二值全息图：原图、掩膜与目标', 'Position', [80 120 1200 380]);
layout = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout);
pictureSizeMm = size(maskRgb, [2 1])/max(size(maskRgb, [1 2]))*holoPixels*pixelPitchUm/1000;
image(ax, pictureSizeMm(1)*[-0.5 0.5], pictureSizeMm(2)*[-0.5 0.5], flipud(maskRgb));
axis(ax, 'xy'); axis(ax, 'image');
xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
title(ax, '原图：只有黄色区域透光');
ax = nexttile(layout);
imagesc(ax, holoAxisMm, holoAxisMm, pictureMask(holoIdx, holoIdx));
axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, gray);
xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
title(ax, sprintf('透光掩膜（%d×%d 像素）：透过 %.1f%% 的面积', holoPixels, holoPixels, ...
    100*sum(pictureMask(:))/incident));
ax = nexttile(layout);
imagesc(ax, thetaDeg(viewCols), thetaDeg(viewRows), double(letters(viewRows, viewCols)));
axis(ax, 'xy'); axis(ax, 'image'); colormap(ax, gray);
xlabel(ax, '\theta_x (°)'); ylabel(ax, '\theta_y (°)');
title(ax, '目标：+1 级处的文字');

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

function mask = fitToHologram(pixels, holoPixels)
% Fit a logical picture into holoPixels-by-holoPixels hologram pixels,
% keeping its aspect ratio (centred, zero outside). Each hologram pixel
% takes the majority of the picture area it covers. Rows are flipped so
% the picture reads upright with axis xy.
pixels = double(flipud(pixels));
scale = max(size(pixels))/holoPixels;  % picture pixels per hologram pixel
k = 2*round((scale - 1)/2) + 1;        % odd box close to one hologram pixel
pixels = conv2(pixels, ones(k)/k^2, 'same');
offsets = ((1:holoPixels) - (holoPixels + 1)/2)*scale;
[xq, yq] = meshgrid((size(pixels, 2) + 1)/2 + offsets, (size(pixels, 1) + 1)/2 + offsets);
mask = interp2(pixels, xq, yq, 'linear', 0) >= 0.5;
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
