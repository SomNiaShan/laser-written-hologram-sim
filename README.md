# GS-star-hologram

MATLAB simulations of binary Fourier holograms, designed with the Gerchberg–Saxton (GS) algorithm, that still show their far-field image after being cut to the shape of a picture.

A Fourier hologram spreads every far-field point over the whole hologram. A hologram masked by a star, or kept only in the yellow parts of *The Starry Night*, therefore still reconstructs its image, though dimmer, with more speckle, and hazy when the mask is fragmented. Rerunning GS with the mask as the hologram-plane support recovers most of the detail. The laser-writing scripts model a sample made by laser direct writing, written either as single dots or as continuous lines.

## Main files

The scripts build on each other, in this order:

- `gs_star_hologram.m`: phase-only hologram (512 × 512 pixels) whose far field reads "STAR", masked by a five-pointed star.
- `gs_star_binary_hologram.m`: the same with a 0/1 amplitude hologram. The GS phase φ is encoded as the carrier grating `cos(2πy/P + φ) > 0`, so the text appears in the +1 diffraction order and its 180°-rotated twin in the −1 order.
- `gs_starry_night_binary_hologram.m`: only the yellow pixels of `starry_night.png` transmit; the far field reads "Vincent van Gogh".
- `gs_starry_night_laser_hologram.m`: laser-written 5 × 5 mm sample with 1000 × 1000 sites at a 5 µm pitch, following `starry_night_v2.png`. The sky is left clear, the black parts block light and the yellow parts carry the hologram as 5 µm round dots, each drawn with 20 × 20 pixels. The 405 nm far field reads "Vincent van Gogh".
- `gs_starry_night_great_wave_hologram.m`: the same sample, but the far field shows *The Great Wave off Kanagawa* (`great_wave.png`) instead of text. Its figures carry English captions that explain the method and each panel.
- `gs_great_wave_scanned_lines.m`: the same holograms as the laser writes them by scanning along x and switching on and off. Each run of neighbouring dot sites in a row becomes one continuous line with round ends (a stadium of width `lineWidthUm`, 5 µm by default), with no gaps between sites. The far field stays exact, because every dot site is one of four tiles (single dot, left end, middle, right end). It prints dots and lines side by side, the number of lines and their lengths, and draws close-ups and far fields. For ③ the efficiency rises from 1.91% (dots) to 2.72% (5 µm lines), with the same picture quality (correlation 0.872 → 0.878 after 5 × 5 averaging); there are 38,191 lines, median 4 sites (20 µm), 2,892 of them a single site.
- `gs_line_width_sweep.m`: the scanned-line model for line widths of 3–5 µm, for ① and ③. The efficiency grows about linearly with the width (③: 1.22% at 3 µm, 1.95% at 4 µm, 2.72% at 5 µm), while the correlation saturates from about 4 µm (③: 0.834, 0.873, 0.878). Widths above the 5 µm pitch, where neighbouring rows overlap, are not modelled yet.

Each script compares three cases: ① a hologram over the whole square, ② hologram ① kept only inside the mask, ③ a hologram designed for the mask.

## Input pictures

- `starry_night.png`: mask for `gs_starry_night_binary_hologram.m`; its yellow pixels transmit.
- `starry_night_v2.png`: three-colour near-field design (blue sky, black, yellow) for the two laser-writing scripts.
- `great_wave.png`: far-field target; its brightness is used, with the dark background set to zero.

## Requirements

- MATLAB only; no toolbox is needed (checked with `matlab.codetools.requiredFilesAndProducts`). Tested with MATLAB R2026a.
- Each script loads its pictures from its own folder, so it can be run from anywhere. The two scanned-line scripts first run `gs_starry_night_great_wave_hologram.m`, so keep them in the same folder (or on the MATLAB path).

## Notes

- The first four scripts label their figures and console output in Chinese; `gs_starry_night_great_wave_hologram.m` uses English. Code comments are in English throughout.
- Model: scalar Fraunhofer diffraction of a uniform, coherent plane wave. In the laser-writing scripts the far field of the full pixel raster is computed exactly with the convolution theorem, `F = C(f)·DFT(clear sites) − D(f)·DFT(dot sites)`, where C and D are the Fourier transforms of one clear site and one dot. The scanned-line scripts replace D by four tile transforms, one per kind of line piece.
