# Starry Night near field (5×5 mm, scan lines) + Great Wave far field: hologram in both yellow and black

Scalar Fraunhofer simulation and laser writing plan for a 5 × 5 mm sample. The near field follows `starry_night_v2.png` in the repo root. Sky is left clear. Both yellow and black carry the same binary hologram (537,194 aperture sites); written sites in both regions are fully absorbing amplitude pixels. The two regions play the same optical role and differ only in how they are written. The far field is *The Great Wave off Kanagawa* (`great_wave.png`).

This is the both-hologram variant of the scanned-line model in `gs_great_wave_scanned_lines.m` (case ③), with the GS aperture extended from yellow-only to yellow + black.

Parameters: 1000×1000 sites, 5 µm pitch and line width, 405 nm, carrier 16.7 µm, 200 GS iterations, seed 1.

## Results (scalar Fraunhofer)

- dots: efficiency 2.73%, corr 0.680 (5×5 0.908)
- scan lines: efficiency 3.95%, corr 0.689 (5×5 0.909) (was 2.72% / 0.647 with black opaque)
- 268,241 written sites; 55,986 lines, median 20 µm, 3,760 at 5 µm, longest 215 µm

## Writing plans

Current file: `writing_plan_v2_starry_wave_5mm_scan_bothholo_zblack0p5.csv`. `make_plan.py` writes this file.

- path, origin (0, 0, 0), 0.03 mm/s for on and off, pause 0, power 40
- Yellow: preset 14, z = 0, written first (37,802 on-segments)
- Black: preset 1, z = +0.5 mm, written second (18,664 on-segments)
- Lines that cross the yellow/black boundary are split there

Earlier file, kept for reference: `writing_plan_v2_starry_wave_5mm_scan_bothholo.csv`. Same path, except black is at z = +1.5 mm (moved to +0.5 mm on 2026-11-02). Do not send that file to the stage.

Self-check on the current plan (`plan_check.txt`): schema/path OK, segments continuous, 1 preset switch, 1759 groups, 112,740 segments. Bounding box x 0–5.000, y 0–4.985 mm, z {0, 0.5}. Estimated time about 58.7 h; the returns to x = 0 between rows are not counted.

## Scripts and figures

- `gs_starry_wave_scan.py`: Python port of the case-③ scanned-line hologram with a yellow + black aperture. Loads the two pictures from the repo root and writes figures, site maps, `scan_segments.csv`, and `summary.json` into this folder. Needs NumPy, SciPy, Pillow, and Matplotlib.
- `make_plan.py`: builds the current writing plan from `site_map_1000x1000.npy` and `nearfield_class_1000x1000.npy`.
- `nearfield.png`, `site_map_1000x1000_black_is_dot.png`: where the laser writes. Yellow sites are preset 14; red sites in the near-field figure are preset 1 (the black paint of the picture).
- `farfield_sim_405nm.png`, `farfield_full_lines.png`: +1 order versus the Great Wave target, and the full field with the zero order and the −1 twin.

## Risks

- The model is a single scalar plane. The black layer sits 0.5 mm higher, which adds defocus and a phase or shear offset between the two sub-holograms. This is not modelled and can lower efficiency or contrast and blur the image (both layers are designed coherently as one hologram).
- Lines as short as 5 µm take only 167 ms to write at 0.03 mm/s, so stream jitter matters.
- Check that the parser accepts the carbide_preset column and the sign of z on the stage.
