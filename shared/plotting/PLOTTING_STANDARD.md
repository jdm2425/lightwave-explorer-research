# Laser Physics PhD plotting standard

The attached reference figure was produced at approximately the original
Matplotlib size used by Jack: 3.4 x 2.6 inches. That is the permanent default
for ordinary single-panel report plots.

## Default
- figure: 3.4 x 2.6 in
- base / axis-label font: 7 pt
- tick + legend font: 6 pt
- line width: 1.5 pt
- marker size: 4 pt
- axis width: 0.8 pt
- inward ticks on all four sides
- dotted grid, 0.5 pt, alpha 0.6
- 300 dpi

## Wide figures
Dense spectral overlays may use 6.8 x 3.6 in, but must retain the same physical
font, line, tick, and axis styling.

## Important rule
Do not set large per-script figure sizes such as 7.2 x 4.6 or 7.4 x 4.8 in for
ordinary summary plots. Those sizes were the reason the labels and lines looked
too small relative to the axes in the pump-GDD comparison.
