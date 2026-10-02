# Data directory

This directory contains the acquisition metadata and the three segmented rock-section fields used by the reference JTPD analysis.

- `Lam_065_metadata.json` — physical calibration and patch-generation metadata.
- `structures/patch_y3800_x3800_c0_mask.png` — R1.
- `structures/patch_y7600_x19000_c0_mask.png` — R2.
- `structures/patch_y7600_x53200_c0_mask.png` — R3.

The three masks are fields from the same source image, not independent lithologies.

The reference loader interprets green-dominant pixels as the solid phase and uses the complement as pore space. At the main 256 × 256 analysis scale, the effective pixel size is 7.04 µm and the field side is 1.80224 mm.

See [../docs/DATA.md](../docs/DATA.md) for full metadata, segmentation convention, and provenance notes.
