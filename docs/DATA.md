# Included data

## Source metadata

The included fields were extracted from the source image described by `Lam_065_metadata.json`:

- source file: `Lam_065.czi`;
- channel: 0;
- source dimensions: 59,208 × 49,805 pixels;
- source pixel size: 0.44 µm in X and Y;
- patch size: 4096 × 4096 pixels;
- stride: 3800 pixels;
- overlap: 296 pixels;
- project metadata: RockFace;
- institution metadata: LCCMat / Petrobras;
- engine metadata: Alfa 1.0.

## Included fields

| Label | File | Patch origin |
|---|---|---|
| R1 | `structures/patch_y3800_x3800_c0_mask.png` | y=3800, x=3800 |
| R2 | `structures/patch_y7600_x19000_c0_mask.png` | y=7600, x=19000 |
| R3 | `structures/patch_y7600_x53200_c0_mask.png` | y=7600, x=53200 |

These are three fields from one source image, not three independent lithologies.

## Mask convention

The supplied masks are RGB images. The reference loader classifies a pixel as **solid** when green exceeds both red and blue by more than 20 intensity units; the pore mask is the complement.

## Reference scale

At 256 × 256 pixels, the effective pixel size is 7.04 µm and the field side is 1.80224 mm.

## Provenance and reuse

The software is MIT-licensed. Image/data reuse may be subject to separate provenance, institutional, contractual, or publication permissions. Verify applicable permissions before redistributing source imagery or derivative datasets.
