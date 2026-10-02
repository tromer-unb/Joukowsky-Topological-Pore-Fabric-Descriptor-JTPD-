# Examples

## Reproduce the included fields

```bash
python jtpd_analysis.py
```

## Different fields from the same acquisition

Edit the `MASKS` and `NAMES` lists in `jtpd_analysis.py` while keeping the acquisition metadata unchanged.

## Different acquisition

Also update the physical calibration used by `source_pixel_um()` and adapt `load_solid_mask()` if your segmentation encoding differs.

For scientific comparisons, control or explicitly report target resolution and field size.
