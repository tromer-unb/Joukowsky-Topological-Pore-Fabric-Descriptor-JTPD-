# Reproducibility protocol

## Environment

Recommended: Python 3.10+ in a clean virtual environment.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
python jtpd_analysis.py
```

On Windows PowerShell activate with `.\.venv\Scripts\Activate.ps1`.

## Expected reference results

| sample | effective pixel (µm) | pore fraction | β₀ | β₁ | χ | X span | Y span |
|---|---:|---:|---:|---:|---:|---:|---:|
| R1 | 7.04 | 0.907211 | 17 | 255 | -238 | 1 | 1 |
| R2 | 7.04 | 0.730698 | 43 | 135 | -92 | 1 | 1 |
| R3 | 7.04 | 0.465149 | 37 | 79 | -42 | 0 | 1 |

Small floating-point differences can occur across library versions; integer topology and spanning results should remain stable.

## Built-in validation

The script validates a disk, annulus, two disjoint disks, and a horizontal spanning channel. It also checks the Joukowsky response under translation, scaling, rotation, reflection, and contour reversal.

## Automated check

`.github/workflows/reproducibility.yml` installs the dependencies and reruns the full reference workflow on pushes and pull requests.

## Checklist for new datasets

Archive the exact masks, physical calibration, code commit/release, generated CSV/JSON outputs, effective analysis scale, and any changed segmentation/connectivity conventions.
