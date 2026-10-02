# Contributing

Contributions improving reproducibility, documentation, testing, portability, or scientifically justified extensions are welcome.

Before proposing a change:

- keep the 2-D scope explicit;
- preserve physical-scale calibration;
- do not silently change digital-connectivity conventions;
- update validation when changing topology or Joukowsky routines;
- keep generated numerical outputs traceable to code and input masks.

Recommended workflow: fork, create a focused branch, install `requirements.txt`, run `python jtpd_analysis.py`, verify validation passes, then open a pull request describing the scientific and numerical effect.
