# Ready-to-run bundle

`JTPD_descriptor_and_rock_fields_v0.1.0.zip` is a compact research bundle containing:

- `jtpd_analysis.py`
- `rock_joukowsky.py`
- `requirements.txt`
- the three R1–R3 segmented masks
- acquisition metadata
- reference CSV/JSON results

After downloading and extracting:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python jtpd_analysis.py
```

For the complete documentation and manuscript files, use the full GitHub repository.
