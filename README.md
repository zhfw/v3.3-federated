# CFPS Manufacturability Predictor v3.3 Integrated

This release merges the trained v3.3 backend with the earlier v2.5 diagnostic/release/reporting workflow.

## Main additions
- v2.5 compatibility/release layer: 70% protein + 30% DNA.
- GO only when combined >=85, protein >=72, DNA >=78.
- Red warning below combined score 85.
- Structured warning -> likely consequence -> possible solution -> validation experiment.
- Downloadable PDF report.
- v3.3 trained CFPS solubility, aggregation, RP3 production, general solubility, and NESG heads.
- Federated network tab retained.

## Important interpretation
The v2.5 compatibility scores are a transparent release layer built from current v3.3 endpoints. They are not a separately retrained v2.5 ML model. Absolute CFPS yield and CHO secretion/titer are not predicted in this release.

## Run
Double-click `run_webpage.command`, or:

```bash
python3 -m pip install -r requirements.txt
python3 -m streamlit run app.py
```
