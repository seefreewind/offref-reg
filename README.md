# OFFREF-REG

A frozen reliability audit of AlphaGenome and Borzoi predictions for complete human non-reference insertions.

## Reproduce figures and check the frozen inputs

Python 3.10+:

```sh
python -m pip install -r requirements.txt
python scripts/verify_release.py
python scripts/build_figures.py
```

The figure builder reads frozen values; it does not recompute primary statistics or confidence intervals. It exports six main figures and two supplementary figures. Arial is used when installed, with DejaVu Sans as a portable fallback. Numeric source files retain full precision.

## Contents

- `source_data/`: all main-panel inputs, 400 core identities, 500 reference identities, robustness/member details, frozen endpoint summaries and panel mapping.
- `tables/`, `supplementary/`: editable numeric main tables and relocated supplementary inputs.
- `scripts/`: portable figure composition and integrity/endpoint checks.
- `analysis_archive/`: historical input reconstruction, scoring, statistics and audit code. Paths retain the executed phase layout. These are provenance source files; rerunning inference requires separately obtained sequences/annotations, official model dependencies, weights and a personal AlphaGenome API key. This release is sufficient for figure and central endpoint verification, not a turnkey re-download of all third-party inputs.
- `SHA256SUMS.json`: hashes of every committed release file except the manifest itself.

## Frozen evidence

400 loci / 387 genes; 318 jointly finite and 82 without common scorable exon support. Off-reference cross-model rho −0.1113027120530779, sign concordance 149/318. Reference panel 500 variants / 473 genes. GNRS biological calibration remains descriptive at 9/33 usable loci. The reference comparator does not isolate reference absence as a causal effect.

## Terms and attribution

**Read [DATA_TERMS.md](DATA_TERMS.md) before using data or figures. AlphaGenome outputs and derivatives are subject to the [AlphaGenome Output Terms of Use](https://deepmind.google.com/science/alphagenome/output-terms), including non-commercial restrictions. The MIT license covers original code only.** No credentials, raw donor data, third-party sequences or model weights are packaged.

The archived creator is the project entity **OFFREF-REG**. Paper authorship is recorded separately in the manuscript.
