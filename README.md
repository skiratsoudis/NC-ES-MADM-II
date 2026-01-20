# NC-ES-MADM-II
Reproducible computational platform for NC-ES-MADM II: entropy–reliability integrated weighting, supervised calibration, robustness diagnostics, and manuscript figure generation.
---

## Relation to the manuscript
This repository constitutes the official computational platform supporting the empirical results, robustness diagnostics, and figure generation reported in the accompanying manuscript.  
All tables and figures referenced in the paper can be reproduced using the scripts in `code/` and `plots/` with the input datasets in `data/raw/` and the corresponding outputs in `results/`.

## Versioning and archival
A citable DOI will be issued through Zenodo upon the first tagged release (v1.0.0) and should be used when citing this software artifact.

## Scientific scope
NC-ES-MADM II is an entropy–reliability integrated multi-attribute decision-making framework designed for robust benchmarking and supervised decision calibration.  
The platform combines objective entropy-based weighting, reliability-aware subjective aggregation, probabilistic decision distributions, and system-level diagnostics for coupling, decisiveness, and stability.

## Computational characteristics
- Deterministic implementation (no stochastic components)
- No hyperparameter tuning
- Simplex-preserving weight aggregation
- Fully auditable intermediate quantities (weights, posteriors, diagnostics)
- Designed for transparency rather than black-box optimization

## Intended use
This repository is intended for:
- Peer review and reproducibility assessment
- Methodological validation and benchmarking
- Educational and research use in decision analytics
- Extension to logistics, policy, and socio-economic evaluation contexts
