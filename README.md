# NC-ES-MADM II — Reproducible Computational Platform  
**Calibrated Neural-Consensus Entropy–Synergy Multicriteria Ranking with Robustness Diagnostics (NC-ES-MADM II)**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18318267.svg)](https://doi.org/10.5281/zenodo.18318267)  
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## Manuscript
This repository provides the **official reproducible computational implementation** of the **NC-ES-MADM II** (Neural-Consensus Entropy–Synergy Multi-Attribute Decision-Making) framework, introduced in the manuscript:

**Calibrated Neural-Consensus Entropy–Synergy Multicriteria Ranking with Robustness Diagnostics: The NC-ES-MADM II Framework**

The platform supports **full reproducibility** of the empirical results, robustness diagnostics, and manuscript figure generation, through a transparent and auditable workflow.

---

## Authors and affiliation
**Sideris Kiratsoudis¹\***, **Vasilis Tsiantos¹**  
¹ Department of Physics, Faculty of Sciences, Democritus University of Thrace, Campus St. Lucas, 65404 Kavala, Greece

**Corresponding author:**  
Lt Col (Ordnance) Dr Sideris Kiratsoudis  
Postdoctoral Researcher, Department of Physics, Faculty of Sciences  
Democritus University of Thrace, Campus St. Lucas, 65404 Kavala, Greece  
Email: skyratso@physics.duth.gr  
ORCID: https://orcid.org/0009-0006-0584-2751

---

## What the framework does
NC-ES-MADM II integrates, within a unified decision-support pipeline:

- **Entropy-based objective weighting**
- **Reliability-aware expert consensus aggregation**
- **Optional supervised calibration** (when a credible target signal exists)
- **Robustness and stability diagnostics** for decision auditability

The design prioritizes **interpretability, traceability, and reproducibility**, rather than black-box optimization.

---

## Repository structure
- `code/` — Core NC-ES-MADM II implementation (R scripts)  
- `plots/` — Manuscript figure-generation scripts (R scripts)  
- `data/raw/` — Input datasets used in the empirical case study  
- `results/` — Reproduced outputs (tables and derived results)

---

## Reproducibility protocol
A formal step-by-step reproduction guide is provided in:

- `REPRODUCIBILITY.md`

In brief, the workflow is:

1. Use input datasets in `data/raw/`  
2. Run the main implementation in `code/` to reproduce tables and decision outputs into `results/`  
3. Run scripts in `plots/` to regenerate manuscript figures from the produced outputs

---

## Computational characteristics
- Deterministic implementation (no stochastic components)
- No hyperparameter tuning
- Simplex-preserving weight aggregation
- Fully auditable intermediate quantities (weights, posteriors, diagnostics)
- Transparency-driven architecture for peer-review reproducibility

---

## Intended use
This repository is intended for:

- Peer review and reproducibility assessment
- Methodological validation and benchmarking
- Educational and research use in decision analytics
- Extensions to logistics, policy, and socio-economic evaluation contexts

---

## How to cite
If you use this repository, please cite the archived software record:

Kiratsoudis, S. (2026). *NC-ES-MADM II: Reproducible Computational Platform* (Version 1.0.0). Zenodo. https://doi.org/10.5281/zenodo.18318267

---

## License
This repository is released under the **MIT License**. See `LICENSE`.

