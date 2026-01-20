[![Release](https://img.shields.io/github/v/release/skiratsoudis/NC-ES-MADM-II)](https://github.com/skiratsoudis/NC-ES-MADM-II/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-100%25-blue)](https://www.r-project.org/)


# Calibrated Neural-Consensus Entropy–Synergy Multicriteria Ranking  
## The NC-ES-MADM II Framework

**Authors**

Sideris Kiratsoudis¹*, Vasilis Tsiantos¹  

¹ Department of Physics, Faculty of Sciences,  
Democritus University of Thrace,  
Campus St. Lucas, 65404 Kavala, Greece  

\* Corresponding author

**Corresponding author**

Lt Col (Ordnance) Dr Sideris Kiratsoudis  
Postdoctoral Researcher, Department of Physics, Faculty of Sciences  
Democritus University of Thrace, Campus St. Lucas, 65404 Kavala, Greece  

Email: skyratso@physics.duth.gr  
ORCID: https://orcid.org/0009-0006-0584-2751

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18318267.svg)](https://doi.org/10.5281/zenodo.18318267)

This repository hosts the reproducible computational implementation of the  
**NC-ES-MADM II (Neural-Consensus Entropy–Synergy Multi-Attribute Decision-Making) framework**,  
as introduced in the manuscript:

*Calibrated Neural-Consensus Entropy–Synergy Multicriteria Ranking with Robustness Diagnostics.*

The framework integrates entropy-based objective weighting, reliability-aware expert consensus,  
optional supervised calibration, and robustness diagnostics within a unified decision-support pipeline.

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

### How to cite

Kiratsoudis, S. (2026). *NC-ES-MADM II: Reproducible Computational Platform* (Version 1.0.0).  
Zenodo. https://doi.org/10.5281/zenodo.18318267


