# Reproducibility protocol

This document describes how to reproduce the computational results reported in the accompanying manuscript using the NC-ES-MADM II platform.

## Software requirements
- R (version ≥ 4.2)
- Packages: ggplot2, dplyr, tidyr, patchwork, scales

## Reproducing the results
1. Run the core model script located in `code/` to generate decision weights, posterior distributions, and diagnostic indices.
2. Use the input datasets provided in `data/raw/`.
3. Verify numerical outputs against the tables provided in `results/`.

## Reproducing the figures
- Figures 4.1–4.6: execute the corresponding script in `plots/`
- Figures 4.7–4.13: execute the corresponding script in `plots/`
- Figure 5.1: execute the corresponding script in `plots/`

All plotting scripts are deterministic and do not perform file input/output operations.

## Notes
This platform prioritizes transparency, interpretability, and auditability over algorithmic complexity.
