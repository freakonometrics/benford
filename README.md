# What Does a Benford Test Actually Test?

Reproducible R/Quarto companion for the paper
**“What Does a Benford Test Actually Test? Marginal Conformity, Sampling Structure, and Forensic Inference.”**

The numerical study separates three questions that are often conflated in Benford applications:

1. whether the logarithmic significand is marginally Benford;
2. whether a reported p-value is calibrated for the actual sampling structure;
3. what the observed digit information can identify about the generating mechanism.

The workflow contains the four controlled Benford constructions, design-aware Pearson calibration,
a continuous-significand CvM robustness check, latent-mixture experiments, power simulations, and the multiplication-game examples.

## Reproduce locally

Required R packages are `ggplot2`, `dplyr`, `tidyr`, and `knitr`.

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "knitr"))
source("test.R")
```

Then render the complete analysis:

```bash
quarto render benford-simulations.qmd --to html --output index.html
```

The first fresh render recomputes any missing cached simulation objects in `output/data/`.
All Monte Carlo seeds are fixed and printed in the rendered document.

## Repository structure

- `benford-simulations.qmd`: master reproducible analysis
- `R/`: generators, statistics, simulation engines, theory, design-aware calibration, and output helpers
- `test.R`: smoke tests
- `output/`: generated data, figures, and LaTeX tables (not versioned by default)
- `.github/workflows/pages.yml`: manual GitHub Pages build/deployment
