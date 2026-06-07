# Host health, not instantaneous pathogen burden, determines survival during acute *Pseudomonas aeruginosa* infection

Analysis code and data for a study of within-host infection dynamics in *Galleria mellonella* larvae infected with *Pseudomonas aeruginosa* (strain PAO1). The central finding is that host survival during acute infection is governed by **cumulative host damage** — captured through a composite health index, the temporal structure of the infection, and a timed antibiotic intervention — rather than by the **instantaneous pathogen burden** a host carries at any single moment.

A single script, `main_code.R`, reproduces every figure and statistical result in the manuscript: Gompertz survival modelling and logistic growth fits, causal (DAG / conditional-independence) analysis, Bayesian structural equation modelling of the health-mediated pathway, and the ciprofloxacin treatment-timing experiment, plus a supplementary cumulative-exposure (Σp) analysis.

## Repository structure

```
Galleria_Survival/
├── main_code.R        # Complete analysis pipeline (statistics + figures)
├── README.md
├── data/              # Raw input data (read-only; not modified by the script)
│   ├── AliveDead.csv
│   ├── bacterial_burden.csv
│   ├── health_assesment.csv
│   ├── time_to_death.csv
│   ├── control_survival.csv
│   ├── bacterial_burden_ab.csv
│   └── Galleria_AB_Data_3rd_trial.csv
└── figures/           # Output figures (created when the script runs)
```

## Data

All inputs are comma-separated files in `data/`. The first five drive the infection-dynamics analysis; the last two are the antibiotic-intervention experiment.

| File | Contents |
|------|----------|
| `AliveDead.csv` | Per-larva survival status across the infection time course (infected cohort). |
| `bacterial_burden.csv` | Destructively sampled bacterial burden (CFU per larva) over time, with live/dead status. |
| `health_assesment.csv` | Per-larva activity and melanization scores; the two components of the composite health index. |
| `time_to_death.csv` | Recorded time of death for larvae that died. |
| `control_survival.csv` | Survival of uninfected / vehicle-control larvae. |
| `bacterial_burden_ab.csv` | Bacterial burden (24 h endpoint, CFU per larva) for the antibiotic-intervention experiment. |
| `Galleria_AB_Data_3rd_trial.csv` | Antibiotic-intervention data: per-larva time from injection to treatment, survival, activity, and melanization. |

**Note on the health index.** In the raw `health_assesment.csv` coding, melanization runs 1 (fully melanized / sickest) to 4 (pale / healthiest), so raw melanization is already oriented *higher = healthier*. Internally the script briefly reorients melanization to a severity scale (higher = worse) for one trajectory plot, but every composite health score is built as `activity + (4 − melanization)`, which restores a consistent **higher = healthier** orientation throughout the analysis and figures.

## Requirements

- **R** ≥ 4.1
- The packages below (grouped by role, mirroring the script header). `ggdist` is required for the posterior path-effect figure; `scales` and `car` are used via `::` (for axis formatting and `vif()`).

```r
install.packages(c(
  # core / modelling
  "minpack.lm", "dplyr", "tidyr", "boot", "MASS", "mgcv", "purrr",
  # plotting
  "ggplot2", "cowplot", "patchwork", "ggdist", "scales",
  # structural equation modelling
  "lavaan", "blavaan", "lavaanPlot", "semPlot", "brglm2", "tidybayes",
  # causal / DAGs
  "dagitty", "ggdag",
  # time handling and survival
  "lubridate", "stringr", "survival",
  # shape-constrained GAMs and misc. statistics
  "scam", "broom", "emmeans", "car"
))
```

**Bayesian SEM backend.** `blavaan` fits models via MCMC and requires a sampling backend — either **Stan** (`rstan`) or **JAGS** (`rjags`, which needs a system JAGS install). Set this up following the `blavaan` documentation before running the SEM section; this is the most computationally intensive step.

## Running the analysis

1. Clone the repository and place the seven data files in `data/`.
2. Open `main_code.R` and point the working directory at the repository root. The script sets it explicitly near the top:

   ```r
   setwd("~/Documents/GitHub/Galleria_Survival")
   ```

   Edit this line to your local path, or comment it out and run from the repo root (e.g. via an RStudio Project).
3. Install the packages listed above and configure the `blavaan` backend.
4. Create the output directory if it does not exist, then run the script:

   ```r
   dir.create("figures", showWarnings = FALSE)
   source("main_code.R")
   ```

   Figures are written to `figures/`.

The script sets seeds at each stochastic step (the growth/T50 bootstraps and the SEM, `set.seed(6789)`), so a clean rerun reproduces the reported estimates. The SEM sampling still takes the longest by far.

## Outputs

### Main figures

| File | Manuscript figure |
|------|-------------------|
| `figure1.pdf` | **Fig 1** — Survival curves and Gompertz mortality scenarios. |
| `figure2.pdf` | **Fig 2** — Logistic pathogen growth across the live–dead threshold. |
| `figure3.pdf` | **Fig 3** — Implied *m(p)* mapping (pole at carrying capacity), three-node DAGs, and conditional-independence diagnostics. |
| `figure4.pdf` | **Fig 4** — Activity and melanization trajectories with T50 timing. |
| `figure5.pdf` | **Fig 5** — Four-node health DAGs and diagnostics. |
| `figure6_posterior.pdf` / `figure6_posterior.png` | **Fig 6** — Bayesian SEM path effects, shown as full posterior distributions (half-eye). |
| `figure7.pdf` | **Fig 7** — Timed ciprofloxacin intervention: grouped summaries (panels A–C) and the continuous per-larva dose–response (panels D–F). |

### Supplementary and supporting outputs

| File | Description |
|------|-------------|
| `figureS1.pdf` | **Fig S1** — Post-mortem bacterial-burden stability (CFU vs. time since death). |
| `figureS3.pdf` | Cumulative-burden estimation methods and the burden–time collinearity that motivates the modelling choices (supports Supplementary Note 2). |
| `cumulative_burden.pdf` | Cumulative pathogen-exposure (Σp) inset used in the damage analysis. |

**Notes on file names.** A few output filenames lag the manuscript numbering and are easy to misread:

- `figure6.pdf` is also written — it is the earlier point-estimate / interval ("bar") rendering of the same SEM and has been **superseded by `figure6_posterior.pdf`** as manuscript Figure 6.
- `figure7.pdf` is manuscript Figure 7, despite an inline `# >>> Manuscript Figure 8` comment in the script (stale label).
- `figureS3.pdf` is the cumulative-burden supplement; the current manuscript cites only Supplementary Figure S1 and folds the cumulative-exposure material into Supplementary Note 2.

## Analysis overview

The script proceeds through the following stages:

1. **Data loading and health index.** Read the survival, burden, health, time-to-death, and control data; construct the activity + melanization composite health score (higher = healthier).
2. **Survival and mortality.** Fit survival curves and Gompertz mortality models, and contrast pathogen- vs. damage-driven mortality scenarios (Fig 1).
3. **Pathogen growth.** Compare logistic, exponential, and linear growth models for bacterial burden — fit to all larvae and to survivors only to expose survivor bias (Fig 2) — check post-mortem burden stability (Fig S1), and derive the implied instantaneous-mortality mapping *m(p)* (panel of Fig 3).
4. **Health dynamics.** Characterise activity and melanization trajectories and their timing (AT50, MT50, LT50) relative to burden (Fig 4).
5. **Causal analysis.** Build candidate DAGs and run conditional-independence tests for the three-node (Fig 3) and four-node health (Fig 5) structures.
6. **Structural equation model.** Fit a Bayesian SEM of the supported causal structure (*t → p → h → s* with *t → h*) and summarise the health-mediated effects of pathogen burden and time on survival as posterior distributions (Fig 6).
7. **Antibiotic intervention.** Analyse the ciprofloxacin treatment-timing experiment, both by treatment group and as a continuous per-larva dose–response against injection-to-treatment delay (Fig 7).
8. **Supplementary cumulative-burden analysis.** Re-express pre-treatment exposure as the integral of fitted logistic growth (Σp), test whether it adds information beyond time, and document the cumulative-burden / sampling-time collinearity (figS3, `cumulative_burden.pdf`; Supplementary Note 2).

## Citation

Karakoç C, O'Sullivan T, Gurney J, Martigoni M, Wollein Waldetoft K, Brown SP. *Host health, not instantaneous pathogen burden, determines survival during acute Pseudomonas aeruginosa infection.* (in preparation).

## License

MIT

## Contact

Corresponding author: Sam P. Brown (sam.brown@biology.gatech.edu).
Code: Canan Karakoç (canankarakoc@gmail.com).