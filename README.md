# Host health, not instantaneous pathogen burden, determines survival during acute *Pseudomonas aeruginosa* infection

Analysis code and data for a study of within-host infection dynamics in *Galleria mellonella* larvae infected with *Pseudomonas aeruginosa* (strain PAO1). The central finding is that host survival during acute infection is governed by **cumulative host damage** — captured through a composite health index, the temporal structure of the infection, and a timed antibiotic intervention — rather than by the **instantaneous pathogen burden** a host carries at any single moment.

The pipeline reproduces every figure and statistical result in the manuscript, from Gompertz survival modelling and logistic growth fits through causal (DAG/conditional-independence) analysis, Bayesian structural equation modelling, and the ciprofloxacin treatment-timing experiment.

## Repository structure

```
Galleria_Survival/
├── main_code.R        # Complete analysis pipeline (statistics + figures)
├── README.md
├── data/              # Raw input data (not modified by the script)
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

All inputs are comma-separated files in `data/`.

| File | Contents |
|------|----------|
| `AliveDead.csv` | Per-larva survival status across the infection time course (infected cohort). |
| `bacterial_burden.csv` | Destructively sampled bacterial burden (CFU per larva) over time. |
| `health_assesment.csv` | Per-larva activity and melanization scores; the two components of the composite health index. |
| `time_to_death.csv` | Recorded time of death for larvae that died. |
| `control_survival.csv` | Survival of uninfected / vehicle-control larvae. |
| `bacterial_burden_ab.csv` | Bacterial burden (24 h endpoint, CFU per larva) for the antibiotic-intervention experiment. |
| `Galleria_AB_Data_3rd_trial.csv` | Antibiotic-intervention data: per-larva time from injection to treatment, survival, activity, and melanization. |

**Note on the health index.** The composite health score combines activity and melanization. In the raw melanization scoring, higher values denote a *healthier* (paler) larva; the code aligns both components so that the composite is oriented with **higher = healthier** throughout.

## Requirements

- **R** ≥ 4.1
- The packages below (grouped by role, as in the script header).

```r
install.packages(c(
  # core / modelling
  "minpack.lm", "dplyr", "tidyr", "boot", "MASS", "mgcv", "purrr",
  # plotting
  "ggplot2", "cowplot", "patchwork", "scales",
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
2. Open `main_code.R` and point the working directory at the repository root. The script currently sets it explicitly near the top:

   ```r
   setwd("~/Documents/GitHub/Galleria_Survival")
   ```

   Edit this line to your local path, or comment it out and run from the repo root (e.g. via an RStudio Project).
3. Install the packages listed above.
4. Run the whole script:

   ```r
   source("main_code.R")
   ```

   Figures are written to `figures/`. The script expects that directory to exist (create it if needed: `dir.create("figures")`).

Because the SEM is sampled with `blavaan`, a full run takes a while and results are stochastic; set a seed if you need bit-for-bit reproducibility of the Bayesian estimates.

## Outputs

### Main figures

| File | Manuscript figure |
|------|-------------------|
| `figure1.pdf` | **Fig 1** — Survival curves and Gompertz mortality scenarios. |
| `figure2.pdf` | **Fig 2** — Logistic pathogen growth. |
| `figure3.pdf` | **Fig 3** — Three-node DAGs and conditional-independence diagnostics. |
| `figure4.pdf` | **Fig 4** — Activity and melanization trajectories with timing summaries. |
| `figure5.pdf` | **Fig 5** — Four-node health DAGs and diagnostics. |
| `figure6.pdf` | **Fig 6** — Bayesian SEM path effects (health-mediated). |
| `figure7.pdf` | **Fig 7** — Timed antibiotic intervention (treatment groups). |
| `figure8.pdf` | **Fig 8** — Continuous treatment-timing dose–response (per-larva). |

### Supplementary figures

| File | Manuscript figure |
|------|-------------------|
| `figureS1.pdf` | **Fig S1** — Post-mortem bacterial burden stability. |
| `figureS2.pdf` | **Fig S2** — Cumulative-burden estimation methods and the collinearity problem. |
| `figureS3.pdf` | **Fig S3** — Implied mortality mapping *m(p)*, with the pole at carrying capacity. |

### Additional outputs

| File | Description |
|------|-------------|
| `figure_overview.pdf` | Combined overview panel: survival, pathogen growth, and composite-health decline. |
| `cumulative_burden.pdf` | Cumulative pathogen-exposure inset used in the supplementary analysis. |

## Analysis overview

The script proceeds through the following stages:

1. **Data loading and health index.** Read the survival, burden, health, time-to-death, and control data; construct the activity + melanization composite health score.
2. **Survival and mortality.** Fit survival curves and Gompertz mortality models, and contrast mortality scenarios (Fig 1).
3. **Pathogen growth.** Compare logistic, exponential, and linear growth models for bacterial burden (Fig 2), check post-mortem burden stability (Fig S1), and derive the implied instantaneous-mortality mapping *m(p)* (Fig S3).
4. **Health dynamics.** Characterise activity and melanization trajectories and their timing relative to burden (Fig 4).
5. **Causal analysis.** Build candidate DAGs and run conditional-independence tests for the three-node (Fig 3) and four-node health (Fig 5) structures.
6. **Structural equation model.** Fit a Bayesian SEM of the supported causal structure, quantifying the health-mediated effects of pathogen burden and time on survival (Fig 6).
7. **Antibiotic intervention.** Analyse the ciprofloxacin treatment-timing experiment, both by treatment group (Fig 7) and as a continuous per-larva dose–response (Fig 8).
8. **Supplementary cumulative-burden analysis.** Document cumulative pathogen exposure and the collinearity between cumulative burden and sampling time that motivates the modelling choices (Fig S2).

## Citation

Karakoç C, O'Sullivan T, Gurney J, Wollein Waldetoft K, Brown SP. *Host health, not instantaneous pathogen burden, determines survival during acute Pseudomonas aeruginosa infection.* (in preparation).

## License

MIT

## Contact

Corresponding author: Sam P. Brown (sam.brown@biology.gatech.edu). 
Code: Canan Karakoç (canankarakoc@gmail.com).