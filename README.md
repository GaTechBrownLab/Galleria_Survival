# Infection history decouples pathogen burden from survival during acute *Pseudomonas aeruginosa* infection

Analysis code and data for a study of within-host infection dynamics in *Galleria mellonella* larvae infected with *Pseudomonas aeruginosa* (strain PAO1). The central finding is that host survival during acute infection is governed by **cumulative host damage** — captured through a composite health index, the temporal structure of the infection, and a timed antibiotic intervention — rather than by the **instantaneous pathogen burden** a host carries at any single moment.

The manuscript frames this as a comparison of five candidate causal structures linking time since infection (*t*), pathogen burden (*p*), host condition (*h*), and survival (*s*): burden-driven, intrinsic decline, immune collapse, instantaneous damage, and cumulative damage. Only cumulative damage survives the tests reported here.

A single script, `main_code.R`, reproduces every figure and statistical result: Gompertz survival modelling and logistic growth fits, causal (DAG / conditional-independence) analysis, Bayesian structural equation modelling of the health-mediated pathway, and the ciprofloxacin treatment-timing experiment, plus a supplementary cumulative-exposure (Σp) analysis.

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

**Note on the health index.** In the raw `health_assesment.csv` coding, melanization runs 1 (fully melanized / sickest) to 4 (no melanization / healthiest), so raw melanization is already oriented *higher = healthier*. Immediately after loading, the script reorients it to a severity scale (`melanization <- 4 - melanization`, higher = worse) because that is the orientation plotted in Figure 4B. Every composite health score is then built as `activity + (4 - melanization)`, applied to the reoriented variable, which recovers the raw value and restores a consistent **higher = healthier** orientation throughout the analysis.

Because the flip and the composite are computed inside the same `mutate()` call, the order of those expressions matters. Do not reorder them.

In the antibiotic dataset the raw `Melanization` column is used directly (`Total_health = Activity + Melanization`), which is already health-oriented and therefore consistent with `health_combined` in the main cohort.

## Requirements

- **R** ≥ 4.1 (results in the manuscript were produced under 4.5.2)
- The packages below, grouped by role.

```r
install.packages(c(
  # core / modelling
  "minpack.lm", "dplyr", "tidyr", "tibble", "boot", "MASS", "mgcv", "purrr",
  # plotting
  "ggplot2", "cowplot", "patchwork", "ggdist", "scales",
  # structural equation modelling
  "lavaan", "blavaan", "lavaanPlot", "semPlot", "brglm2", "tidybayes", "coda",
  # causal / DAGs
  "dagitty", "ggdag",
  # time handling and survival
  "lubridate", "stringr", "survival",
  # shape-constrained GAMs and misc. statistics
  "scam", "broom", "emmeans", "car", "logistf"
))
```

`logistf` supplies the penalised (Firth) likelihood used where host health near-perfectly separates survival; `coda` is used via `::` when combining MCMC chains; `car` provides `vif()` and `scales` the axis formatters.

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

The script must be run top to bottom in a clean session: several later objects (the fitted logistic parameters `K`, `p0`, `r`, and the Gompertz coefficients) are reused by downstream sections and are not recomputed.

A seed is set at each stochastic step — the growth and T50 bootstraps each seed their own resampling, and the SEM uses `set.seed(6789)` — so a clean rerun reproduces the reported estimates. The SEM sampling takes the longest by far.

## Outputs

### Main figures

| File | Manuscript figure |
|------|-------------------|
| `figure1.pdf` | **Fig 1** — The five candidate causal structures (burden-driven, intrinsic decline, immune collapse, instantaneous damage, cumulative damage). |
| `figure2.pdf` | **Fig 2** — Survival curves with Gompertz fit and instantaneous mortality *m(t)* inset, plus the growth and burden-to-mortality mapping the burden-driven reading requires. |
| `figure3.pdf` | **Fig 3** — Logistic pathogen growth across the live–dead threshold (all larvae vs. survivors only). |
| `figure4.pdf` | **Fig 4** — Activity and melanization against pathogen burden, with AT50 / MT50 / LT50 timing. |
| `figure5.pdf` | **Fig 5** — Bayesian SEM path effects for the supported structure, as full posterior distributions, with the supported DAG inset. Assembled manually from `figure5_base.pdf` and `figure5_inset.pdf` (see below). |
| `figure6.pdf` | **Fig 6** — Timed ciprofloxacin intervention: grouped summaries (panels A–C) and the continuous per-larva dose–response against injection-to-treatment delay (panels D–F). |

### Supplementary figures

| File | Manuscript figure |
|------|-------------------|
| `figureS1.pdf` | **Fig S1** — Post-mortem bacterial-burden stability (CFU vs. time since death). |
| `figureS2.pdf` | **Fig S2** — Cumulative-burden estimation methods and the burden–time collinearity (supports Supplementary Note S2). |
| `figureS3.pdf` | **Fig S3** — Implied *m(p)* mapping with its pole at carrying capacity, three-node pathogen DAGs, and conditional-independence diagnostics (supports Note S3). |
| `figureS4.pdf` | **Fig S4** — Four-node health DAGs and conditional-independence diagnostics (supports Note S3). |

The script may also emit intermediate outputs that are **not** part of the manuscript — an alternative point-estimate/interval rendering of the SEM (`p5B`), a `semPaths` diagram, and a standalone cumulative-exposure (Σp) inset — which can be ignored.

### Figure 5 is assembled manually

Every figure except Figure 5 is written directly by `main_code.R` and needs no further editing. Figure 5 is the exception: the `patchwork::inset_element` placement of the supported-DAG inset did not render acceptably, so the published figure was composed by hand from two script outputs — `figure5_base.pdf` (the posterior half-eye plot) and `figure5_inset.pdf` (the supported-model DAG).

`main_code.R` therefore does **not** write `figure5.pdf`. The published composite is committed to `figures/` and is not regenerated by a rerun, so a full rerun will not overwrite it.

## Analysis overview

The script proceeds through the following stages:

1. **Candidate causal structures.** Draw the five hypotheses as DAGs (Fig 1). Schematic only; no statistics.
2. **Data loading and health index.** Read the survival, burden, health, time-to-death, and control data; construct the activity + melanization composite health score (higher = healthier).
3. **Survival and mortality.** Fit survival curves and Gompertz mortality models, and show what a strictly burden-driven account would require (Fig 2).
4. **Pathogen growth.** Compare logistic, exponential, and linear growth models for bacterial burden — fit to all larvae and to survivors only to expose survivor bias (Fig 3) — check post-mortem burden stability (Fig S1), and derive the implied instantaneous-mortality mapping *m(p)* (a panel of Fig S3).
5. **Health dynamics.** Characterise activity and melanization trajectories and their timing (AT50, MT50, LT50) relative to burden (Fig 4).
6. **Causal analysis.** Build candidate DAGs and run conditional-independence tests for the three-node pathogen (Fig S3) and four-node health (Fig S4) structures, including the melanization-only sensitivity analysis.
7. **Structural equation model.** Fit a Bayesian SEM of the supported causal structure (*t → p → h → s* with *t → h*) and summarise the health-mediated effects of pathogen burden and time on survival as posterior distributions (Fig 5).
8. **Antibiotic intervention.** Analyse the ciprofloxacin treatment-timing experiment, both by treatment group and as a continuous per-larva dose–response against injection-to-treatment delay (Fig 6).
9. **Supplementary cumulative-burden analysis.** Re-express pre-treatment exposure as the integral of fitted logistic growth (Σp), test whether it adds information beyond time, and document the cumulative-burden / sampling-time collinearity (Fig S2; Supplementary Note S2).

Note that the script's internal section order does not match the manuscript figure order: Figure 5 (SEM) is produced after Figure 6 (antibiotics), and the supplementary figures are interleaved. The `ggsave` filenames are authoritative.

## Citation

Karakoç C, O'Sullivan T, Gurney J, Martigoni M, Wollein Waldetoft K, Brown SP. *Infection history decouples pathogen burden from survival during acute Pseudomonas aeruginosa infection.* (in preparation).

## License

MIT

## Contact

Corresponding author: Sam P. Brown (sam.brown@biology.gatech.edu).
Code: Canan Karakoç (canankarakoc@gmail.com).