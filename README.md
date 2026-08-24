# Infection history decouples pathogen burden from survival during acute *Pseudomonas aeruginosa* infection

Analysis code and data for a study of within-host infection dynamics in *Galleria mellonella* larvae infected with *Pseudomonas aeruginosa* (strain PAO1). The central finding is that host survival during acute infection is governed by **cumulative host damage** — captured through a composite health index, the temporal structure of the infection, and a timed antibiotic intervention — rather than by the **instantaneous pathogen burden** a host carries at any single moment.

The manuscript frames this as a comparison of five candidate causal structures linking time since infection (*t*), pathogen burden (*p*), host condition (*h*), and survival (*s*): burden-driven, intrinsic decline, immune collapse, instantaneous damage, and cumulative damage. Only cumulative damage survives the tests reported here.

A single script, `main.R`, reproduces every figure and statistical result: Gompertz survival modelling and logistic growth fits, causal (DAG / conditional-independence) analysis, Bayesian structural equation modelling of the health-mediated pathway, and the ciprofloxacin treatment-timing experiment, plus a supplementary cumulative-exposure (Σp) analysis. It also writes a single labelled statistics file that every number quoted in the manuscript can be checked against.

## Repository structure

```
Galleria_Survival/
├── main.R             # Complete analysis pipeline (statistics + figures)
├── README.md
├── data/              # Raw input data (read-only; not modified by the script)
│   ├── survival.csv
│   ├── bacterial_burden.csv
│   ├── health_assesment.csv
│   ├── time_to_death.csv
│   ├── control_survival.csv
│   ├── bacterial_burden_ab.csv
│   └── health_assesment_ab.csv
├── figures/           # Output figures (created when the script runs)
└── results/           # statistics_for_manuscript.txt (created when the script runs)
```

The working directory also contains an `old/` folder of superseded script versions. It is listed in `.gitignore` and is not part of the repository.

## Data

All inputs are comma-separated files in `data/`. The first five drive the infection-dynamics analysis; the last two are the antibiotic-intervention experiment. Note that the two `*survival*.csv` files are **cohort counts per timepoint**, while everything else is **one row per larva**.

| File | Grain | Contents |
|------|-------|----------|
| `survival.csv` | one row per hour | Infected cohort: `time`, `total`, `dead`, `alive`. Cohort counts, not per-larva records. Drives the Gompertz fit (Fig 2). |
| `control_survival.csv` | one row per hour | Same four columns for the uninfected / vehicle-control cohort. |
| `bacterial_burden.csv` | one row per plate | Replicate colony counts (`rep1`–`rep4`) for destructively sampled larvae over the time course, with the `Dilution`, `Volume_ul` and `Buffer_dilution_factor` needed to convert them to CFU per larva. `Countable` marks unreadable plates. Carries **no** live/dead status — that is joined in from `health_assesment.csv`. |
| `health_assesment.csv` | one row per larva | `activity`, `melanization` and `survival` for the time-course cohort. The first two are the components of the composite health index; `survival` supplies the live/dead status used throughout. |
| `time_to_death.csv` | one row per larva | Recorded time of death, used only for the post-mortem burden check (Fig S1). |
| `bacterial_burden_ab.csv` | one row per plate | Replicate colony counts (`rep1`–`rep5`) at the 24 h endpoint of the antibiotic experiment, with the same dilution columns. |
| `health_assesment_ab.csv` | one row per larva | Antibiotic experiment: `Treatment`, the injection / treatment / sampling clock times, and `Activity`, `Melanization`, `Survival`. The injection-to-treatment delay used in Fig 6D–F is derived from the clock times. |

Some columns are present but unused by the analysis: `Contamination` and `Batch` in `bacterial_burden.csv`, `cocoon_formation` and `total_score` in `health_assesment.csv`, and `Coccoon` in `health_assesment_ab.csv`.

### The health index

Both component scores are used **exactly as recorded, with no reorientation**:

| Score | Range | Direction |
|-------|-------|-----------|
| `activity` | 0 (no movement) – 3 (fully active) | higher = healthier |
| `melanization` | 0 (fully melanized) – 4 (no melanization) | higher = healthier |

Both are already oriented *higher = healthier*, so both decline as infection progresses, and the composite is simply their sum:

```r
health_combined = activity + melanization      # range 0–7, higher = healthier
```

The same convention applies in the antibiotic dataset (`Total_health = Activity + Melanization`), so the two cohorts are directly comparable. Both shape-constrained smooths in Figure 4 are therefore fitted as monotone **decreasing** in burden (`bs = "mpd"`).

Earlier versions of this code flipped melanization to a severity scale and rebuilt the composite as `activity + (4 - melanization)`. That is no longer done anywhere, and no `4 - melanization` term should appear in the script.

### Detection limit

A larva plated with no colonies is a measurement at the limit of detection, not a missing value. The lowest-dilution plates were 100 µl from a 250 µl homogenate at 10×, so one colony corresponds to

```
1 colony × 10 (dilution) × (250 µl / 100 µl) = 25 CFU per larva
```

Larvae yielding no colonies on any plate therefore enter the analysis at **LOD/2 = 12.5 CFU** rather than being dropped. Plates with no readable count (`Countable == "no"`) remain excluded. This gives **n = 86** larvae in the main cohort over the 4–36 h window (49 alive, 37 dead at sampling), one of which sits at LOD/2.

In the antibiotic cohort the same rule applies per larva. One larva (G4 L20) had a dedicated 100 µl plate reading zero alongside a 10 µl spot on a shared plate reading 2 colonies; because multi-sample spot plates are prone to carryover, that spot is excluded and the larva enters at LOD/2.

All standardised variables (`scaled_time`, `scaled_cfu`, `scaled_health`, …) are computed **after** the `Time < 37` filter, so they are standardised on the sample the models actually use.

## Requirements

- **R** ≥ 4.1 (results in the manuscript were produced under 4.5.2)
- The packages below, grouped by role.

```r
install.packages(c(
  # core / modelling
  "minpack.lm", "dplyr", "tidyr", "tibble", "boot", "MASS", "mgcv", "purrr",
  # plotting
  "ggplot2", "patchwork", "ggdist", "scales",
  # structural equation modelling
  "lavaan", "blavaan", "semPlot", "tidybayes", "coda",
  # causal / DAGs
  "dagitty", "ggdag",
  # time handling and survival
  "lubridate", "stringr", "survival",
  # shape-constrained GAMs and misc. statistics
  "scam", "broom", "car", "logistf"
))
```

`logistf` supplies the penalised (Firth) likelihood used where host health near-perfectly separates survival; `survival` provides `survreg()` for the censored-health robustness check; `car` provides `vif()` and `scales` the axis formatters. `coda` and `tibble` are needed but never attached: `coda` is called as `coda::is.mcmc()`, and `tibble()` reaches the script through dplyr's re-export. `MASS` **is** attached, and it masks `dplyr::select()` — which is why every `select()` in the script is written as `dplyr::select()`. If you drop `library(MASS)` (only `MASS::dose.p()` is used), those prefixes become optional but must not be removed carelessly.

**Bayesian SEM backend.** `blavaan` fits models via MCMC and requires a sampling backend — either **Stan** (`rstan`) or **JAGS** (`rjags`, which needs a system JAGS install). Set this up following the `blavaan` documentation before running the SEM section; this is the most computationally intensive step.

## Running the analysis

1. Clone the repository and place the seven data files in `data/`.
2. Run from the repository root — most simply by opening `Galleria_Survival.Rproj` in RStudio. The script does not set a working directory of its own; it checks that `data/` is visible and stops with a clear message if it is not. If you need to point it somewhere else, uncomment and edit the line near the top:

   ```r
   # setwd("~/Documents/GitHub/Galleria_Survival")
   ```
3. Install the packages listed above and configure the `blavaan` backend.
4. Run the script:

   ```r
   source("main.R")
   ```

   `figures/` and `results/` are created automatically. A full run takes a few minutes, dominated by the three SEM fits.

**Run it top to bottom in a clean session.** This matters more than it sounds. Later sections reuse fitted objects from earlier ones — the logistic parameters `K`, `p0`, `r`, the Gompertz coefficients, and the SEM posterior draws — and none of them are recomputed. Running sections out of order, or into a session left over from a previous version of the script, can silently produce figures built from stale objects rather than failing. Restarting R (or running with `Rscript`, which does not restore `.RData`) avoids this entirely.

A seed is set at each stochastic step — the growth and T50 bootstraps each seed their own resampling, the plotted jitter is seeded, and the SEMs use `set.seed(6789)` and `set.seed(6790)` — so a clean rerun reproduces the reported estimates.

## Outputs

### Main figures

| File | Manuscript figure |
|------|-------------------|
| `figure1.pdf` | **Fig 1** — The five candidate causal structures (burden-driven, intrinsic decline, immune collapse, instantaneous damage, cumulative damage). |
| `figure2.pdf` | **Fig 2** — Survival curves with Gompertz fit and instantaneous mortality *m(t)* inset, plus the growth and burden-to-mortality mapping the burden-driven reading requires. |
| `figure3.pdf` | **Fig 3** — Logistic pathogen growth across the live–dead threshold (all larvae vs. survivors only). |
| `figure4.pdf` | **Fig 4** — Activity (A) and melanization (B) against pathogen burden, coloured by time, with the AT50 / MT50 / LT50 transition timing (C). |
| `figure5.pdf` | **Fig 5** — Bayesian SEM path effects as full posterior distributions, faceted by outcome equation, with the fitted DAG inset. The time route is shown as a marginal effect at three sampling times. |
| `figure6.pdf` | **Fig 6** — Timed ciprofloxacin intervention. Top row: grouped summaries of survival (A), burden (B), health (C). Bottom row: the continuous per-larva dose–response against injection-to-treatment delay, in the same order (D–F). |

### Supplementary figures

| File | Manuscript figure |
|------|-------------------|
| `figureS1.pdf` | **Fig S1** — Post-mortem bacterial-burden stability (CFU vs. time since death). |
| `figureS2.pdf` | **Fig S2** — Cumulative-burden estimation methods and the burden–time collinearity (supports Supplementary Note S2). |
| `figureS3.pdf` | **Fig S3** — Implied *m(p)* mapping with its pole at carrying capacity, three-node pathogen DAGs, and conditional-independence diagnostics (supports Note S3). |
| `figureS4.pdf` | **Fig S4** — Four-node health DAGs and conditional-independence diagnostics (supports Note S3). |

All ten figures are written directly by `main.R` and need no further editing.

### DAG rendering

Every DAG in the paper is drawn by one function (`draw_dag()`) with one node layout, one geometry, and one type scale, so panels are directly comparable across Figures 1, S3, S4 and the Figure 5 inset. Structural DAGs are drawn in **grey**: they are possibilities under test, not results. Colour appears in one place only — the Figure 5 inset — where each edge carries the sign of its fitted path (teal arrowhead for positive, red blunt bar for negative), read from the SEM posterior.

### Statistics file

`results/statistics_for_manuscript.txt` collects every value quoted in the manuscript in one place, grouped by where it appears, and compares each against the value currently in the text:

```
  Gompertz b (h^-1)                            0.4                     [ok]
  standardised beta (z_time:z_cfu)             -0.939                  [ok]
    text claims p < 0.001                      4.45e-08                [ok: < 0.001]
  MT50 as printed in the fig 4 caption         21                      [ok]
```

| Marker | Meaning |
|--------|---------|
| `[ok]` | recomputed value matches the manuscript |
| `[CHECK: ms = x]` | they differ — reconcile before submission |
| `[ok: < x]` | a claim stated as a bound, and the bound holds |
| `MISSING` | object not in the workspace (section not run?) |

The reference column is meant to mirror **what the manuscript currently says**, not what the code computes. Nine of those values were updated after the detection-limit correction; the file header lists them with their previous values, and each is annotated `# reconciled 2026-08-23 (was …)` at its call site, so an `[ok]` can always be distinguished from a number silently compared with itself.

## Analysis overview

The script proceeds through the following stages:

1. **Candidate causal structures.** Draw the five hypotheses as DAGs (Fig 1). Schematic only; no statistics.
2. **Data loading and health index.** Read the survival, burden, health, time-to-death, and control data; apply the detection limit; construct the activity + melanization composite health score (higher = healthier).
3. **Survival and mortality.** Fit survival curves and Gompertz mortality models, and show what a strictly burden-driven account would require (Fig 2).
4. **Pathogen growth.** Compare logistic, exponential, and linear growth models for bacterial burden — fit to all larvae and to survivors only to expose survivor bias (Fig 3) — check post-mortem burden stability (Fig S1), and derive the implied instantaneous-mortality mapping *m(p)* (a panel of Fig S3).
5. **Health dynamics.** Characterise activity and melanization trajectories and their timing (AT50, MT50, LT50) relative to burden (Fig 4).
6. **Causal analysis.** Build candidate DAGs and run conditional-independence tests for the three-node pathogen (Fig S3) and four-node health (Fig S4) structures, including the melanization-only sensitivity analysis and a censored-regression robustness check for the 0–7 floor and ceiling of the health score.
7. **Structural equation model.** Fit Bayesian SEMs of the supported causal structure (*t → p → h → s* with *t → h*) in a linear and a quadratic-in-time form, and summarise the health-mediated effects as posterior distributions (Fig 5). Because the quadratic model makes the *t → h* path depend on when the larva was sampled, that path is reported as a marginal effect at three sampling times rather than as a single coefficient.
8. **Antibiotic intervention.** Analyse the ciprofloxacin treatment-timing experiment, both by treatment group and as a continuous per-larva dose–response against injection-to-treatment delay (Fig 6).
9. **Supplementary cumulative-burden analysis.** Re-express pre-treatment exposure as the integral of fitted logistic growth (Σp), test whether it adds information beyond time, and document the cumulative-burden / sampling-time collinearity (Fig S2; Supplementary Note S2).
10. **Manuscript number report.** Recompute every quoted value and write `results/statistics_for_manuscript.txt`.

A note on model comparison in stage 7: the DIC of the no-mediation model is **not** comparable with the two mediation models, because it has no health outcome and so sums deviance over two modelled variables rather than three. Only the linear and quadratic mediation models can be compared on DIC; the no-mediation model is judged on its posterior predictive *p* alone. The statistics file flags this where the numbers are printed.

Note that the script's internal section order does not match the manuscript figure order: Figure 5 (SEM) is produced after Figure S4, and the supplementary figures are interleaved. The `ggsave` filenames are authoritative.

## Citation

Karakoç C, O'Sullivan T, Gurney J, Martigoni M, Wollein Waldetoft K, Brown SP. *Infection history decouples pathogen burden from survival during acute Pseudomonas aeruginosa infection.* (in preparation).

## License

MIT

## Contact

Corresponding author: Sam P. Brown (sam.brown@biology.gatech.edu).
Code: Canan Karakoç (canankarakoc@gmail.com).
