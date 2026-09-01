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
| `time_to_death.csv` | one row per larva | Recorded time of death, used only for the post-mortem burden check (Fig S2). |
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

### The analysis window

`filter(Time < 37)` is applied at the **top** of the data-preparation pipeline, on the raw burden table. The 48 h sample is therefore absent from every downstream object, and the **4–36 h window and n = 86 apply to the whole paper** — the conditional-independence tests, the T50 fits, the GAMs and the SEM alike, not only to the SEM. Every analysis frame (`dat_fig3`, `dat_fig3b`, `dat_fig4d`, `dat_fig5`, `dat_cs`, `dat_fig5_mel`) carries the same 86 larvae with a maximum sampling time of 36 h; the statistics file prints this for each one.

This matters for two reasons. It is **not** a standardisation choice — standardisation merely happens after it, so describing the restriction as something the SEM does would imply the conditional-independence tests used a wider window. And any figure legend binned by time must stop at 36 h: bin labels are derived from the data rather than hard-coded, because a literal "24–48 h" label described a range the top bin never contains.

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
| `figure3.pdf` | **Fig 3** — (A) Logistic pathogen growth across the live–dead threshold (all larvae vs. survivors only). (B) Survival against burden at three sampling times, with the burden-driven null. |
| `figure4.pdf` | **Fig 4** — Activity (A) and melanization (B) against pathogen burden, coloured by time; composite health against burden at three sampling times (C); AT50 / MT50 / LT50 transition timing (D). |
| `figure5.pdf` | **Fig 5** — Bayesian SEM path effects as full posterior distributions, faceted by outcome equation, with the fitted DAG inset. The time route is shown as a marginal effect at three sampling times. |
| `figure6.pdf` | **Fig 6** — Timed ciprofloxacin intervention. Top row: grouped summaries of survival (A), burden (B), health (C). Bottom row: the continuous per-larva dose–response against injection-to-treatment delay, in the same order (D–F). |

### Supplementary figures

| File | Manuscript figure |
|------|-------------------|
| `figureS1.pdf` | **Fig S1** — Implied *m(p)* mapping with its pole at carrying capacity, three-node pathogen DAGs, and conditional-independence diagnostics. |
| `figureS2.pdf` | **Fig S2** — Post-mortem bacterial-burden stability (CFU vs. time since death). |
| `figureS3.pdf` | **Fig S3** — Four-node health DAGs and conditional-independence diagnostics. Panels D–F are the transposes of Figures 3B and 4C, read off the same two fits. |
| `figureS4.pdf` | **Fig S4** — Cumulative-burden estimation methods and the burden–time collinearity (supports Supplementary Note S2). |

Supplementary figures are numbered by float order in the manuscript, which is **not** the order the script produces them in. The filenames were previously offset from the printed numbers; they now agree, so `figureS1.pdf` is Figure S1. If the supplement is reordered, change the `ggsave` targets to match rather than renaming files by hand.

All ten figures are written directly by `main.R` and need no further editing.

### Paired panels

Several panels show the same fit more than once, and they do it in two different ways. Each pair is read off **one** model — `main.R` fits it once and reuses it, and the statistics file carries an identity check — so a coefficient repeated across two captions is deliberate. Say so in the ESM captions, or the repetition reads as a copy-paste error.

**Figures 3B and S1F are the same plot**, not transposes: both show survival against log₁₀ CFU with curves by sampling time, from the same fit (`m_cond`). What differs is the null overlaid on them, and therefore which independence is being rejected:

| Panel | Null overlaid | Rejects | Quoted |
|-------|---------------|---------|--------|
| 3B | `alive ~ log_CFU` (burden only) | *s* ⊥ *t* \| *p* | β_t = −0.22 |
| S1F | `survival ~ Time` (time only) | *s* ⊥ *p* \| *t* | β_p = −1.64 |

So it is one fit carrying two tests, not redundancy — but a reader still meets the same points and the same solid curves twice, so S1's caption needs a clause identifying F as figure 3B annotated for the reciprocal test. **S1E** is the transpose of 3B (survival against time, curves by burden); it is not the duplicate.

Note the two are not drawn identically: 3B bins time into tertiles (curves at 8/20/32 h, each clipped to its band's burden range) while S1F bins at 0–12/12–24/24–48 h (curves at 6/18/30 h, unclipped). Same fit, different display. If both stay in the paper, they are worth harmonising.

**Figures 4C and S3D are genuine transposes** of one fit (`m_D`): 4C plots health against burden at fixed times (β_p = −0.66), S3D plots health against time at fixed burdens (β_t = −0.17). The same cross-reference clause applies to the S3 caption.

Solid lines are the tested model, dashed the null — but the two panels reject their nulls differently, and the captions should not use the same phrasing. In **3B** the null has no time term, so its three band-curves coincide into one grey reference and the rejection is the *displacement* of the solid curves from it, in time order. In **4C** the null has no burden term, so each band's null is a *flat* line and the rejection is that the solid lines *slope*. A fitted line and a flat line at the same conditional mean must cross exactly once, so the gap between them is zero somewhere in every band and means nothing: do not describe 4C as "the gap between the curves".

Time is encoded as the same three sampling-time tertiles (Early 4–12 h, Mid 16–24 h, Late 28–36 h) in Figures 3B and 4A–C, with one palette defined once as `pal_tband`. Panels 4A and 4B previously used a continuous colourbar for the same variable, which put two encodings of time in one figure.

In 4C each line is drawn only where its own prediction lies inside the observable 0–7 score range. Early larvae sit at the ceiling, so the common-slope model predicts above the maximum possible score at low burden for early times; the Early solid line is consequently short (0.8 of the 7 log units on the x axis). That is the bounded-index attenuation the ESM describes, not a plotting artefact, and the censored `survreg` check is the robustness analysis for it. Curves are clipped to the range of the x variable observed within each band, because the time bands barely overlap in burden (Early 1.1–3.7, Mid 3.6–7.6, Late 4.6–8.2 log₁₀ CFU) and an unclipped curve would be mostly extrapolation. In 3B the null has no time term, so its three band-curves coincide and it is drawn once, in grey.

The bands barely share an x-range, and that is biology rather than a plotting choice: early larvae never exceeded log₁₀ CFU 3.7 and late larvae never fell below 4.6, so **no burden is simultaneously supported by all three curves**. Mid and Late do share 3.0 log units, which is a genuine matched-density comparison. For all three, the defensible reading is that each solid curve is displaced from the *same* grey null, in time order — that needs no mutual overlap, and it is what the regression actually tests, since it conditions on burden across the full range. The statistics file prints the per-band support and every pairwise overlap.

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
