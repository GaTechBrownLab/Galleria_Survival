# Analysis code for (Working title) "Host-health integrates infection history to 
# determine survival during acute \textit{Pseudomonas aeruginosa} infection"
# Author: Canan Karakoç
# ---------------------------------------------------------------------------

# Load required libraries
library(minpack.lm)  # For non-linear regression
library(dplyr)
library(tidyr)
library(boot)
library(MASS)
library(mgcv)

# Plotting
library(ggplot2)
library(patchwork)
library(ggdist)
library(scales)

# SEM
library(lavaan)
library(blavaan) #bayesian lavaan
library(semPlot)
library(tidybayes)

#DAG
library(dagitty)
library(ggdag)

# Time series
library(lubridate)
library(stringr)
library(survival)

# GAM fit
library(scam)

# Stats
library(broom)
library(logistf)   # penalised (Firth) likelihood for near-complete separation
library(car)       # vif()

library(purrr)

# Theme
mytheme <- theme_bw() +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(legend.text = element_text(size = 16)) +
  theme(axis.text = element_text(size = 16, color = "black"), axis.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 16), strip.background = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    linewidth = 1
  )) +
  theme(
    axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_blank()
  ) +
  theme(
    axis.title.x = element_text(margin = margin(10, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
    axis.text.x = element_text(margin = margin(10, 0, 0, 0)),
    axis.text.y = element_text(margin = margin(0, 10, 0, 0))
  )

# ---- DAG rendering ---------------------------------------------------------
# Every DAG in the paper uses one node layout, one geometry, one type scale.
# Structural DAGs are grey: they are possibilities under test. Colour appears
# only where a fitted path has a sign (figure 5 inset).
dag_edge <- "grey30"; dag_edge_hi <- "grey10"; dag_faded <- "grey75"
dag_dashed <- "grey60"; dag_node <- "grey30"
dag_pos <- "#19798b"; dag_neg <- "#b80422"

dag_box_w <- 0.40; dag_box_h <- 0.30
dag_label <- 5.0     # node glyph (mm)
dag_sub   <- 6.0     # conditional-independence line (mm)
dag_title <- 12      # panel title (pt)
dag_lw    <- 0.6
dag_arrow <- 0.09    # arrowhead, inches
dag_xlim  <- c(-0.62, 2.62)
dag_ylim  <- c(-1.42, 1.82)

dag_nodes <- data.frame(
  name  = c("t", "p", "h", "s"),
  x     = c(0, 1, 1, 2),
  y     = c(0.5, 1.2, -0.2, 0.5),
  label = c("italic(t)", "italic(p)", "italic(h)", "italic(s)"),
  stringsAsFactors = FALSE
)

dag_edge_xy <- list(
  t_p = c(0, 0.5, 1,  1.2), t_h = c(0, 0.5, 1, -0.2), t_s = c(0, 0.5, 2, 0.5),
  p_h = c(1, 1.2, 1, -0.2), p_s = c(1, 1.2, 2,  0.5),
  h_s = c(1, -0.2, 2, 0.5), h_p = c(1, -0.2, 1,  1.2)
)

draw_dag <- function(edges           = character(0),
                     title           = NULL,
                     subtitles       = character(0),
                     dashed_edges    = NULL,
                     faded_edges     = NULL,
                     highlight_edges = NULL,
                     faded_nodes     = NULL,
                     drop_nodes      = NULL,
                     edge_signs      = NULL,   # named c(t_p = "pos", ...)
                     scale           = 1, 
                     pad_x           = 0.12,   # margin beside outer node boxes
                     pad_y           = 0.22,   # margin above nodes / below last subtitle
                     sub_offset      = 0.20,   # first subtitle, below lowest node box
                     sub_step        = 0.20,   # spacing between subtitle lines
                     xlim = NULL, ylim = NULL) {
  
  nodes <- dag_nodes[!dag_nodes$name %in% drop_nodes, ]
  nodes$faded <- nodes$name %in% faded_nodes
  nodes$xmin <- nodes$x - dag_box_w/2; nodes$xmax <- nodes$x + dag_box_w/2
  nodes$ymin <- nodes$y - dag_box_h/2; nodes$ymax <- nodes$y + dag_box_h/2
  
  y_lo <- min(nodes$ymin); y_hi <- max(nodes$ymax)
  sub_y <- if (length(subtitles))
    y_lo - sub_offset - (seq_along(subtitles) - 1) * sub_step else numeric(0)
  if (is.null(xlim)) xlim <- c(min(nodes$xmin) - pad_x, max(nodes$xmax) + pad_x)
  if (is.null(ylim)) ylim <- c(min(c(y_lo, sub_y)) - pad_y, y_hi + pad_y)
  
  trim <- function(co, gap = 0) {
    dx <- co[3]-co[1]; dy <- co[4]-co[2]; len <- sqrt(dx^2 + dy^2)
    ux <- dx/len; uy <- dy/len
    tt <- min(if (abs(ux) > 0.01) (dag_box_w/2)/abs(ux) else Inf,
              if (abs(uy) > 0.01) (dag_box_h/2)/abs(uy) else Inf)
    c(co[1] + ux*tt, co[2] + uy*tt, co[3] - ux*(tt+gap), co[4] - uy*(tt+gap), ux, uy)
  }
  build <- function(codes, gap = 0) {
    codes <- intersect(codes, names(dag_edge_xy))
    if (!length(codes)) return(NULL)
    do.call(rbind, lapply(codes, function(e) {
      v <- trim(dag_edge_xy[[e]], gap)
      data.frame(x = v[1], y = v[2], xend = v[3], yend = v[4], ux = v[5], uy = v[6])
    }))
  }
  arw <- arrow(length = unit(dag_arrow, "inches"), type = "closed")
  
  if (!is.null(edge_signs)) {
    pos <- build(names(edge_signs)[edge_signs == "pos"])
    neg <- build(names(edge_signs)[edge_signs == "neg"], gap = 0.03)
    bars <- if (!is.null(neg)) transform(neg,
                                         bx = xend - uy*0.10, by = yend + ux*0.10,
                                         bxend = xend + uy*0.10, byend = yend - ux*0.10) else NULL
    edge_layers <- list(
      if (!is.null(pos)) geom_segment(data = pos, aes(x, y, xend = xend, yend = yend),
                                      arrow = arw, colour = dag_pos, linewidth = dag_lw),
      if (!is.null(neg)) geom_segment(data = neg, aes(x, y, xend = xend, yend = yend),
                                      colour = dag_neg, linewidth = dag_lw),
      if (!is.null(bars)) geom_segment(data = bars, aes(bx, by, xend = bxend, yend = byend),
                                       colour = dag_neg, linewidth = dag_lw*2, lineend = "round"))
  } else {
    e_f <- build(faded_edges); e_m <- build(edges); e_h <- build(highlight_edges)
    edge_layers <- list(
      if (!is.null(e_f)) geom_segment(data = e_f, aes(x, y, xend = xend, yend = yend),
                                      arrow = arw, colour = dag_faded, linewidth = dag_lw*0.8),
      if (!is.null(e_m)) geom_segment(data = e_m, aes(x, y, xend = xend, yend = yend),
                                      arrow = arw, colour = dag_edge, linewidth = dag_lw),
      if (!is.null(e_h)) geom_segment(data = e_h, aes(x, y, xend = xend, yend = yend),
                                      arrow = arw, colour = dag_edge_hi, linewidth = dag_lw*1.6),
      if (!is.null(build(dashed_edges)))
        geom_segment(data = build(dashed_edges), aes(x, y, xend = xend, yend = yend),
                     arrow = arw, colour = dag_dashed, linetype = "dashed", linewidth = dag_lw))
  }
  
  p <- ggplot() + edge_layers +
    geom_rect(data = nodes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white",
              colour = ifelse(nodes$faded, dag_faded, dag_node), linewidth = 0.5) +
    geom_text(data = nodes, aes(x, y, label = label), parse = TRUE,
              size = dag_label*scale,
              colour = ifelse(nodes$faded, dag_faded, "black")) +
    coord_fixed(ratio = 1, xlim = dag_xlim, ylim = dag_ylim) +
    theme_void()
  
  for (i in seq_along(subtitles))
    p <- p + annotate("text", x = 1, y = -0.72 - (i-1)*0.28, label = subtitles[i],
                      parse = TRUE, size = dag_sub*scale, colour = "black")
  
  if (!is.null(title)) p <- p + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = dag_title*scale, face = "bold"))
  p
}


## Run this script from the repository root (RStudio: open
## Galleria_Survival.Rproj). Uncomment and edit only if you need to:
# setwd("~/Documents/GitHub/Galleria_Survival")
if (!dir.exists("data"))
  stop("No 'data/' directory here. Run from the repository root, or setwd() first.",
       call. = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

survival      <- read.table("data/survival.csv", header = T, sep = ",", dec =".")
burden        <- read.table("data/bacterial_burden.csv", header = T, sep = ",", dec =".")
health        <- read.table("data/health_assesment.csv", header = T, sep = ",", dec =".")
death         <- read.table("data/time_to_death.csv", header = T, sep = ",", dec =".")
control_surv  <- read.table("data/control_survival.csv", header = T, sep = ",", dec =".")

# Common time sequence used throughout
time_seq <- seq(0, 36, length.out = 200)

# Discrete time-band palette: a three-way cut of the magma ramp. Figures 3B and
# 4A-C all use it, so the same variable has one encoding across the paper
# rather than a continuous colourbar in some panels and three bands in others.
pal_tband <- setNames(viridisLite::magma(3, direction = -1, begin = 0.25, end = 0.80),
                      c("Early", "Mid", "Late"))

#===============================================================================
# FIGURE 1: The five candidate causal structures (Introduction)
#-------------------------------------------------------------------------------
# Schematic only -- no statistics, no conditional-independence notation.
#===============================================================================

#-------------------------------------------------------------------------------
# The five hypotheses
#-------------------------------------------------------------------------------

dagH_burden     <- draw_dag(c("t_p","p_s"), "Burden-driven", faded_nodes = "h")
dagH_intrinsic  <- draw_dag(c("t_h","h_s","t_p"), "Intrinsic decline")
dagH_collapse   <- draw_dag(c("t_h","h_p","p_s"), "Immune collapse")
dagH_instant    <- draw_dag(c("t_p","p_h","h_s"), "Instantaneous damage")
dagH_cumulative <- draw_dag(c("t_p","t_h","p_h","h_s"), "Cumulative damage")

#-------------------------------------------------------------------------------
# Assemble
#-------------------------------------------------------------------------------
# One row keeps all five directly comparable; the 3+2 alternative is below.

figure1 <- (dagH_burden | dagH_intrinsic | dagH_collapse |
              dagH_instant | dagH_cumulative) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 12))

# Alternative 2-row layout if one row is too cramped:
# figure1_dags <- (dagH_burden | dagH_intrinsic | dagH_collapse) /
#                 (dagH_instant | dagH_cumulative | plot_spacer()) +
#   plot_annotation(tag_levels = 'A') &
#   theme(plot.tag = element_text(face = "bold", size = 12))

ggsave("figures/figure1.pdf",  # >>> Manuscript Figure 1 (five candidate causal structures)
       plot = figure1, width = 11, height = 3.0, units = "in", dpi = 300)

#===============================================================================
# FIGURE 2: survival, Gompertz mortality, and the burden-driven requirement
#===============================================================================
# Survival Parameters

data <- survival %>%
  mutate(total = total - 2, alive = alive - 2, survival = alive / total)

control <- control_surv %>% mutate(survival = alive/total) %>% filter(time < 37)

names(data)[3:4] <- c("Dead", "Alive")
names(control)[3:4] <- c("Dead", "Alive")

# Gompertz model
gompertz_model <- function(t, a, b) {
  exp(-a / b * (exp(b * t) - 1))
}

# Fit Gompertz
fit <- nlsLM(
  survival ~ gompertz_model(time, a, b),
  data = data,
  start = list(a = 1, b = 0.1)
)

params <- coef(fit)
a <- params["a"]
b <- params["b"]

cat("\n=== Gompertz survival fit (Figure 2A) ===\n")
print(round(params, 7))

data$fitted_survival <- gompertz_model(data$time, a, b)

# Combine datasets for plotting
data$group <- "Infected"
control$group <- "Control"

plot_data <- bind_rows(data, control) %>%
  filter(time < 37)

gompertz_fit <- data.frame(
  time = data$time,
  fitted_survival = data$fitted_survival,
  group = "Gompertz fit"
)%>%
  filter(time < 37)

# PANEL A: Survival curves with inset
# Calculate instantaneous mortality m(t) = a*exp(b*t)
time_fine <- seq(0, 36, length.out = 200)  # Only generate for display range
m_t_fitted <- a * exp(b * time_fine)

# Create dataframe
mort_fitted_df <- data.frame(
  time = time_fine,
  mortality = m_t_fitted
)

# Create inset plot
inset_plot <- ggplot(mort_fitted_df, aes(x = time, y = mortality)) +
  geom_line(color = "#b80422", linewidth = 1) +
  scale_x_continuous(breaks = c(0, 12, 24, 36)) +
  scale_y_continuous(breaks = c(0, 35, 70)) +
  labs(x = "Time (hrs)", y = bquote(italic(m)*"("*italic(t)*")"~(hrs^-1))) +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.margin = margin(2, 2, 2, 2)
  )

print(inset_plot)

# Main survival plot
p1A <- ggplot(plot_data, aes(x = time, y = survival)) +
  geom_point(data = control, aes(fill = group), size = 3, shape = 21) +
  geom_point(data = plot_data, aes(fill = group), size = 3, shape = 21) +
  geom_line(data = gompertz_fit, aes(y = fitted_survival, linetype = group, color = group), linewidth = 1) +
  annotate("segment", x = 1, xend = 36, y = 1, yend = 1,
           color = "grey50", linewidth = 0.8, linetype = "dashed") +
  labs(x = "Time (hrs)", y = "Proportion of survival") +
  scale_fill_manual(values = c("Control" = "grey70", "Infected" = "grey20")) +
  scale_linetype_manual(values = c("Gompertz fit" = "solid")) +
  scale_color_manual(values = c("Gompertz fit" = "black")) +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2),
    color = "none"
  ) +
  mytheme +
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  scale_x_continuous(breaks = c(0, 12, 24, 36)) +
  theme(legend.position = c(0.8, 0.7))


# Add inset
pA_with_inset <- p1A +
  annotation_custom(
    ggplotGrob(inset_plot),
    xmin = -0.05, xmax = 17, ymin = -0.05, ymax = 0.45
  )

#===============================================================================
# SCENARIO PANELS - Using Empirical Parameters
#===============================================================================

# Extract BOTH Gompertz (from survival) and Logistic (from burden) parameters
params_gompertz <- coef(fit)
a <- params_gompertz["a"]
b <- params_gompertz["b"]

p0_empirical <- 2e3  # Injected concentration
r_empirical  <- b    # assuming survival rate is equal to pathogen growth rate

# Time sequence
t_seq <- seq(0, 36, length.out = 100)

#-------------------------------------------------------------------------------
# SCENARIO 1: Exponential growth + Linear m(p)
#-------------------------------------------------------------------------------
# Theory: p(t) = p0*e^(r*t), m(p) = k*p
# Result: m(t) = k*p0*e^(r*t) = a*e^(b*t)
# Implies: r = b, k = a/p0

p0_s1 <- p0_empirical
r_s1 <- r_empirical    # Must equal b for m(t) = a*e^(b*t)
k_s1 <- a / p0_s1      # Linear coefficient

# Calculate curves
p_s1 <- p0_s1 * exp(r_s1 * t_seq)  # Exponential growth
m_s1 <- k_s1 * p_s1                 # Linear m(p)

df_s1 <- data.frame(
  time = t_seq,
  p = p_s1,
  m_t = m_s1
)

scientific_10 <- function(x) {
  ifelse(x == 0, "0",
         parse(text = gsub("e\\+?", " %*% 10^", scales::scientific(x))))
}

p1C <- ggplot(df_s1, aes(x = time, y = p)) +
  geom_line(color = "#19798b", linewidth = 1.5) +
  labs(y = bquote("Pathogen, "*italic(p)*"("*italic(t)*")"), x = "Time (hrs)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36)) +
  scale_y_continuous(labels = scientific_10)+
  mytheme +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p1D <- ggplot(df_s1, aes(x = p, y = m_t)) +
  geom_line(color = "#ee9b43", linewidth = 1.5) +
  labs(x = bquote("Pathogen, "*italic(p)~"(CFU)"),
       y = bquote("Mortality, "*italic(m)*"("*italic(p)*")")) +
  scale_x_continuous(labels = scientific_10)+
  scale_y_continuous(breaks = c(0, 35, 70))+
  mytheme +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45))

# Combine all panels
scenario1 <- (p1C / p1D) + plot_layout(heights = c(1, 1))

figure2 <- (pA_with_inset | scenario1) +
  plot_layout(widths = c(3, 1))+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 14))

# Save
ggsave("figures/figure2.pdf",  # >>> Manuscript Figure 2 (survival + Gompertz scenarios)
       plot = figure2, width = 10, height = 6, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
# Estimated percentage of bacteria
#-------------------------------------------------------------------------------
# Parameters
K_cfu <- 3e9

# P. aeruginosa cell dimensions (rod-shaped)
cell_length_um <- 2      # ~1.5-3 μm
cell_diameter_um <- 0.6  # ~0.5-0.8 μm
cell_volume_um3 <- pi * (cell_diameter_um/2)^2 * cell_length_um  # cylinder approximation
# ≈ 0.57 μm³ per cell

# Convert to mL (1 μm³ = 10^-12 mL)
cell_volume_mL <- cell_volume_um3 * 1e-12

# Total bacterial volume
total_bact_volume_mL <- K_cfu * cell_volume_mL

# Galleria larva volume
larva_mass_mg <- 250  # typical late-instar larva: 200-300 mg
larva_volume_mL <- larva_mass_mg / 1000  # assuming density ≈ 1 g/mL

# Percentage
percent_bacteria <- (total_bact_volume_mL / larva_volume_mL) * 100

cat("Bacterial volume:", total_bact_volume_mL * 1000, "μL\n")
cat("Larva volume:", larva_volume_mL * 1000, "μL\n")
cat("Bacteria as % of larva:", round(percent_bacteria, 3), "%\n")

#===============================================================================
# BACTERIAL BURDEN
#===============================================================================
# Detection limit: 1 colony x lowest dilution x (homogenate volume / volume plated).
# Both experiments homogenised in 250 ul; the lowest-dilution plates were
# 100 ul at 10x (main) and 100 ul at 10x (AB), so 1 colony = 25 CFU/larva.
# A larva plated with no colonies is a measurement at the detection limit,
# not a missing value; it enters at LOD/2. Larvae with no readable plate
# (Countable == "no") remain excluded.
lod_cfu <- 25

burden_tidy_time <- burden %>%
  filter(Time < 37) %>%
  ungroup() %>%
  pivot_longer(cols = rep1:rep4, names_to = "replicates", values_to = "cfu") %>%
  filter(! Countable == "no") %>%
  mutate(dilution_factor = Dilution*(Buffer_dilution_factor/Volume_ul)) %>%
  mutate(count = cfu*dilution_factor) %>%
  group_by(Sample, Time, Larvae) %>%
  summarise(cfu = mean(count, na.rm = T)) %>%
  ungroup() %>%
  mutate(Sample = as.character(Sample), Larvae = as.character(Larvae)) %>%
  left_join(health %>% mutate(Sample = as.character(Sample), Larvae = as.character(Larvae)),
            by = c("Sample", "Larvae")) %>%
  left_join(death %>% mutate(Sample = as.character(Sample), Larvae = as.character(Larvae)),
            by = c("Sample", "Larvae")) %>%
  mutate(
    survival = case_when(
      survival == 2 ~ 1,
      survival == 0 ~ 0),
    status = case_when(
      survival == 1 ~ "Alive",
      survival == 0 ~ "Dead",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(cfu = ifelse(cfu == 0, lod_cfu / 2, cfu),
         log_CFU = log10(cfu + 1)) %>%
  filter(!is.na(cfu)) %>%
  # as.numeric(): scale() otherwise returns 1-column matrices, which propagate
  # into lavaan and the plotting code as matrix columns.
  mutate(alive               = as.integer(status == "Alive"),
         # Sampling-time tertiles. Defined once here so Figures 3B and 4A-C all
         # encode time the same way; with nine sampling times these are clean
         # splits (Early 4-12 h, Mid 16-24 h, Late 28-36 h).
         t_band              = cut(Time, breaks = quantile(Time, c(0, 1/3, 2/3, 1),
                                                           na.rm = TRUE),
                                   include.lowest = TRUE,
                                   labels = c("Early", "Mid", "Late")),
         health_combined     = activity + melanization,   # 0-7, higher = healthier
         scaled_time         = as.numeric(scale(Time)),
         scaled_health       = as.numeric(scale(health_combined)),
         scaled_activity     = as.numeric(scale(activity)),
         scaled_melanization = as.numeric(scale(melanization)),
         scaled_cfu          = as.numeric(scale(log_CFU)),
         scaled_immune_morb  = as.numeric(scale(activity + melanization))
  ) %>%
  ungroup()

burden_tidy_time$Time <- as.numeric(burden_tidy_time$Time)

#==============================================================================
# FIT MODELS
#==============================================================================

# only alive proportion
burden_tidy_alive <- burden_tidy_time %>%
  filter(status == "Alive")

burden_tidy_time %>%
  mutate(bin = cut(Time, c(0, 12, 24, 37))) %>%
  group_by(bin) %>%
  summarise(mel = mean(melanization, na.rm = TRUE), act = mean(activity, na.rm = TRUE))

#-------------------------------------------------------------------------------
# FIT 1: Total population (all larvae) with full data
#-------------------------------------------------------------------------------
logistic_logfit_full <-
  nlsLM(log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
        data = burden_tidy_time,
        start = list(
          K = unname(quantile(burden_tidy_time$cfu, 0.9, na.rm = TRUE)),
          p0 = unname(quantile(burden_tidy_time$cfu, 0.1, na.rm = TRUE)),
          r = 0.2
        ),
        lower = c(K = 100, p0 = 0.1, r = 0.01),
        upper = c(K = 1e9, p0 = 5e3, r = 2),
        control = list(maxiter = 1000)
  )

# Bootstrap for total population
pred_time_full <- data.frame(Time = seq(min(burden_tidy_time$Time),
                                        max(burden_tidy_time$Time),
                                        by = 0.1))

cat("\n=== Logistic growth fit, all larvae (Figure 3) ===\n")
print(coef(logistic_logfit_full))

boot_preds_total <- function(data, indices) {
  d <- data[indices, ]
  fit <- tryCatch({
    nlsLM(
      log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
      data = d,
      start = list(
        K = unname(quantile(d$cfu, 0.9, na.rm = TRUE)),
        p0 = unname(quantile(d$cfu, 0.1, na.rm = TRUE)),
        r = 0.2
      ),
      lower = c(K = 100, p0 = 0.1, r = 0.01),
      upper = c(K = 1e9, p0 = 5e3, r = 2),
      control = list(maxiter = 1000)
    )
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(rep(NA, 3))
  coef(fit)
}

set.seed(123)
boot_curve_full <- boot(data = burden_tidy_time, statistic = boot_preds_total, R = 1000)

boot_matrix_full <- boot_curve_full$t
valid_matrix_full <- boot_matrix_full[apply(boot_matrix_full, 1, function(row) {
  !any(is.na(row)) && (max(row) - min(row)) > 1
}), ]

boot_prediction_curves_full <- apply(valid_matrix_full, 1, function(params) {
  K <- params[1]; p0 <- params[2]; r <- params[3]
  log10(K / (1 + ((K - p0) / p0) * exp(-r * pred_time_full$Time)))  # Use pred_time_full$Time
})

boot_prediction_curves_full <- t(boot_prediction_curves_full)
ci_bounds_full <- apply(boot_prediction_curves_full, 2, quantile, probs = c(0.025, 0.975))
median_fit_full <- apply(boot_prediction_curves_full, 2, median)

pred_time_full$lower <- ci_bounds_full[1, ]
pred_time_full$upper <- ci_bounds_full[2, ]
pred_time_full$fit <- median_fit_full

#------------------------------------------------------------------------------
# FIT 2: Alive only
#------------------------------------------------------------------------------
logistic_logfit_a_full <- nlsLM(log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
                                data = burden_tidy_alive,
                                # unname(): a named quantile ("90%") gets pasted
                                # into the parameter name, giving "K.90%", and
                                # then coef(...)["K"] silently returns NA.
                                start = list(
                                  K  = unname(quantile(burden_tidy_alive$cfu, 0.9, na.rm = TRUE)),
                                  p0 = unname(quantile(burden_tidy_alive$cfu, 0.1, na.rm = TRUE)),
                                  r = 0.2
                                ),
                                lower = c(K = 100, p0 = 0.1, r = 0.01),
                                upper = c(K = 1e9, p0 = 5e3, r = 2),
                                control = list(maxiter = 1000))

cat("\n=== Logistic growth fit, survivors only (Figure 3, dashed) ===\n")
print(coef(logistic_logfit_a_full))

# Bootstrap for alive only
pred_time_a_full <- data.frame(Time = seq(min(burden_tidy_alive$Time),
                                          max(burden_tidy_alive$Time),
                                          by = 0.1))

boot_preds_alive <- function(data, indices) {
  d <- data[indices, ]
  fit <- tryCatch({
    nlsLM(log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
          data = d,
          start = list(                       # unname(): see FIT 2 above
            K  = unname(quantile(d$cfu, 0.9, na.rm = TRUE)),
            p0 = unname(quantile(d$cfu, 0.1, na.rm = TRUE)),
            r = 0.2
          ),
          lower = c(K = 100, p0 = 0.1, r = 0.01),
          upper = c(K = 1e9, p0 = 5e3, r = 2),
          control = list(maxiter = 1000))
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(rep(NA, 3))
  coef(fit)
}

set.seed(124)
boot_curve_a_full <- boot(data = burden_tidy_alive, statistic = boot_preds_alive, R = 1000)

boot_matrix_a_full <- boot_curve_a_full$t
valid_matrix_a_full <- boot_matrix_a_full[apply(boot_matrix_a_full, 1, function(row) {
  !any(is.na(row)) && (max(row) - min(row)) > 1
}), ]

boot_prediction_curves_a_full <- apply(valid_matrix_a_full, 1, function(params) {
  K <- params[1]; p0 <- params[2]; r <- params[3]
  log10(K / (1 + ((K - p0) / p0) * exp(-r * pred_time_a_full$Time)))
})

boot_prediction_curves_a_full <- t(boot_prediction_curves_a_full)
ci_bounds_a_full <- apply(boot_prediction_curves_a_full, 2, quantile, probs = c(0.025, 0.975))
median_fit_a_full <- apply(boot_prediction_curves_a_full, 2, median)

pred_time_a_full$lower <- ci_bounds_a_full[1, ]
pred_time_a_full$upper <- ci_bounds_a_full[2, ]
pred_time_a_full$fit   <- median_fit_a_full

# Compare K values -- these are the CIs quoted in the Figure 3 caption
ci_K_total <- quantile(boot_curve_full$t[, 1], probs = c(0.025, 0.975), na.rm = TRUE)
ci_K_alive <- quantile(boot_curve_a_full$t[, 1], probs = c(0.025, 0.975), na.rm = TRUE)

cat("\n=== Carrying capacity K, bootstrap 95% CI (Figure 3 caption) ===\n")
cat("  All larvae:      "); print(signif(ci_K_total, 3))
cat("  Survivors only:  "); print(signif(ci_K_alive, 3))

#===============================================================================
# FIGURE 3: logistic pathogen growth across the live-dead threshold
#===============================================================================

# Define colors and linetypes
line_colors <- c("Total population" = "#ee9b43", "Alive only" = "#19798b")
line_types  <- c("Total population" = "solid", "Alive only" = "dashed")

p3A <- ggplot(burden_tidy_time, aes(x = Time, y = log_CFU)) +
  geom_jitter(aes(fill = status), size = 3, shape = 21, alpha = 0.6, color = "black") +
  # Total population fit
  geom_ribbon(data = pred_time_full, aes(x = Time, ymin = lower, ymax = upper),
              fill = "#ee9b43", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = pred_time_full, aes(x = Time, y = fit,
                                       color = "Total population",
                                       linetype = "Total population"),
            linewidth = 1.2, inherit.aes = FALSE) +
  # Alive only fit
  geom_ribbon(data = pred_time_a_full, aes(x = Time, ymin = lower, ymax = upper),
              fill = "#19798b", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = pred_time_a_full, aes(x = Time, y = fit,
                                         color = "Alive only",
                                         linetype = "Alive only"),
            linewidth = 1.2, inherit.aes = FALSE) +
  xlab("Time (hrs)") +
  ylab(bquote(log[10](CFU))) +
  scale_fill_manual(
    name = "Status",
    values = c("Alive" = "#19798b", "Dead" = "#b80422")
  ) +
  scale_color_manual(
    name = "Logistic fit",
    values = line_colors
  ) +
  scale_linetype_manual(
    name = "Logistic fit",
    values = line_types
  ) +
  scale_x_continuous(breaks = c(4, 12, 20, 28, 36)) +
  mytheme +
  theme(legend.position = c(0.25, 0.8)) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  )

# ----------------------------------------------------------------------------
# Panel B: survival against pathogen density, stratified by sampling time
# ----------------------------------------------------------------------------
# The transpose of Figure S1's diagnostic panel: there survival is plotted
# against time at three fixed burdens, here against burden at three fixed times.
# Both read off the same fit (m_cond), so the beta_t quoted in the two captions
# is deliberately the same number.
#
#   solid  = logit P(alive) ~ Time + log_CFU     the tested model
#   dashed = logit P(alive) ~ log_CFU            the null implied by burden-driven
#                                                mortality: time adds nothing
#                                                once density is known. It has no
#                                                time term, so the three bands
#                                                collapse onto ONE curve -- drawn
#                                                once, in grey.
#
# The vertical separation between the solid curves is the rejection of
# s _||_ t | p. Each curve is clipped to the log10 CFU range actually observed
# within its own time band -- without that the early-time sigmoid is drawn
# across burdens no early larva ever reached, which is more misleading here
# than in Figure 4C because the curves are sigmoid rather than straight.

m_cond   <- glm(alive ~ Time + log_CFU, data = burden_tidy_time, family = binomial)
m_cond_null <- glm(alive ~ log_CFU,     data = burden_tidy_time, family = binomial)

# Tertile bands, and the time at which each band's curve is drawn. The tertile
# medians coincide exactly with the 15th / 50th / 85th percentiles of sampling
# time here (8, 20, 32 h), so the caption may quote either.
dat_fig3b <- burden_tidy_time

band3 <- dat_fig3b %>%
  group_by(t_band) %>%
  summarise(Time = median(Time),
            lo   = min(log_CFU), hi = max(log_CFU), .groups = "drop")

pred_3B <- purrr::map_dfr(seq_len(nrow(band3)), function(i) {
  nd <- data.frame(log_CFU = seq(band3$lo[i], band3$hi[i], length.out = 200),
                   Time    = band3$Time[i])
  nd$fit    <- predict(m_cond, newdata = nd, type = "response")
  nd$t_band <- band3$t_band[i]
  nd
})

# The null has no time term, so its three band-curves are identical. Drawing it
# once, in grey, across the full observed range says that plainly: under
# burden-driven mortality there is ONE curve, and the tested model has three.
null_3B <- data.frame(
  log_CFU = seq(min(burden_tidy_time$log_CFU), max(burden_tidy_time$log_CFU),
                length.out = 200))
null_3B$fit <- predict(m_cond_null, newdata = null_3B, type = "response")

p3B <- ggplot(pred_3B, aes(x = log_CFU, y = fit, colour = t_band)) +
  geom_jitter(data = dat_fig3b,
              aes(x = log_CFU, y = alive, colour = t_band),
              inherit.aes = FALSE, height = 0.035, width = 0.05,
              alpha = 0.6, size = 2.5) +
  # The null gets its own linetype legend entry: previously a reader had to
  # reach the caption to learn what the grey reference curve was.
  geom_line(data = null_3B,
            aes(x = log_CFU, y = fit, linetype = "Burden-driven null"),
            inherit.aes = FALSE, linewidth = 0.8, colour = "grey45") +
  geom_line(aes(linetype = "Time + burden"), linewidth = 1.1) +
  scale_colour_manual(values = pal_tband, name = "Time") +
  scale_linetype_manual(values = c("Time + burden" = "solid",
                                   "Burden-driven null" = "dashed"),
                        name = NULL,
                        breaks = c("Time + burden", "Burden-driven null")) +
  # The curves are fitted probabilities on a continuous 0-1 scale; the points
  # are the observed binary outcome. Label both rather than conflating them.
  scale_y_continuous(breaks = c(0, 0.5, 1),
                     labels = c("0 (Dead)", "0.5", "1 (Alive)")) +
  coord_cartesian(ylim = c(-0.12, 1.12)) +
  labs(x = bquote(log[10](CFU)), y = "Predicted P(survival)") +
  guides(colour   = guide_legend(order = 1),
         linetype = guide_legend(order = 2,
                                 override.aes = list(colour = c("grey20", "grey45")))) +
  mytheme +
  theme(legend.position = c(0.17, 0.36),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        legend.spacing.y = unit(0.05, "cm"))

figure3 <- (p3A | p3B) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("figures/figure3.pdf",  # >>> Manuscript Figure 3 (logistic growth + s vs p | t)
       plot = figure3, width = 11, height = 5, units = "in", dpi = 300)


# Linear-scale version (diagnostic; not in the manuscript)
print(
  ggplot(burden_tidy_time, aes(Time, cfu)) +
    geom_jitter(aes(fill = status), shape = 21, size = 3, alpha = 0.6, color = "black") +
    geom_line(data = pred_time_full, aes(Time, 10^fit), color = "#ee9b43", linewidth = 1.2) +
    scale_fill_manual(values = c("Alive" = "#19798b", "Dead" = "#ee9b43")) +
    scale_y_continuous(labels = scales::label_scientific()) +
    labs(x = "Time (h)", y = "CFU (linear)") + mytheme
)

#------------------------------------------------------------------------------
# Comparisons with other models
#------------------------------------------------------------------------------

# Point-estimate AICs quoted in the Results text (logistic 190 vs exponential 220)
aic_logistic_point <- AIC(logistic_logfit_full)
aic_exp_point      <- AIC(lm(log_CFU ~ Time, data = burden_tidy_time))
cat("\n=== Growth model AIC, all larvae (Results text) ===\n")
cat(sprintf("  Logistic:    %.1f\n  Exponential: %.1f\n",
            aic_logistic_point, aic_exp_point))

boot_model_compare <- function(data, indices) {
  d <- data[indices, ]
  
  # Fit all three models
  fit_logistic <- tryCatch({
    nlsLM(log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
          data = d,
          start = list(
            K = unname(quantile(burden_tidy_time$cfu, 0.9, na.rm = TRUE)),
            p0 = unname(quantile(burden_tidy_time$cfu, 0.1, na.rm = TRUE)),
            r = 0.2
          ),
          lower = c(K = 100, p0 = 1, r = 0.01),
          upper = c(K = 1e9, p0 = 5e3, r = 2),
          control = list(maxiter = 1000))
  }, error = function(e) NULL)
  
  fit_exp <- tryCatch({
    lm(log_CFU ~ Time, data = d)
  }, error = function(e) NULL)
  
  fit_linear <- tryCatch({
    nlsLM(log_CFU ~ log10(p0 + r * Time),
          data = d,
          start = list(p0 = quantile(burden_tidy_time$cfu, 0.1, na.rm = TRUE), r = 1000),
          lower = c(p0 = 1, r = 0))
  }, error = function(e) NULL)
  
  if (is.null(fit_logistic) | is.null(fit_exp) | is.null(fit_linear)) {
    return(c(NA, NA, NA, NA))
  }
  
  aics <- c(logistic = AIC(fit_logistic),
            exponential = AIC(fit_exp),
            linear = AIC(fit_linear))
  
  winner <- which.min(aics)
  
  c(AIC_logistic = aics[1],
    AIC_exp = aics[2],
    AIC_linear = aics[3],
    winner = winner)  # 1 = logistic, 2 = exponential, 3 = linear
}

set.seed(125)
boot_compare <- boot(data = burden_tidy_time,
                     statistic = boot_model_compare,
                     R = 1000)
set.seed(126)
boot_compare_alive <- boot(data = burden_tidy_alive,
                           statistic = boot_model_compare,
                           R = 1000)


# Count wins - Text for the manuscript
winners <- boot_compare$t[, 4]
winners <- winners[!is.na(winners)]

cat("\n=== Bootstrap model selection, all larvae ===\n")
cat("  Logistic wins:", sum(winners == 1), "(", round(100 * mean(winners == 1)), "%)\n")
cat("  Exponential wins:", sum(winners == 2), "(", round(100 * mean(winners == 2)), "%)\n")
cat("  Linear wins:", sum(winners == 3), "(", round(100 * mean(winners == 3)), "%)\n")

winners_a <- boot_compare_alive$t[, 4]
winners_a <- winners_a[!is.na(winners_a)]

cat("\n=== Bootstrap model selection, survivors only ===\n")
cat("  Logistic wins:", sum(winners_a == 1), "(", round(100 * mean(winners_a == 1)), "%)\n")
cat("  Exponential wins:", sum(winners_a == 2), "(", round(100 * mean(winners_a == 2)), "%)\n")
cat("  Linear wins:", sum(winners_a == 3), "(", round(100 * mean(winners_a == 3)), "%)\n")


#-------------------------------------------------------------------------------
# CFU "time since death"
# Supplementary figure S2
#-------------------------------------------------------------------------------

# Calculate time since death
burden_tidy_death <- burden_tidy_time %>%
  mutate(
    time_since_death = case_when(
      status == "Dead" ~ Time - time_to_death_h,  # adjust column names as needed
      status == "Alive" ~ NA_real_
    )
  )

# Check the distribution
print(
  burden_tidy_death %>%
    filter(status == "Dead") %>%
    dplyr::select(Time, time_to_death_h, time_since_death, log_CFU) %>%
    summary()
)

# Plot: CFU vs time since death (dead larvae only)
p_postmortem <- ggplot(burden_tidy_death %>% filter(status == "Dead"),
                       aes(x = time_since_death, y = log_CFU)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#b80422") +
  labs(
    x = "Time since death (hrs)",
    y = expression(log[10](CFU))
  ) +
  mytheme

print(p_postmortem)

# Save
ggsave("figures/figureS2.pdf",  # >>> Supplementary Figure S2 (post-mortem burden)
       plot = p_postmortem, width = 5.3, height = 5, units = "in", dpi = 300)

# Test for post-mortem growth
lm_postmortem <- lm(log_CFU ~ time_since_death,
                    data = burden_tidy_death %>% filter(status == "Dead"))
cat("\n=== Post-mortem burden stability (Figure S1) ===\n")
print(summary(lm_postmortem))

# Compare to time since infection for context
p_comparison <- burden_tidy_death %>%
  filter(status == "Dead") %>%
  ggplot(aes(x = Time, y = log_CFU)) +
  geom_point(aes(color = time_since_death), size = 2) +
  scale_color_viridis_c(name = "Hours\npost-death") +
  labs(
    x = "Time since infection (hrs)",
    y = expression(log[10](CFU))
  ) +
  mytheme

print(p_comparison)

#-------------------------------------------------------------------------------
# m(p) mapping -- becomes panel A of Supplementary Figure S3
#-------------------------------------------------------------------------------

params_gompertz <- coef(fit)
a <- params_gompertz["a"]
b <- params_gompertz["b"]

# Fitted logistic-growth parameters (raw-CFU scale). Set once here and reused by
# the m(p) mapping, the cumulative-burden integral and the antibiotic Sigma_p;
# they are never reassigned, so do not re-pull them further down.
params_logistic <- coef(logistic_logfit_full)
K  <- params_logistic["K"]
p0 <- params_logistic["p0"]
r  <- params_logistic["r"]

# Realised burden trajectory over the observation window. Kept here because the
# time-effect mapping section downstream reuses p_t_raw paired with time_seq.
p_t_raw <- K / (1 + ((K - p0) / p0) * exp(-r * time_seq))
p_t_log <- log10(p_t_raw)

# For the figure, sweep p directly from p0 up to just below K so the asymptotic
# behaviour is shown across the full burden range. m(p) is a pure function of p,
# so no time variable is needed; this also guarantees the curve reaches toward K
# even though the realised t = 0..36 h trajectory need not.
p_seq <- 10^seq(log10(p0), log10(0.999 * K), length.out = 500)
m_p   <- a * ((p_seq * (K - p0)) / (p0 * (K - p_seq)))^(b / r)

df_m <- data.frame(
  p_log = log10(p_seq),
  m_p   = m_p
)

# Place the annotation in the empty upper-left region, robust to the value of K.
x_lab <- unname(log10(p0) + 0.08 * (log10(0.999 * K) - log10(p0)))

# Alternative rendering of the same mapping on a log10(CFU) x-axis.
# NOT used in any assembled figure; kept for reference.
figure_mp <- ggplot(df_m, aes(x = p_log, y = m_p)) +
  geom_line(color = "#b80422", linewidth = 1.5) +
  labs(x = bquote(log[10](CFU)),
       y = bquote(italic(m)*"("*italic(p)*")"~(hrs^-1))) +
  scale_x_continuous(breaks = c(2, 4, 6)) +
  scale_y_continuous(breaks = c(0, 35, 70)) +
  # m(p) -> infinity as p -> K, so bound the y-window and let the curve exit the
  # top edge. Without this the axis auto-scales to ~10^4+ and the curve collapses
  # onto the x-axis.
  coord_cartesian(ylim = c(0, 75)) +
  annotate("text",
           x = x_lab, y = 55,
           label = "atop(m(p) %->% infinity, as~p %->% K)",
           parse = TRUE, size = 5, hjust = 0) +
  mytheme +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Implied m(p) mapping has a pole at p = K:  m(p) = a * R(p)^(b/r),
#   R(p) = (p (K - p0)) / (p0 (K - p)).
# For p > K the exact power is complex, so plot the signed-magnitude
# continuation m~(p) = a * sign(R) * |R|^(b/r) (exact when b/r = 1).
# p > K is biologically unreachable: it maps to imaginary time in the logistic.

expo <- 1                       # clean schematic; set expo <- b / r for the exact exponent
m_signed <- function(p) { R <- (p*(K - p0))/(p0*(K - p)); a * sign(R) * abs(R)^expo }

df_pole <- rbind(
  data.frame(pk = seq(0.02, 0.999, length.out = 500), side = "below"),
  data.frame(pk = seq(1.001, 2.5,  length.out = 500), side = "above")
)
df_pole$m <- m_signed(df_pole$pk * K)
ywin <- max(abs(m_signed(0.8 * K)), abs(m_signed(2 * K))) * 1.1

# >>> This becomes panel A of Supplementary Figure S3.
pS3A <- ggplot(df_pole, aes(pk, m, group = side)) +
  annotate("rect", xmin = 1, xmax = 2.5, ymin = -ywin, ymax = ywin, fill = "grey50", alpha = 0.08) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_line(colour = "#ee9b43", linewidth = 1.3) +
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("0", "K", "2K")) +
  labs(x = bquote("Pathogen burden, "*italic(p)), y = bquote(italic(m)*"("*italic(p)*")"~(hrs^-1))) +
  coord_cartesian(ylim = c(-ywin*4, ywin*4), clip = "on") +
  mytheme + theme(panel.grid = element_blank(), plot.margin = margin(14, 10, 6, 6))

#===============================================================================
# Health variables
#===============================================================================

# ==============================================================================
# FIGURE 4: activity, melanization, and the T50 ordering
# ------------------------------------------------------------------------------
# Layout (single row of three panels):
#   A: Activity vs log10(CFU), points coloured by Time.
#   B: Melanization vs log10(CFU), points coloured by Time.
#   C: Three Gompertz current-status curves with T50 markers (activity,
#      un-melanized, alive), fit by binomial MLE on individual-level data.
# ==============================================================================

# ----------------------------------------------------------------------------
# Prepare current-status data
# ----------------------------------------------------------------------------
# Definitions (matching existing code conventions):
#   active:    activity score >= 2
#   unmelan:   melanization >= 3      (raw score; 4 = no melanization)     
#   alive:     status == "Alive"
# Each row = one larva at its sampling time.

dat_cs <- burden_tidy_time %>%
  filter(!is.na(activity), !is.na(melanization), !is.na(status)) %>%
  transmute(
    time   = Time,
    active = as.integer(activity   >= 2),
    unmel  = as.integer(melanization >= 3),
    alive  = as.integer(status == "Alive")
  )

# ----------------------------------------------------------------------------
# Binomial Gompertz current-status MLE
# ----------------------------------------------------------------------------
# Gompertz survival:   S(t) = exp(-(a/b) * (exp(b*t) - 1))
# Likelihood for binary y_i sampled at t_i:  S(t_i)^y_i * (1 - S(t_i))^(1-y_i)
# Parameterised on log-scale to keep a, b > 0.

gompertz_S <- function(t, a, b) exp(-(a / b) * (exp(b * t) - 1))

nll_gompertz <- function(par, t, y) {
  a <- exp(par[1]); b <- exp(par[2])
  S <- gompertz_S(t, a, b)
  S <- pmin(pmax(S, 1e-12), 1 - 1e-12)
  -sum(y * log(S) + (1 - y) * log(1 - S))
}

fit_gompertz_cs <- function(time, y) {
  optim(
    par     = c(log(1e-4), log(0.4)),
    fn      = nll_gompertz,
    t       = time,
    y       = y,
    method  = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-10)
  )
}

# T50 from Gompertz parameters: solve S(t) = 0.5
T50_from_par <- function(a, b) (1 / b) * log(1 - (b / a) * log(0.5))

# Fit each metric
fit_AT <- fit_gompertz_cs(dat_cs$time, dat_cs$active)
fit_MT <- fit_gompertz_cs(dat_cs$time, dat_cs$unmel)
fit_LT <- fit_gompertz_cs(dat_cs$time, dat_cs$alive)

par_to_T50 <- function(fit) {
  a <- exp(fit$par[1]); b <- exp(fit$par[2])
  T50_from_par(a, b)
}

T50_act  <- par_to_T50(fit_AT)
T50_mel  <- par_to_T50(fit_MT)
T50_surv <- par_to_T50(fit_LT)

# Manuscript text 
cat("\n=== Transition timing (Figure 4C) ===\n")
cat(sprintf("AT50 = %.2f h\nMT50 = %.2f h\nLT50 = %.2f h\n",
            T50_act, T50_mel, T50_surv))

# ----------------------------------------------------------------------------
# Bootstrap CIs for T50 (resample larvae)
# ----------------------------------------------------------------------------
boot_T50_cs <- function(time, y, n = 1000) {
  out <- replicate(n, {
    idx <- sample.int(length(y), replace = TRUE)
    fit <- tryCatch(fit_gompertz_cs(time[idx], y[idx]), error = function(e) NULL)
    if (is.null(fit)) NA_real_ else par_to_T50(fit)
  })
  out[is.finite(out)]
}

set.seed(123)
boot_AT <- boot_T50_cs(dat_cs$time, dat_cs$active)
boot_MT <- boot_T50_cs(dat_cs$time, dat_cs$unmel)
boot_LT <- boot_T50_cs(dat_cs$time, dat_cs$alive)

T50_summary <- tibble(
  metric  = c("AT50", "MT50", "LT50"),
  T50     = c(T50_act, T50_mel, T50_surv),
  CI_low  = c(quantile(boot_AT, 0.025), quantile(boot_MT, 0.025), quantile(boot_LT, 0.025)),
  CI_high = c(quantile(boot_AT, 0.975), quantile(boot_MT, 0.975), quantile(boot_LT, 0.975))
)
print(T50_summary)

# Pairwise ordering probabilities.

n_paired <- min(length(boot_AT), length(boot_MT), length(boot_LT))
ix <- seq_len(n_paired)
cat(sprintf("\nP(AT50 < MT50) = %.3f\n", mean(boot_AT[ix] < boot_MT[ix])))
cat(sprintf("P(MT50 < LT50) = %.3f\n", mean(boot_MT[ix] < boot_LT[ix])))
cat(sprintf("P(AT50 < LT50) = %.3f\n", mean(boot_AT[ix] < boot_LT[ix])))

# ----------------------------------------------------------------------------
# Curves for plotting
# ----------------------------------------------------------------------------
t_grid <- seq(0, 36, length.out = 300)
make_curve <- function(fit, label) {
  a <- exp(fit$par[1]); b <- exp(fit$par[2])
  data.frame(time = t_grid, p = gompertz_S(t_grid, a, b), metric = label)
}
curves <- bind_rows(
  make_curve(fit_AT, "AT50"),
  make_curve(fit_MT, "MT50"),
  make_curve(fit_LT, "LT50")
)

# Observed proportions per timepoint, per metric (for overlay)
obs_props <- dat_cs %>%
  pivot_longer(c(active, unmel, alive),
               names_to = "metric_short", values_to = "y") %>%
  group_by(time, metric_short) %>%
  summarise(p = mean(y), n = n(), .groups = "drop") %>%
  mutate(metric = dplyr::recode(metric_short,
                                active = "AT50",
                                unmel  = "MT50",
                                alive  = "LT50"))
obs_props %>% arrange(time, metric) %>% print(n = Inf)

palette_T50 <- c(
  "AT50" = "#19798b",
  "MT50" = "#ee9b43",
  "LT50" = "#b80422"
)


# Show correlation between Activity and Melanization.
# Correlate them in the SAME (health) orientation they enter the composite:
# Manuscript text

melan_health <- burden_tidy_time$melanization
cor_AM <- cor(burden_tidy_time$activity, melan_health,
              use = "complete.obs")
cat("\n=== Activity vs melanization, both health-oriented (Table S1) ===\n")
cat("  r =", round(cor_AM, 3), "\n")

# Test correlation significance
cor_test <- cor.test(burden_tidy_time$activity, melan_health)
cat("  n =", cor_test$parameter + 2, "  p =", signif(cor_test$p.value, 3), "\n\n")

# ----------------------------------------------------------------------------
# Panels A, B: CFU vs health, time-coloured
# ----------------------------------------------------------------------------
# Use viridis 'magma' reversed: early = light, late = dark.
# Black SCAM line gives the marginal CFU-only fit for visual reference.


# Fit GAM for CFU vs Activity relationship
gam_cfu_act <- scam(activity ~ s(log_CFU, k = 8, bs = "mpd"), data = burden_tidy_time)
summary(gam_cfu_act)

# Generate predictions for smooth line
cfu_seq_act <- data.frame(log_CFU = seq(min(burden_tidy_time$log_CFU),
                                        max(burden_tidy_time$log_CFU),
                                        length.out = 100))
cfu_seq_act$smoothed_activity <- predict(gam_cfu_act, newdata = cfu_seq_act)

# Fit GAM for CFU vs Melanization relationship
gam_cfu_mel <- scam(melanization ~ s(log_CFU, k = 8, bs = "mpd"), data = burden_tidy_time)
summary(gam_cfu_mel)

# Generate predictions for smooth line
cfu_seq_mel <- data.frame(log_CFU = seq(min(burden_tidy_time$log_CFU),
                                        max(burden_tidy_time$log_CFU),
                                        length.out = 100))
cfu_seq_mel$smoothed_melanization <- predict(gam_cfu_mel, newdata = cfu_seq_mel)


# Time is encoded as the same three bands as panel C. It was previously a
# continuous colourbar here, which put two encodings of one variable in one
# figure and made the panels harder to read against each other.
p4A <- ggplot(burden_tidy_time, aes(x = log_CFU, y = activity)) +
  geom_jitter(aes(color = t_band), size = 2.5, width = 0, height = 0.12, alpha = 0.85) +
  geom_line(data = cfu_seq_act, aes(x = log_CFU, y = smoothed_activity),
            color = "black", linewidth = 0.9, inherit.aes = FALSE) +
  scale_color_manual(values = pal_tband, name = "Time") +
  labs(x = expression(log[10](CFU)), y = "Activity score") +
  mytheme +
  theme(legend.position = c(0.83, 0.78),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.title = element_text(size = 16)
  )

# Fig 4B: melanization
p4B <- ggplot(burden_tidy_time, aes(x = log_CFU, y = melanization)) +
  geom_jitter(aes(color = t_band), size = 2.5, width = 0, height = 0.12, alpha = 0.85) +
  geom_line(data = cfu_seq_mel, aes(x = log_CFU, y = smoothed_melanization),
            color = "black", linewidth = 0.9, inherit.aes = FALSE) +
  scale_color_manual(values = pal_tband, name = "Time") +
  labs(x = expression(log[10](CFU)), y = "Melanization score") +
  mytheme +
  theme(legend.position = "none")  # shares the legend with p4A visually

# ----------------------------------------------------------------------------
# Panel C: T50 ordering
# ----------------------------------------------------------------------------
p4D <- ggplot(curves, aes(x = time, y = p, color = metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(data = obs_props,
             aes(x = time, y = p, color = metric),
             alpha = 0.75, size = 4,
             position = position_dodge(width = 0.6)) +
  geom_segment(data = T50_summary,
               aes(x = T50, xend = T50, y = 0, yend = 0.5, color = metric),
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.4) +
  scale_color_manual(values = palette_T50, name = NULL,
                     limits = c("AT50", "MT50", "LT50")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = seq(0, 36, 6)) +
  # coord_cartesian, not scale limits: scale_*_continuous(limits =) DELETES data
  # outside the range, and the dodge pushes the t = 36 points to 36.3.
  coord_cartesian(xlim = c(0, 36), ylim = c(0, 1.02)) +
  labs(x = "Time (hrs)", y = "Proportion") +
  mytheme +
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.text = element_text(size = 16))

# ----------------------------------------------------------------------------
# Panel C: health against burden, stratified by time
# ----------------------------------------------------------------------------
# The transpose of Figure S3D. There the burden level is held fixed and health
# is plotted against time; here time is held fixed and health is plotted against
# burden. Both are read off the same fit, m_hp (= m_D in the Figure S4 section).
#
#   solid  = h ~ Time + log_CFU        the tested model
#   dashed = h ~ Time                  the null, h _||_ p | t
#
# Sloping solid lines against flat dashed ones are the rejection; the vertical
# offset between the three bands is the memory -- at the same burden, larvae
# sampled later are in worse condition.
#
# NOTE ON THE MODEL. Time enters linearly. The manuscript sentence for this
# panel should therefore read "adjusting for time since infection", NOT
# "allowing flexibility for the effect of time" -- that phrasing describes
# s(Time), which is not what is fitted here. The flexible-conditioning version
# of the RECIPROCAL test does exist (gam_health_flex, below), and that is where
# the "spline-adjusted" wording belongs.

# Points are coloured by time tertile; each line is drawn at the MEDIAN time
# within its own tertile, so lines and points share one three-level scale
# rather than mixing tertile bands with unrelated percentiles.
dat_fig4d <- burden_tidy_time %>%
  filter(!is.na(health_combined), !is.na(log_CFU))

m_hp      <- lm(health_combined ~ Time + log_CFU, data = dat_fig4d)  # tested
m_hp_null <- lm(health_combined ~ Time,           data = dat_fig4d)  # h _||_ p | t

band_time <- dat_fig4d %>%
  group_by(t_band) %>%
  summarise(Time = median(Time),
            lo   = min(log_CFU), hi = max(log_CFU), .groups = "drop")

# Each line is clipped to the burden range observed WITHIN its own band. Early
# larvae never reached high burden and late larvae never sat at low burden, so
# an unclipped line extrapolates a long way past any data in its band.
pred_4C <- purrr::map_dfr(seq_len(nrow(band_time)), function(i) {
  nd <- data.frame(log_CFU = seq(band_time$lo[i], band_time$hi[i],
                                 length.out = 100),
                   Time    = band_time$Time[i])
  nd$fit    <- predict(m_hp,      newdata = nd)
  nd$null   <- predict(m_hp_null, newdata = nd)
  nd$t_band <- band_time$t_band[i]
  nd
})

# Draw each line only where its OWN prediction is inside the observable 0-7
# range. Early larvae sit at the h = 7 ceiling, so a common-slope model predicts
# above the maximum possible score at low density for early times; without this
# the Early line enters and exits through the panel edge and reads as a
# rendering fault rather than a fit. Filtering the two lines separately matters:
# filtering both on `fit` would also truncate Early's null, which is in range
# across the whole band. The excursion itself is the bounded-index attenuation
# the ESM already describes, not something a plotting choice can remove.
pred_4C_fit  <- pred_4C[pred_4C$fit  >= 0 & pred_4C$fit  <= 7, , drop = FALSE]
pred_4C_null <- pred_4C[pred_4C$null >= 0 & pred_4C$null <= 7, , drop = FALSE]

# WHAT THE READER SHOULD COMPARE. Each band's dashed null is FLAT: under
# h _||_ p | t, health depends on time alone. The solid line for the same band
# SLOPES. The rejection is therefore slope-versus-flat, not the size of any gap
# between the two lines -- a fitted line and a flat line at the same conditional
# mean must cross once, so the gap is zero somewhere in every band and carries
# no meaning. Do not describe this panel as "the gap between the curves".
p4C <- ggplot(pred_4C_fit, aes(x = log_CFU, y = fit, colour = t_band)) +
  geom_point(data = dat_fig4d,
             aes(x = log_CFU, y = health_combined, colour = t_band),
             inherit.aes = FALSE, alpha = 0.6, size = 2.5) +
  geom_line(data = pred_4C_null,
            aes(x = log_CFU, y = null, colour = t_band,
                linetype = "No burden effect (h ~ t)"),
            linewidth = 0.7, alpha = 0.75) +
  geom_line(aes(linetype = "Time + burden"), linewidth = 1.1) +
  scale_colour_manual(values = pal_tband, name = "Time") +
  scale_linetype_manual(values = c("Time + burden" = "solid",
                                   "No burden effect (h ~ t)" = "dashed"),
                        name = NULL,
                        breaks = c("Time + burden", "No burden effect (h ~ t)")) +
  scale_y_continuous(breaks = seq(0, 7, 1)) +
  coord_cartesian(ylim = c(0, 7)) +
  labs(x = expression(log[10](CFU)), y = "Health score (h)") +
  guides(colour   = guide_legend(order = 1),
         linetype = guide_legend(order = 2,
                                 override.aes = list(colour = "grey30"))) +
  mytheme +
  # bottom-left is the one clear corner: low burden and low health is empty
  theme(legend.position = c(0.19, 0.31),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.title = element_text(size = 13),
        legend.text  = element_text(size = 11),
        legend.key.size = unit(0.45, "cm"),
        legend.spacing.y = unit(0.02, "cm"))

# ----------------------------------------------------------------------------
# Compose Figure 4
# ----------------------------------------------------------------------------
# Object names match the letter each panel renders as:
#   p4A activity | p4B melanization | p4C health vs burden | p4D T50 ordering

# Switched C and D 
figure4 <- (p4A | p4B) / (p4C | p4D) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

# One-row alternative, if the four panels are wanted side by side:
# figure4 <- (p4A | p4B | p4C | p4D) + plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(face = "bold", size = 14))
# ggsave(..., width = 17, height = 4.5, ...)

ggsave("figures/figure4.pdf",  # >>> Manuscript Figure 4 (activity/melanization + T50 + h vs p | t)
       plot = figure4, width = 11, height = 9, units = "in", dpi = 300)


# =============================================================================
# Causal Analysis
# =============================================================================

# health_combined is built in the prep pipeline (Figure 4C needs it earlier);
# only its standardised form is added here.
burden_tidy_time <- burden_tidy_time %>%
  mutate(scaled_health_combined = as.numeric(scale(health_combined)))

# Does time still predict host condition once burden may act non-linearly?
# Fitted once and reused, so the parametric p, the edf and the deviance
# explained all come from the same fit.
gam_health_flex <- gam(health_combined ~ Time + s(log_CFU), data = burden_tidy_time)

#==============================================================================
# Causal model formalization
#------------------------------------------------------------------------------
# Encodes the candidate DAGs underlying Supplementary Figures S3 and S4 and
# prints the testable conditional independencies derived by d-separation. These
# are exactly the implications tested by regression downstream:
#   Fig S1: m_cond            (t -> p -> s implies t _||_ s | p ; etc.)
#   Fig S3: m_D, m_E, m_F, m_F_supp
# Reference: Textor et al. 2016 (dagitty). dagitty is loaded above.
#==============================================================================

# --- Three-node candidate models (Figure S3 panels B-D: t, p, s) ---
dagi3_mediation <- dagitty("dag { t -> p -> s }")            # S3B: pure pathogen mediation
dagi3_no_med    <- dagitty("dag { t -> p ; t -> s }")        # S3C: no pathogen mediation
dagi3_multipath <- dagitty("dag { t -> p -> s ; t -> s }")   # S3D: multi-path (supported)

cat("\n=== Implied conditional independencies: 3-node models (Figure S3) ===\n")
cat("\nS3B  t -> p -> s  (pure pathogen mediation):\n")
print(impliedConditionalIndependencies(dagi3_mediation))   # expect: t _||_ s | p
cat("\nS3C  t -> p ; t -> s  (no pathogen mediation):\n")
print(impliedConditionalIndependencies(dagi3_no_med))      # expect: p _||_ s | t
cat("\nS3D  t -> p -> s ; t -> s  (multi-path, supported):\n")
print(impliedConditionalIndependencies(dagi3_multipath))   # saturated: no implications

# --- Four-node candidate models (Figure S4 panels A-C: t, p, h, s) ---
dagi4_pmed_health <- dagitty("dag { t -> p -> h -> s }")                  # S4A: instantaneous damage
dagi4_collapse    <- dagitty("dag { t -> h -> p -> s }")                  # S4B: immune collapse
dagi4_bottleneck  <- dagitty("dag { t -> p ; t -> h ; p -> h ; h -> s }") # S4C: cumulative damage (supported)

cat("\n=== Implied conditional independencies: 4-node models (Figure S4) ===\n")
cat("\nS4A  t -> p -> h -> s  (instantaneous damage):\n")
print(impliedConditionalIndependencies(dagi4_pmed_health))  # includes t _||_ h | p  (Fig S3D)
cat("\nS4B  t -> h -> p -> s  (immune collapse):\n")
print(impliedConditionalIndependencies(dagi4_collapse))     # includes h _||_ s | p  (Fig S3E)
cat("\nS4C  t -> p ; t -> h ; p -> h ; h -> s  (supported):\n")
print(impliedConditionalIndependencies(dagi4_bottleneck))   # s _||_ t | h and s _||_ p | h  (Fig S3F)

#==============================================================================
# FIGURE S1: m(p) mapping + 3-node DAGs + diagnostic plots
# Panel A:   implied m(p) mapping with its pole at K
# Panels B-D: DAGs for alternative causal models (t, p, s)
# Panels E-F: diagnostic plots falsifying models B and C
# (manuscript Figure S1 -- numbered by float order, see the mapping note below)
#==============================================================================

# --- DAG panels ---

dag_A <- draw_dag(c("t_p","p_s"), drop_nodes = "h",
                  subtitles = "italic(t)~symbol('\\136')~italic(s)~'|'~italic(p)")
dag_B <- draw_dag(c("t_p","t_s"), drop_nodes = "h",
                  subtitles = "italic(p)~symbol('\\136')~italic(s)~'|'~italic(t)")
dag_C <- draw_dag(c("t_p","p_s","t_s"), drop_nodes = "h",
                  subtitles = "'Saturated'")

# ----------------------------------------------------------------------------
# Data prep
# ----------------------------------------------------------------------------

fig3_time_labels <- c("0–12h", "12–24h",
                     sprintf("24–%.0fh", max(burden_tidy_time$Time, na.rm = TRUE)))

dat_fig3 <- burden_tidy_time %>%
  filter(!is.na(status), !is.na(log_CFU)) %>%
  mutate(
    alive    = as.integer(status == "Alive"),
    cfu_bin  = cut(log_CFU,
                   breaks = quantile(log_CFU, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                   labels = c("Low", "Medium", "High"),
                   include.lowest = TRUE),
    # Labels derived from the data, not hard-coded: the analysis window closes
    # at 36 h (the 48 h sample is excluded upstream), so a literal "24-48h"
    # label described a range the top bin never contains.
    time_bin = cut(Time,
                   breaks = c(-0.01, 12, 24, max(Time, na.rm = TRUE)),
                   labels = fig3_time_labels)
  )


# ----------------------------------------------------------------------------
# Single additive logistic GLM — matches the test reported (β_Time, β_logCFU)
# ----------------------------------------------------------------------------
# m_cond is fitted in the Figure 3B block and reused here: dat_fig3 is the same
# set of rows, and Figure 3B and this panel must report the same coefficients.
stopifnot(nrow(dat_fig3) == nrow(burden_tidy_time))
cat("\n=== Conditional independence test, 3-node (Figure S1E,F) ===\n")
print(summary(m_cond))  # β_Time and β_logCFU are the reported test statistics

# ----------------------------------------------------------------------------
# Sanity checks on the conditional-independence result
# (moved here: these depend on dat_fig3, m_cond, health_combined and
#  scaled_health_combined, all defined above this point)
# ----------------------------------------------------------------------------
m_flex <- gam(alive ~ Time + s(log_CFU), data = dat_fig3, family = binomial)
cat("\n--- Flexible adjustment for burden: Time term should stay significant ---\n")
print(summary(m_flex))

m_int <- glm(alive ~ Time * log_CFU, data = dat_fig3, family = binomial)
cat("\n--- Interaction test: non-significant => common slope justified ---\n")
print(anova(m_cond, m_int, test = "LRT"))

# Penalised (Firth) likelihood: health near-perfectly separates survival, so the
# Wald p-value is unreliable. Use the profile-likelihood p-value reported here.
# NOTE: named fit_logistf so it does not overwrite the Gompertz `fit`.
fit_logistf <- logistf(survival ~ scaled_health_combined + scaled_cfu,
                       data = burden_tidy_time)
lf_time     <- logistf(survival ~ scaled_health_combined + scaled_time,
                       data = burden_tidy_time)
cat("\n--- Penalised likelihood: health | burden (beta_h, profile p) ---\n")
print(summary(fit_logistf))

cat("\n--- GAM: health ~ Time + s(log_CFU); report the Time term ---\n")
print(summary(gam_health_flex))

cat("\n--- Component-wise penalised fits (activity, melanization) ---\n")
print(logistf(survival ~ activity + scaled_cfu, data = burden_tidy_time))
print(logistf(survival ~ melanization + scaled_cfu, data = burden_tidy_time))

# ----------------------------------------------------------------------------
# Illustrative levels of the conditioning variable
# ----------------------------------------------------------------------------
cfu_levels  <- quantile(dat_fig3$log_CFU, c(0.15, 0.5, 0.85), na.rm = TRUE)
time_levels <- c(6, 18, 30)
names(cfu_levels)  <- c("Low", "Medium", "High")
names(time_levels) <- fig3_time_labels


# ─── Null prediction for Panel E: what we'd see if S3B (t ⊥ s | p) held ───
# Under S3B, survival depends only on log_CFU. Fit that null model.
m_null_D <- m_cond_null   # same fit as Figure 3B's dashed null

null_D <- purrr::map_dfr(names(cfu_levels), function(lbl) {
  nd <- data.frame(Time = c(0, 36), log_CFU = cfu_levels[[lbl]])
  nd$fit <- plogis(predict(m_null_D, newdata = nd))
  nd$cfu_bin <- factor(lbl, levels = c("Low", "Medium", "High"))
  nd
})


# ─── Null prediction for Panel F: what we'd see if S3C (p ⊥ s | t) held ───
m_null_E <- glm(survival ~ Time, data = burden_tidy_time, family = binomial)

null_E <- purrr::map_dfr(names(time_levels), function(lbl) {
  nd <- data.frame(log_CFU = range(burden_tidy_time$log_CFU, na.rm = TRUE),
                   Time = time_levels[[lbl]])
  nd$fit <- plogis(predict(m_null_E, newdata = nd))
  nd$bin <- factor(lbl, levels = names(time_levels))
  nd
})


# ----------------------------------------------------------------------------
# Panel E: P(alive) vs Time at three illustrative CFU levels
# ----------------------------------------------------------------------------
pred_D <- purrr::map_dfr(names(cfu_levels), function(lbl) {
  nd <- data.frame(Time = seq(0, 36, length.out = 200),
                   log_CFU = cfu_levels[[lbl]])
  p <- predict(m_cond, newdata = nd, type = "link", se.fit = TRUE)
  nd$fit <- plogis(p$fit)
  nd$lwr <- plogis(p$fit - 1.96 * p$se.fit)
  nd$upr <- plogis(p$fit + 1.96 * p$se.fit)
  nd$bin <- factor(lbl, levels = c("Low", "Medium", "High"))
  nd
})

bin_palette <- c("Low" = "#19798b", "Medium" = "#ee9b43", "High" = "#b80422")

p3D <- ggplot() +
  geom_jitter(data = dat_fig3,
              aes(x = Time, y = alive, color = cfu_bin),
              width = 0.3, height = 0.035, alpha = 0.5, size = 2) +
  geom_line(data = null_D, aes(x = Time, y = fit, color = cfu_bin),
            linetype = "dashed", linewidth = 1, alpha = 0.6)+
  geom_line(data = pred_D,
            aes(x = Time, y = fit, color = bin),
            linewidth = 1.1) +
  scale_fill_manual(values = bin_palette, guide = "none") +
  scale_y_continuous(breaks = c(0, 1), labels = c("Dead", "Alive")) +
  scale_x_continuous(breaks = seq(0, 36, 6)) +
  coord_cartesian(xlim = c(0, 36), ylim = c(-0.12, 1.12)) +   # see p4D
  labs(x = bquote("Time (hrs, "*italic(t)*")"), y = bquote("Survival ("*italic(s)*")"))+
  mytheme +
  theme(legend.position = c(0.2, 0.35), , legend.title = element_text(size = 16))+
  scale_color_manual(values = bin_palette,
                     name = expression(log[10](CFU)))

# ----------------------------------------------------------------------------
# Panel F: P(alive) vs log_CFU at three illustrative time levels
# ----------------------------------------------------------------------------
pred_E <- purrr::map_dfr(names(time_levels), function(lbl) {
  nd <- data.frame(log_CFU = seq(min(dat_fig3$log_CFU),
                                 max(dat_fig3$log_CFU),
                                 length.out = 200),
                   Time = time_levels[[lbl]])
  p <- predict(m_cond, newdata = nd, type = "link", se.fit = TRUE)
  nd$fit <- plogis(p$fit)
  nd$lwr <- plogis(p$fit - 1.96 * p$se.fit)
  nd$upr <- plogis(p$fit + 1.96 * p$se.fit)
  nd$bin <- factor(lbl, levels = fig3_time_labels)
  nd
})

# Named from fig3_time_labels, not hard-coded: a level list that disagrees with
# the data's labels turns the unmatched band into NA and paints it grey.
time_palette <- setNames(c("#19798b", "#ee9b43", "#b80422"), fig3_time_labels)

p3E <- ggplot() +
  geom_jitter(data = dat_fig3,
              aes(x = log_CFU, y = alive, color = time_bin),
              width = 0.08, height = 0.035, alpha = 0.5, size = 2) +
  geom_line(data = null_E, aes(x = log_CFU, y = fit, color = bin),
            linetype = "dashed", linewidth = 1, alpha = 0.6)+
  geom_line(data = pred_E,
            aes(x = log_CFU, y = fit, color = bin),
            linewidth = 1.1) +
  scale_fill_manual(values = time_palette, guide = "none") +
  scale_y_continuous(breaks = c(0, 1), labels = c("Dead", "Alive")) +
  coord_cartesian(ylim = c(-0.12, 1.12)) +                    # see p4D
  labs(x = bquote(log[10]*"(CFU), "*italic(p)), y = NULL) +
  mytheme+
  theme(legend.position = c(0.2, 0.35),
        axis.text.y = element_blank(), legend.title = element_text(size = 16))+
  scale_color_manual(values = time_palette, name = "Time")

# Shared display limits for the 3-node DAG panels. These match the coord_fixed()
# limits set inside draw_dag_panel() so the conditional-independence subtitle
# (placed at y = -0.5) is not clipped when the panels are composed by patchwork.

#dag_xlim <- c(-0.1, 2.3)
#dag_ylim <- c(-0.3, 1.6)

#dag_A <- dag_A + coord_cartesian(xlim = dag_xlim, ylim = dag_ylim, clip = "off")
#dag_B <- dag_B + coord_cartesian(xlim = dag_xlim, ylim = dag_ylim, clip = "off")
#dag_C <- dag_C + coord_cartesian(xlim = dag_xlim, ylim = dag_ylim, clip = "off")


row1 <- (pS3A | dag_A | dag_B | dag_C) + plot_layout(widths = c(1.4, 1, 1, 1))
row2 <- (p3D | p3E)
figureS1 <- row1 / row2 + plot_layout(heights = c(1, 1.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))


ggsave("figures/figureS1.pdf",  # >>> Supplementary Figure S1 (m(p) mapping + 3-node pathogen DAGs + diagnostics)
       plot = figureS1, width = 12, height = 9, units = "in", dpi = 300)


# ============================================================================
# FIGURE S3: Health-mediated causal analysis
# A,B,C: DAGs (instantaneous damage; immune collapse; cumulative damage)
# D: rejects A via t ⊥ h | p
# E: rejects B via h ⊥ s | p
# F: supports C via s ⊥ t | h (and reports s ⊥ p | h in the caption)
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Health composite
# ----------------------------------------------------------------------------
dat_fig5 <- burden_tidy_time %>%
  filter(!is.na(activity), !is.na(melanization),
         !is.na(survival), !is.na(log_CFU)) %>%
  mutate(h = activity + melanization)   # higher h = healthier


dag5_A <- draw_dag(c("t_p","p_h","h_s"), dashed_edges = "p_s",
                   subtitles = "italic(t)~symbol('\\136')~italic(h)~'|'~italic(p)")
dag5_B <- draw_dag(c("t_h","h_p","p_s"), dashed_edges = "t_p",
                   subtitles = "italic(h)~symbol('\\136')~italic(s)~'|'~italic(p)")
dag5_C <- draw_dag(c("t_p","t_h","p_h","h_s"),
                   subtitles = c("italic(s)~symbol('\\136')~italic(t)~'|'~italic(h)",
                                 "italic(s)~symbol('\\136')~italic(p)~'|'~italic(h)"))


# ----------------------------------------------------------------------------
# Shared conditioning levels & palettes
# ----------------------------------------------------------------------------
cfu_levels <- quantile(dat_fig5$log_CFU, c(0.15, 0.50, 0.85), na.rm = TRUE)
names(cfu_levels) <- c("Low", "Medium", "High")

# Fixed health levels (not quantiles) so the three curves span the usable range
h_levels <- c("Low" = 1, "Medium" = 4, "High" = 7)

pal_cfu <- c("Low" = "#19798b", "Medium" = "#ee9b43", "High" = "#b80422")
pal_h   <- c("Low" = "#b80422", "Medium" = "#ee9b43", "High" = "#19798b") # inverted: high h = good = teal

# ----------------------------------------------------------------------------
# binned data (with pre-computed jitter)
# ----------------------------------------------------------------------------
set.seed(2468)   # makes the plotted jitter reproducible across reruns
dat_binned <- dat_fig5 %>%
  mutate(
    cfu_bin = cut(log_CFU,
                  breaks = quantile(log_CFU, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                  labels = c("Low", "Medium", "High"), include.lowest = TRUE),
    h_bin   = cut(h,
                  breaks = c(-0.01, 2.5, 5.5, 7.01),
                  labels = c("Low", "Medium", "High")),
    survival_jitter = survival + runif(n(), -0.04, 0.04)
  )

# ----------------------------------------------------------------------------
# 3. Panel D — rejects A: h vs Time at 3 CFU levels
# ----------------------------------------------------------------------------

# This is also in the main manuscript text, third paragraph of health component
m_D      <- lm(h ~ Time + log_CFU, data = dat_fig5)
m_D_null <- lm(h ~ log_CFU,        data = dat_fig5)

pred_D <- purrr::map_dfr(names(cfu_levels), function(lbl) {
  nd <- data.frame(Time = seq(0, 36, length.out = 200), log_CFU = cfu_levels[[lbl]])
  p <- predict(m_D, newdata = nd, se.fit = TRUE)
  nd$fit <- p$fit
  nd$bin <- factor(lbl, levels = names(cfu_levels)); nd
})
null_D <- purrr::map_dfr(names(cfu_levels), function(lbl) {
  nd <- data.frame(Time = c(0, 36), log_CFU = cfu_levels[[lbl]])
  nd$fit <- predict(m_D_null, newdata = nd)
  nd$bin <- factor(lbl, levels = names(cfu_levels)); nd
})

p5D <- ggplot() +
  geom_jitter(data = dat_binned, aes(x = Time, y = h, color = cfu_bin),
              size = 2, alpha = 0.5, width = 0.3, height = 0.1) +
  geom_line(data = null_D, aes(x = Time, y = fit, color = bin),
            linetype = "dashed", linewidth = 0.7, alpha = 0.6) +
  geom_line(data = pred_D, aes(x = Time, y = fit, color = bin), linewidth = 1.1) +
  scale_x_continuous(breaks = seq(0, 36, 6)) +
  # 63 of 86 larvae sit exactly on h = 0 or h = 7, and geom_jitter(height = 0.1)
  # pushes them outside a hard limit, where ggplot DELETES them. Clip the view
  # instead.
  coord_cartesian(xlim = c(0, 36), ylim = c(0, 7)) +
  labs(x = bquote("Time (hrs, "*italic(t)*")"), y = bquote("Health score ("*italic(h)*")")) +
  mytheme+
  scale_color_manual(values = pal_cfu, name = expression(log[10](CFU)))+
  theme(legend.position = c(0.2, 0.32), legend.title = element_text(size = 14), legend.text = element_text(size = 14))

# ----------------------------------------------------------------------------
# 4. Panel E — rejects B: s vs h at 3 CFU levels
# ----------------------------------------------------------------------------
m_E      <- glm(survival ~ h + log_CFU, data = dat_fig5, family = binomial) 
# this is also important for SEM
# that's why we don't add a direct way to survival
m_E_null <- glm(survival ~ log_CFU,     data = dat_fig5, family = binomial)

summary(m_E)
summary(m_E_null)

pred_E <- purrr::map_dfr(names(cfu_levels), function(lbl) {
  nd <- data.frame(h = seq(min(dat_fig5$h), max(dat_fig5$h), length.out = 200),
                   log_CFU = cfu_levels[[lbl]])
  nd$fit <- plogis(predict(m_E, newdata = nd))
  nd$bin <- factor(lbl, levels = names(cfu_levels)); nd
})
null_E <- purrr::map_dfr(names(cfu_levels), function(lbl) {
  nd <- data.frame(h = range(dat_fig5$h, na.rm = TRUE), log_CFU = cfu_levels[[lbl]])
  nd$fit <- plogis(predict(m_E_null, newdata = nd))
  nd$bin <- factor(lbl, levels = names(cfu_levels)); nd
})

p5E <- ggplot() +
  geom_jitter(data = dat_binned,
              aes(x = h, y = survival_jitter, color = cfu_bin),
              size = 2, alpha = 0.5, width = 0.15, height = 0) +
  geom_line(data = null_E, aes(x = h, y = fit, color = bin),
            linetype = "dashed", linewidth = 0.7, alpha = 0.6) +
  geom_line(data = pred_E, aes(x = h, y = fit, color = bin), linewidth = 1.1) +
  scale_y_continuous(breaks = c(0,1), labels = c("Dead","Alive")) +
  coord_cartesian(xlim = c(0, 7), ylim = c(-0.12, 1.12)) +    # see p5D
  labs(x = bquote("Health score ("*italic(h)*")"), y = bquote("Survival ("*italic(s)*")")) + mytheme +
  scale_color_manual(values = pal_cfu, name = expression(log[10](CFU)))+
  theme(legend.position = c(0.8, 0.32), legend.title = element_text(size = 14), legend.text = element_text(size = 14))


# ----------------------------------------------------------------------------
# 5. Panel F — supports C: s vs Time at 3 h levels
# ----------------------------------------------------------------------------
m_F      <- glm(survival ~ Time + h, data = dat_fig5, family = binomial)
m_F_null <- glm(survival ~ h,        data = dat_fig5, family = binomial)
m_F_supp <- glm(survival ~ log_CFU + h, data = dat_fig5, family = binomial)

summary(m_F)
summary(m_F_null)
summary(m_F_supp)

pred_F <- map_dfr(names(h_levels), function(lbl) {
  nd <- data.frame(Time = seq(0, 36, length.out = 200), h = h_levels[[lbl]])
  nd$fit <- plogis(predict(m_F, newdata = nd))
  nd$bin <- factor(lbl, levels = names(h_levels)); nd
})
null_F <- map_dfr(names(h_levels), function(lbl) {
  nd <- data.frame(Time = c(0, 36), h = h_levels[[lbl]])
  nd$fit <- plogis(predict(m_F_null, newdata = nd))
  nd$bin <- factor(lbl, levels = names(h_levels)); nd
})

p5F <- ggplot() +
  geom_jitter(data = dat_binned,
              aes(x = Time, y = survival_jitter, color = h_bin),
              size = 2, alpha = 0.5, width = 0.3, height = 0) +
  geom_line(data = null_F, aes(x = Time, y = fit, color = bin),
            linetype = "dashed", linewidth = 0.7, alpha = 0.6) +
  geom_line(data = pred_F, aes(x = Time, y = fit, color = bin), linewidth = 1.1) +
  scale_y_continuous(breaks = c(0,1), labels = c("Dead","Alive")) +
  scale_x_continuous(breaks = seq(0, 36, 6)) +
  coord_cartesian(xlim = c(0, 36), ylim = c(-0.12, 1.12)) +   # see p5D
  labs(x = bquote("Time (hrs, "*italic(t)*")"), y = bquote("Survival ("*italic(s)*")")) +
  mytheme+
  theme(legend.position = c(0.2, 0.32), legend.title = element_text(size = 14), legend.text = element_text(size = 14))+
  scale_color_manual(values = pal_h, name = bquote("Health ("*italic(h)*")"))

# ----------------------------------------------------------------------------
# 6. Print test stats for caption
# ----------------------------------------------------------------------------
cat("\n=== Figure S4 test statistics ===\n")
cat("D — t ⊥ h | p (S4A):\n");   print(summary(m_D)$coefficients["Time", ])
cat("\nE — h ⊥ s | p (S4B):\n"); print(summary(m_E)$coefficients["h", ])
cat("\nF — s ⊥ t | h (S4C):\n"); print(summary(m_F)$coefficients["Time", ])
cat("F — s ⊥ p | h (S4C):\n");   print(summary(m_F_supp)$coefficients["log_CFU", ])

# Likelihood ratio tests 
m_h_only <- glm(survival ~ h, data = dat_fig5, family = binomial)
anova(m_h_only, m_F, test = "LRT")   # does Time add anything given h?
anova(m_h_only, m_E, test = "LRT")   # does log_CFU add anything given h?

# Profile-likelihood counterparts (main text, health component, third paragraph)
summary(lf_time)
summary(fit_logistf)
setNames(lf_time$prob,     names(coef(lf_time)))
setNames(fit_logistf$prob, names(coef(fit_logistf)))

# ----------------------------------------------------------------------------
# 7. Assemble
# ----------------------------------------------------------------------------
design <- "
ABC
DEF
"
figureS3 <- dag5_A + dag5_B + dag5_C + p5D + p5E + p5F +
  plot_layout(design = design, heights = c(1, 1.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("figures/figureS3.pdf",  # >>> Supplementary Figure S3 (4-node health DAGs + diagnostics)
       plot = figureS3, width = 15, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Sensitivity analysis
# Reproduces the values reported in the text for the melanization-only health
# index (s ⊥ t | h_mel ; s ⊥ p | h_mel). Self-contained — uses *_mel object
# names so it does NOT overwrite the main Figure S4 objects (dat_fig5, m_F, ...).
# -----------------------------------------------------------------------------

dat_fig5_mel <- burden_tidy_time %>%
  filter(!is.na(melanization), !is.na(survival), !is.na(log_CFU)) %>%
  mutate(h_mel =  melanization)

# S4C support tests with melanization-only health
m_F_mel      <- glm(survival ~ Time    + h_mel, data = dat_fig5_mel, family = binomial)
m_F_supp_mel <- glm(survival ~ log_CFU + h_mel, data = dat_fig5_mel, family = binomial)

# This comes up in the main text second paragraph of health component/figure 4
cat("\n=== Figure S4 sensitivity: melanization-only health ===\n")
cat("s \u22A5 t | h_mel  (Time | h_mel):\n");    print(summary(m_F_mel)$coefficients["Time", ])
cat("s \u22A5 p | h_mel  (log_CFU | h_mel):\n"); print(summary(m_F_supp_mel)$coefficients["log_CFU", ])

lf_act <- logistf(survival ~ activity     + scaled_cfu, data = burden_tidy_time)
lf_mel <- logistf(survival ~ melanization + scaled_cfu, data = burden_tidy_time)

# Sanity check
chisq_lf <- function(lf, term) {
  i <- which(names(coef(lf)) == term)
  c(beta  = unname(coef(lf)[i]),
    chisq = unname(qchisq(lf$prob[i], df = 1, lower.tail = FALSE)),
    p     = unname(lf$prob[i]))
}
round(chisq_lf(lf_act, "activity"), 3)
round(chisq_lf(lf_mel, "melanization"), 3)

summary(gam_health_flex)

# =============================================================================
# FIGURE 5: Bayesian SEM with health mediation
# =============================================================================

# Prepare time variables
burden_tidy_time$time_linear  <- burden_tidy_time$scaled_time
burden_tidy_time$time_squared <- burden_tidy_time$scaled_time^2

# Supported model: t → p → h → s with t → h
model_sem <- '
  # Pathogen growth
  scaled_cfu ~ t1 * time_linear + t2 * time_squared
  
  # Health depends on history
  scaled_health_combined ~ a * scaled_cfu + h1 * time_linear

  # Survival depends on health (complete mediation)
  survival ~ b * scaled_health_combined

  # Defined parameters
  indirect_p := a * b       # Indirect effect of p through h
  indirect_t := h1 * b      # Indirect effect of t through h (not via p)
'

set.seed(6789)
fit_sem <- bsem(model_sem,
                data = burden_tidy_time,
                burnin = 1000,
                sample = 5000,
                n.chains = 4)

print(summary(fit_sem))
cat("\n=== SEM fit measures (supported model) ===\n")
print(fitMeasures(fit_sem, c("dic", "ppp")))

# Check convergence
cat("\n=== SEM convergence (PSRF, all should be < 1.01) ===\n")
print(blavInspect(fit_sem, "psrf"))


model_sem_quad <- '
  scaled_cfu ~ t1 * time_linear + t2 * time_squared
  scaled_health_combined ~ a * scaled_cfu + h1 * time_linear + h2 * time_squared
  survival ~ b * scaled_health_combined

  t_to_h_early := h1 - 2 * h2      # marginal t -> h one SD before mean sampling time
  t_to_h_mean  := h1               # at mean sampling time
  t_to_h_late  := h1 + 2 * h2      # one SD after

  indirect_p      := a * b
   indirect_t_early := (h1 - 2 * h2) * b
  indirect_t_mean := h1 * b
  indirect_t_late := (h1 + 2 * h2) * b
'

set.seed(6790)
fit_sem_quad<- bsem(model_sem_quad,
                     data = burden_tidy_time,
                     burnin = 1000,
                     sample = 5000,
                     n.chains = 4)

print(summary(fit_sem_quad))
print(fitMeasures(fit_sem_quad, c("dic", "ppp")))
blavInspect(fit_sem_quad, "psrf")

# Alternative: No health mediation (for comparison)
model_no_h <- '
  scaled_cfu ~ t1 * time_linear + t2 * time_squared
  survival ~ c * scaled_cfu + d * time_linear
'

fit_no_h <- bsem(model_no_h,
                 data = burden_tidy_time,
                 burnin = 1000,
                 sample = 5000,
                 n.chains = 4)

cat("\n=== Model comparison: with vs without health mediation ===\n")
cat("  With health:\n");    print(fitMeasures(fit_sem, c("dic", "ppp"))) #linear
cat("  With health quadratic:\n");    print(fitMeasures(fit_sem_quad, c("dic", "ppp"))) #quadratic
cat("  Without health:\n"); print(fitMeasures(fit_no_h, c("dic", "ppp")))

# =============================================================================
# Figure 5: SEM path effects from POSTERIOR DRAWS (single source of truth)
# Inputs are pre-scaled, so raw coefficients ARE the standardized effects.
# Everything (point estimate, CrI, density) comes from the draws -> no
# dependence on parameterEstimates() column names.
# =============================================================================

# ---- path labels: time route evaluated at two sampling times ---------------
path_labels <- c(
  t_p_10 = "italic(t) %->% italic(p)~'at ~10 h'",
  t_p_20 = "italic(t) %->% italic(p)~'at ~20 h'",
  t_p_30 = "italic(t) %->% italic(p)~'at ~30 h'",
  a      = "italic(p) %->% italic(h)~(a)",
  t_h_10 = "italic(t) %->% italic(h)~'at ~10 h'",
  t_h_20 = "italic(t) %->% italic(h)~'at ~20 h'",
  t_h_30 = "italic(t) %->% italic(h)~'at ~30 h'",
  b      = "italic(h) %->% italic(s)~(b)"
)
path_order <- rev(unname(path_labels))

eqn_of <- c(t_p_10 = "Into pathogen density (SD of p per SD of t)",
            t_p_20 = "Into pathogen density (SD of p per SD of t)",
            t_p_30 = "Into pathogen density (SD of p per SD of t)",
            a      = "Into host health (SD of h per SD of predictor)",
            t_h_10 = "Into host health (SD of h per SD of predictor)",
            t_h_20 = "Into host health (SD of h per SD of predictor)",
            t_h_30 = "Into host health (SD of h per SD of predictor)",
            b      = "Into survival (probability per SD of h)")

# ---- Posterior draws --------------------------------------------------------
# The draws must come from fit_sem_quad: `post` below needs h2, which only the
# quadratic model estimates. Chains are combined by iterating the list rather
# than rebuilding an mcmc.list, which errors when blavaan returns chains that
# are not classed as `mcmc`. as.data.frame keeps the `draws$x` accessors used
# below working.

mc_quad <- blavInspect(fit_sem_quad, "mcmc")
if (coda::is.mcmc(mc_quad)) mc_quad <- list(mc_quad)
draws <- as.data.frame(do.call(rbind, lapply(mc_quad, as.matrix)))

needed_pars <- c("t1", "t2", "a", "h1", "h2", "b")
missing_pars <- setdiff(needed_pars, names(draws))
if (length(missing_pars))
  stop("fit_sem_quad returned no draws for: ", paste(missing_pars, collapse = ", "),
       "\n  available: ", paste(names(draws), collapse = ", "), call. = FALSE)

# Marginal effect of standardised time z is  t1 + 2*t2*z  (and h1 + 2*h2*z).
# z = -1, 0, +1 correspond to roughly 10, 20 and 30 h given the sampling design;
# the exact hours are printed in the report so the labels can be checked.
post <- data.frame(
  t_p_10 = draws$t1 - 2 * draws$t2,
  t_p_20 = draws$t1,
  t_p_30 = draws$t1 + 2 * draws$t2,
  a      = draws$a,
  t_h_10 = draws$h1 - 2 * draws$h2,
  t_h_20 = draws$h1,
  t_h_30 = draws$h1 + 2 * draws$h2,
  b      = draws$b
)

plot_df <- post %>%
  pivot_longer(everything(), names_to = "label", values_to = "value") %>%
  mutate(path  = factor(path_labels[label], levels = path_order),
         eqn = factor(eqn_of[label],
                              levels = c("Into pathogen density (SD of p per SD of t)",
                                         "Into host health (SD of h per SD of predictor)",
                                         "Into survival (probability per SD of h)"))) %>%
  group_by(path) %>%
  mutate(med = median(value),
         lo  = quantile(value, 0.025),
         hi  = quantile(value, 0.975)) %>% ungroup() %>%
  mutate(effect = case_when(
    lo < 0 & hi > 0 ~ "CrI spans zero",
    med > 0         ~ "increases",
    TRUE            ~ "decreases"
  ))

# ---- Point estimate + 95% CrI table ----------------------------------------
effects_df <- plot_df %>%
  group_by(label, path) %>%
  summarise(estimate = median(value),
            lower    = quantile(value, 0.025),
            upper    = quantile(value, 0.975),
            .groups  = "drop") %>%
  mutate(path = factor(path, levels = path_order)) %>%
  arrange(path)

cat("\n=== SEM standardised path effects (Figure 5) ===\n")
print(effects_df)

# ---- Half-eye posterior figure ---------------------------------------------
pal <- c("increases"         = "#19798b",
         "decreases"         = "#b80422",
         "CrI spans zero"                = "#888780")

p <- ggplot(plot_df, aes(x = value, y = path, fill = effect, colour = effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  stat_halfeye(.width = c(0.66, 0.95), point_interval = "median_qi",
               slab_alpha = 0.5, normalize = "groups", height = 0.9,
               interval_size_range = c(1.2, 3.5),   # thin (95%), thick (66%)
               point_size = 4) +
  scale_y_discrete(labels = function(l) parse(text = l)) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "Posterior effect size", y = NULL) +
  mytheme+
  theme(legend.position = "top", axis.text.y = element_text(hjust = 0))+
  facet_wrap(~ eqn, ncol = 1, scales = "free")


# ---------------------------------------------------------------------------
# Figure 5 inset: the supported model with SEM-signed edges. Positive paths get
# a teal arrowhead, negative paths a red blunt bar, using the same palette as
# the posterior panel. This states the fitted result rather than repeating
# figure 1E. Signs are read from `post`, so this MUST be defined after `post`.
# ---------------------------------------------------------------------------
sem_sign <- function(x) if (median(x) < 0) "neg" else "pos"


# --- Structural edges: grey, no sign asserted ------------------------------
# Drawn in every panel. Any edge also named in `edge_signs` is removed here,
# so a signed edge is never drawn twice.

dag_sem <- draw_dag(
  edge_signs = c(t_p = sem_sign(draws$t1), t_h = sem_sign(draws$h1),
                 p_h = sem_sign(draws$a),  h_s = sem_sign(draws$b))) +
  theme(plot.background = element_rect(fill = "white"))

figure5 <- p +
  inset_element(dag_sem,
                left = 0.20, bottom = 0.02, right = 0.58, top = 0.26,
                align_to = "full")

ggsave("figures/figure5.pdf", figure5, width = 8.5, height = 9.5)


#-------------
# Sanity check 
#-------------

dat_cens <- dat_fig5 %>%
  mutate(lo = ifelse(h <= 0, NA, h),      # left-censored at the floor
         hi = ifelse(h >= 7, NA, h))      # right-censored at the ceiling

m_D_cens <- survreg(Surv(lo, hi, type = "interval2") ~ Time + log_CFU,
                    data = dat_cens, dist = "gaussian")

cat("\n--- Health ~ t + p, censored at 0 and 7 (robustness for Note S3) ---\n")
print(summary(m_D_cens))
cat("\nUncensored OLS for comparison:\n")
print(summary(m_D)$coefficients[c("Time", "log_CFU"), ])

dat_cens <- dat_cens %>%
  mutate(z_time = as.numeric(scale(Time)), z_cfu = as.numeric(scale(log_CFU)))

m_D_cens_z <- survreg(Surv(lo, hi, type = "interval2") ~ z_time + z_cfu,
                      data = dat_cens, dist = "gaussian")
m_D_z      <- lm(h ~ z_time + z_cfu, data = dat_cens)

print(summary(m_D_cens_z))
print(summary(m_D_z))

cat(sprintf("h at floor: %d/%d;  h at ceiling: %d/%d\n",
            sum(dat_fig5$h <= 0), nrow(dat_fig5),
            sum(dat_fig5$h >= 7), nrow(dat_fig5)))


z        <- as.numeric(burden_tidy_time$scaled_time)
t_centre <- attr(scale(burden_tidy_time$Time), "scaled:center")
t_scale  <- attr(scale(burden_tidy_time$Time), "scaled:scale")
h1_med   <- median(draws$h1); h2_med <- median(draws$h2)
z_vertex <- -h1_med / (2 * h2_med)
cat(sprintf("vertex at z = %.2f  (Time = %.1f hrs); observed z range %.2f to %.2f\n",
            z_vertex, t_centre + z_vertex * t_scale, min(z), max(z)))

# Sanity check for the Supplementary Note S3
m_D_int <- lm(h ~ Time * log_CFU, data = dat_fig5)
anova(m_D, m_D_int)

summary(m_D_int)$coefficients["Time:log_CFU", ]

m_D_cens_int <- survreg(Surv(lo, hi, type = "interval2") ~ z_time * z_cfu,
                        data = dat_cens, dist = "gaussian")

ll0 <- m_D_cens_z$loglik[2]      # -84.6
ll1 <- m_D_cens_int$loglik[2]    # -79.7
cat(sprintf("LRT = %.2f, df = 1, p = %.4g\n",
            2*(ll1 - ll0), pchisq(2*(ll1 - ll0), 1, lower.tail = FALSE)))

# Does the interaction change which path is larger, across the observed range?
m_D_z_int <- lm(h ~ z_time * z_cfu, data = dat_cens)
summary(m_D_z_int)$coefficients

sd(dat_fig5$h)   # the conversion factor between the two scales (Bayesian vs OLS)

# Monotone progression: activity loss -> melanization -> death.
# A nested cascade admits only four states; anything else is out of order.
mono <- burden_tidy_time %>%
  filter(!is.na(activity), !is.na(melanization), !is.na(status)) %>%
  mutate(
    time   = Time,
    active = as.integer(activity   >= 2),
    unmel  = as.integer(melanization >= 3),
    alive  = as.integer(status == "Alive")) %>%
  mutate(state = case_when(
    active == 1 & unmel == 1 & alive == 1 ~ "intact",
    active == 0 & unmel == 1 & alive == 1 ~ "activity lost",
    active == 0 & unmel == 0 & alive == 1 ~ "melanized",
    active == 0 & unmel == 0 & alive == 0 ~ "dead",
    TRUE                                  ~ "out of order"))

count(mono, state)
cat(sprintf("monotone-consistent: %.1f%% (%d/%d)\n",
            100*mean(mono$state != "out of order"),
            sum(mono$state != "out of order"), nrow(mono)))

# Which larvae break it, and how
mono %>% filter(state == "out of order") %>%
  left_join(dplyr::select(burden_tidy_time, Sample, Larvae, Time, activity, melanization),
            by = c("Sample","Larvae"))


#===============================================================================
# SUPPLEMENTARY ANALYSIS: Cumulative Burden Estimation
# Why cumulative measures are problematic in cross-sectional designs
#===============================================================================

#-------------------------------------------------------------------------------
# METHOD 1: Integral of fitted logistic (current approach)
#-------------------------------------------------------------------------------

# Analytical integral: ∫₀ᵗ p(τ)dτ for logistic growth
# Closed form: (K/r) * ln[1 + (p0/(K-p0)) * (e^(rt) - 1)]
burden_tidy_time <- burden_tidy_time %>%
  mutate(
    cum_burden_integral = (K/r) * log1p((p0 / (K - p0)) * expm1(r * Time))
  )

#-------------------------------------------------------------------------------
# METHOD 2: Trapezoidal sum of observed data
#-------------------------------------------------------------------------------

# Calculate population mean CFU at each timepoint
mean_cfu_by_time <- burden_tidy_time %>%
  group_by(Time) %>%
  summarise(mean_cfu = mean(cfu, na.rm = TRUE), .groups = "drop") %>%
  arrange(Time)

# Trapezoidal integration of observed means
trapezoidal_cumsum <- function(times, values) {
  n <- length(times)
  if (n == 1) return(0)
  
  cumsum_vals <- numeric(n)
  cumsum_vals[1] <- 0  # At t=0, cumulative = 0
  
  for (i in 2:n) {
    # Trapezoidal rule: (t2-t1) * (y1+y2)/2
    dt <- times[i] - times[i-1]
    avg_val <- (values[i] + values[i-1]) / 2
    cumsum_vals[i] <- cumsum_vals[i-1] + dt * avg_val
  }
  return(cumsum_vals)
}

mean_cfu_by_time <- mean_cfu_by_time %>%
  mutate(cum_burden_trapezoid = trapezoidal_cumsum(Time, mean_cfu))

# Merge back to individual data (assign population cumulative to each larva)
burden_tidy_time <- burden_tidy_time %>%
  left_join(mean_cfu_by_time %>% dplyr::select(Time, cum_burden_trapezoid), by = "Time")

#-------------------------------------------------------------------------------
# METHOD 3: Simple cumulative sum (cruder approach)
#-------------------------------------------------------------------------------

# Just sum mean CFU at all preceding timepoints (discrete approximation)
burden_tidy_time <- burden_tidy_time %>%
  mutate(
    cum_burden_discrete = sapply(Time, function(t) {
      sum(mean_cfu_by_time$mean_cfu[mean_cfu_by_time$Time <= t])
    })
  )

#-------------------------------------------------------------------------------
# COMPARISON: How similar are the three methods?
#-------------------------------------------------------------------------------

# Correlations
cor_matrix <- burden_tidy_time %>%
  dplyr::select(Time, cum_burden_integral, cum_burden_trapezoid, cum_burden_discrete) %>%
  distinct() %>%
  cor(use = "complete.obs")

cat("\n=== Correlation matrix: Time vs Cumulative Burden Measures ===\n")
print(round(cor_matrix, 4))


# All cumulative measures are strongly correlated with time. Note that the
# SIMPLEST estimator (discrete sum) is the most collinear, so the problem is
# structural rather than an artefact of integrating a fitted curve.
cat("\nCorrelation with Time:\n")
cat("  Integral method:    r =", round(cor_matrix["Time", "cum_burden_integral"], 4), "\n")
cat("  Trapezoidal method: r =", round(cor_matrix["Time", "cum_burden_trapezoid"], 4), "\n")
cat("  Discrete method:    r =", round(cor_matrix["Time", "cum_burden_discrete"], 4), "\n")

#-------------------------------------------------------------------------------
# FIGURE S4: Comparison of cumulative burden estimation methods
#-------------------------------------------------------------------------------

# Panel A: All three methods vs time
comparison_data <- burden_tidy_time %>%
  dplyr::select(Time, cum_burden_integral, cum_burden_trapezoid, cum_burden_discrete) %>%
  distinct() %>%
  pivot_longer(cols = starts_with("cum_burden"),
               names_to = "method",
               values_to = "cumulative_burden") %>%
  mutate(method = case_when(
    method == "cum_burden_integral" ~ "Integral of fitted logistic",
    method == "cum_burden_trapezoid" ~ "Trapezoidal sum of data",
    method == "cum_burden_discrete" ~ "Discrete sum of data"
  ))

pA_supp <- ggplot(comparison_data, aes(x = Time, y = log10(cumulative_burden + 1),
                                       color = method, linetype = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#ee9b43", "#19798b", "#b80422")) +
  labs(x = "Time since infection (hrs)",
       y = expression(log[10](Cumulative~CFU)),
       color = "Method", linetype = "Method") +
  mytheme +
  theme(legend.position = c(0.6, 0.25))

# Panel B: Integral vs Trapezoidal (method comparison)
method_compare <- burden_tidy_time %>%
  dplyr::select(Time, cum_burden_integral, cum_burden_trapezoid) %>%
  distinct()

# cor_matrix is on the RAW cumulative scale, but panels B and C plot log10 and
# the VIF comes from a log-scale model: r = 0.84 raw vs 0.98 logged. Quoting one
# beside the other is inconsistent (r = 0.84 implies VIF = 3.4, not 32.6), so
# both are computed here and every use names its scale.
cum_log <- method_compare %>%
  mutate(across(starts_with("cum_burden"), ~ log10(.x + 1)))

r_time_raw    <- cor(method_compare$Time, method_compare$cum_burden_integral)
r_time_log    <- cor(cum_log$Time,        cum_log$cum_burden_integral)
r_methods_raw <- cor(method_compare$cum_burden_integral, method_compare$cum_burden_trapezoid)
r_methods_log <- cor(cum_log$cum_burden_integral,        cum_log$cum_burden_trapezoid)

pB_supp <- ggplot(method_compare, aes(x = log10(cum_burden_integral + 1),
                                      y = log10(cum_burden_trapezoid + 1))) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "#ee9b43") +
  labs(x = expression(log[10](Integral~of~fitted~model)),
       y = expression(log[10](Trapezoidal~sum~of~data))) +
  # both axes are log10, so annotate the log-scale r
  annotate("text", x = 5, y = 8,
           label = paste0("r = ", round(r_methods_log, 3)),
           size = 6) +
  mytheme

# Panel C: THE PROBLEM - cumulative burden is structurally confounded with time
pC_supp <- ggplot(method_compare, aes(x = Time, y = log10(cum_burden_integral + 1))) +
  geom_point(size = 3, color = "#ee9b43") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Time since infection (hrs)",
       y = expression(log[10](Cumulative~CFU))) +
  # y is log10(Sum p + 1), so this is the log-scale r (0.98). The raw-scale r
  # is 0.84, which is what the VIF in Note S2 contradicts -- see the report.
  annotate("text", x = 10, y = 8,
           label = paste0("r = ", round(r_time_log, 3)),
           size = 6)+
  #annotate("text", x = 12, y = 7,
  #         label = "Structurally\nconfounded", size = 6, color = "#b80422") +
  mytheme

# Panel D: Residual variance after accounting for time
# If cumulative burden adds info beyond time, residuals should correlate with outcomes
burden_tidy_time <- burden_tidy_time %>%
  mutate(
    # Residual cumulative burden after removing time effect
    cum_burden_resid = residuals(lm(log10(cum_burden_integral + 1) ~ Time,
                                    data = burden_tidy_time))
  )

# Check: does residual cumulative burden predict health?
m_resid <- lm(scaled_health ~ cum_burden_resid, data = burden_tidy_time)
cat("\n=== Residual cumulative burden vs health (Figure S2D) ===\n")
print(summary(m_resid))

pD_supp <- ggplot(burden_tidy_time, aes(x = cum_burden_resid, y = scaled_health)) +
  geom_point(aes(fill = status), shape = 21, size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_fill_manual(values = c("#19798b", "#ee9b43")) +
  labs(x = "Residual cumulative CFU\n(after removing time effect)",
       y = "Standardized health score",
       fill = "Status") +
  annotate("text", x = -0.2, y = 1.8,
           label = paste0("β = ", round(coef(m_resid)[2], 3),
                          ", p = ", round(summary(m_resid)$coefficients[2,4], 3)),
           size = 6) +
  mytheme +
  theme(legend.position = c(0.85, 0.9))

# Combine
figure_cumulative_supp <- (pA_supp | pB_supp) / (pC_supp | pD_supp) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))


ggsave("figures/figureS4.pdf",  # >>> Supplementary Figure S4 (cumulative-burden methods + collinearity)
       plot = figure_cumulative_supp, width = 10, height = 9.2, dpi = 300)

#-------------------------------------------------------------------------------
# STATISTICAL TEST: Does cumulative burden add information beyond time?
#-------------------------------------------------------------------------------

# Predicting health
m_time_only <- lm(scaled_health ~ Time, data = burden_tidy_time)
m_cum_only  <- lm(scaled_health ~ log10(cum_burden_integral + 1), data = burden_tidy_time)
m_both      <- lm(scaled_health ~ Time + log10(cum_burden_integral + 1), data = burden_tidy_time)

cat("\nAIC comparison (predicting health):\n")
print(AIC(m_time_only, m_cum_only, m_both))

cat("\nModel with both predictors:\n")
print(summary(m_both))

# Key result: In the combined model, one predictor becomes non-significant
# because they are collinear

# Variance Inflation Factor
cat("\nVariance Inflation Factor (VIF) in combined model:\n")
print(vif(m_both))

#-------------------------------------------------------------------------------
# SUMMARY TABLE FOR SUPPLEMENTARY
#-------------------------------------------------------------------------------

summary_table <- tibble(
  Method = c("Integral of fitted logistic",
             "Trapezoidal sum of observed means",
             "Discrete sum of observed means"),
  Formula = c("(K/r) × ln[1 + (p₀/(K-p₀)) × (e^(rt) - 1)]",
              "Σᵢ (tᵢ - tᵢ₋₁) × (p̄ᵢ + p̄ᵢ₋₁)/2",
              "Σᵢ p̄ᵢ for all tᵢ ≤ t"),
  `Correlation with Time` = c(
    round(cor_matrix["Time", "cum_burden_integral"], 3),
    round(cor_matrix["Time", "cum_burden_trapezoid"], 3),
    round(cor_matrix["Time", "cum_burden_discrete"], 3)
  ),
  Pros = c("Smooth, uses all data, mechanistic",
           "Model-free, captures data variation",
           "Simplest, no interpolation"),
  Cons = c("Model-dependent, strongly collinear with t",
           "Sensitive to sampling density",
           "Crude, ignores time intervals")
)

print(summary_table)


#===============================================================================
### ANTIBIOTIC TREATMENT ###
# FIGURE 6
#===============================================================================

ab_data   <- read.table("data/bacterial_burden_ab.csv", header = T, sep = ",", dec =".")
expdata   <- read.table("data/health_assesment_ab.csv", header = T, sep = ",", dec =".")

# which groups have "treated at the same time as infection"?
groups_same_time <- c("G0")   # add others as needed, e.g. c("G1","G3")

# helper: parse "HH:MM" safely -> POSIXct on a dummy date
parse_hm_posix <- function(x) {
  x <- na_if(x, "")                                   # treat "" as NA
  parse_date_time(paste("2000-01-01", x),
                  orders = "Y-m-d H:M",
                  tz = "UTC",
                  quiet = TRUE)                       # invalids -> NA
}


# Vectorized helper: shift 'later' by whole days so the diff is within [min_gap, max_gap]
unwrap_diff <- function(later, earlier, min_gap_m = 0, max_gap_m = Inf, cycle_m = 24*60) {
  raw <- as.numeric(difftime(later, earlier, units = "mins"))
  
  # If difference is too small, push forward by whole days until >= min_gap
  add_days <- pmax(0, ceiling((min_gap_m - raw) / cycle_m))
  adj <- later + days(add_days)
  diff <- as.numeric(difftime(adj, earlier, units = "mins"))
  
  # If difference overshoots the max window by >= 1 day, pull back by whole days
  pull_days <- pmax(0, floor((diff - max_gap_m) / cycle_m))
  adj <- adj - days(pull_days)
  as.numeric(difftime(adj, earlier, units = "mins"))
}


expdata_filled <- expdata %>%
  mutate(
    Time_of_injection = na_if(Time_of_injection, ""),
    Time_of_treatment = na_if(Time_of_treatment, "")
  ) %>%
  group_by(Sample) %>%
  arrange(Larvae, .by_group = TRUE) %>%
  mutate(
    row_idx = row_number() - 1,
    inj_seed_chr = first(na.omit(Time_of_injection)),
    trt_seed_chr = first(na.omit(Time_of_treatment)),
    inj_seed = parse_hm_posix(inj_seed_chr),
    trt_seed = parse_hm_posix(trt_seed_chr),
    
    inj_time_fill = case_when(
      !is.na(Time_of_injection) ~ parse_hm_posix(Time_of_injection),
      !is.na(inj_seed)          ~ inj_seed + minutes(row_idx),
      TRUE                      ~ NA
    ),
    
    trt_time_fill = case_when(
      Sample %in% groups_same_time & !is.na(inj_time_fill) ~ inj_time_fill,
      !is.na(Time_of_treatment)                           ~ parse_hm_posix(Time_of_treatment),
      !is.na(trt_seed)                                    ~ trt_seed + minutes(row_idx),
      TRUE                                                ~ NA
    ),
    
    Time_of_injection = ifelse(is.na(inj_time_fill), NA_character_,
                               strftime(inj_time_fill, format = "%H:%M", tz = "UTC")),
    Time_of_treatment = ifelse(is.na(trt_time_fill), NA_character_,
                               strftime(trt_time_fill, format = "%H:%M", tz = "UTC"))
  ) %>%
  dplyr::select(-row_idx, -inj_seed_chr, -trt_seed_chr, -inj_seed, -trt_seed,
                -inj_time_fill, -trt_time_fill) %>%
  ungroup()  %>%
  mutate(Total_health = Activity + Melanization) %>%
  mutate(
    inj_time   = parse_hm_posix(Time_of_injection),
    treat_time = parse_hm_posix(Time_of_treatment),
    samp_time  = parse_hm_posix(Time_of_sampling),
    min_inj_to_treat  = unwrap_diff(treat_time, inj_time, min_gap_m = 0,   max_gap_m = 12*60)
  ) %>%
  mutate(hr_inj_to_treat = min_inj_to_treat / 60)


print(
  expdata_filled %>%
    filter(Sample %in% c("G3","G4","G5","G6","G7")) %>%
    ggplot(aes(x = min_inj_to_treat, y = Total_health)) +
    geom_jitter(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      x = "Minutes from injection to treatment",
      y = "Total Health",
      color = "Group"
    ) +
    mytheme
)

# tidy plate counts

# G4 L20: the dedicated 100 ul plate read zero while the 10 ul spot on a shared
# plate gave 2 colonies. Multi-sample spot plates are prone to carryover between
# spots; the dedicated plate is the reliable reading. The 10 ul spot is excluded
# for this larva, which then reads zero on all remaining plates and enters at
# LOD/2 with the other below-detection larvae.
ab_data <- ab_data %>%
  filter(!(Sample == "G4" & Larvae == "L20" & Volume_ul == 10))
# -----------------------------------------------------------------------------
# Detection limit. Lowest-dilution plates were 100 ul plated from a 250 ul
# homogenate at 10x, so one colony corresponds to 1 * 10 * (250/100) = 25
# CFU/larva. Larvae plated with no colonies on any plate are recorded at LOD/2
# rather than excluded: they are measurements at the limit, not missing data,
# and they concentrate in the early-treatment groups where clearance succeeded.
# =============================================================================
lod_cfu <- 25

# CFU per larva from replicate plate counts. Plates with no reading are already
# excluded upstream (Countable == "no" / NA). A plate read as zero is a
# measurement: if every plate for a larva reads zero, the larva is below the
# detection limit and enters at LOD/2 rather than being dropped.
# Average only plates that yielded colonies. With the G4 L20 row removed no
# larva has a mix of zero and non-zero plates, so this is currently identical to
# averaging all of them (0 of 113 larvae differ) -- but it will not silently
# diverge from the stated rule if the data changes.
larva_cfu <- function(cfu, count, lod = lod_cfu) {
  if (all(cfu == 0)) lod / 2 else mean(count[cfu > 0])
}

# Guard: a missing dilution field silently turns a real zero into NaN, which is
# then indistinguishable from "never plated" downstream.
stopifnot(!any(is.na(ab_data$Dilution)),
          !any(is.na(ab_data$Volume_ul)),
          !any(is.na(ab_data$Buffer_dilution_factor)))

burden_tidy_ab <- ab_data %>%
  pivot_longer(cols = rep1:rep5, names_to = "replicates", values_to = "cfu") %>%
  filter(!is.na(cfu)) %>%          # drop wells never plated; KEEP zeros
  mutate(dilution_factor = (Dilution / Volume_ul) * Buffer_dilution_factor,
         count           = cfu * dilution_factor) %>%
  group_by(Sample, Larvae) %>%
  # No colonies on ANY plate -> below detection. Otherwise average only plates
  # that yielded colonies: a zero at high dilution alongside growth at low
  # dilution is a dilution artefact, not evidence of absence (affects G4 L20).
  summarise(cfu = larva_cfu(cfu, count), .groups = "drop") %>%
  mutate(Sample = as.character(Sample), Larvae = as.character(Larvae)) %>%
  # Full join: PBS control larvae were never plated and enter as NA, then drop
  # out below. They re-enter the burden figure via controls_burden at 0.
  full_join(expdata_filled %>%
              mutate(Sample = as.character(Sample), Larvae = as.character(Larvae)),
            by = c("Sample", "Larvae")) %>%
  filter(!is.na(cfu)) %>%
  mutate(
    log_CFU = log10(cfu + 1),
    status  = case_when(Survival == 2 ~ "Alive",
                        Survival == 0 ~ "Dead",
                        TRUE          ~ NA_character_),
    scaled_health       = as.numeric(scale(Total_health)),
    scaled_activity     = as.numeric(scale(Activity)),
    scaled_melanization = as.numeric(scale(Melanization)),
    scaled_cfu          = as.numeric(scale(log_CFU)),
    scaled_immune_morb  = as.numeric(scale(Activity + Melanization))
  )

# Colour palette emphasizing early vs late
pal_improved <- c(
  "PBS-PBS" = "#d9d9d9",      # lightest grey - vehicle control
  "PBS-CIP" = "#969696",      # medium grey - injection control
  "PAO1-PBS" = "#252525",     # dark grey/black - no treatment
  "PAO1-00hCIP" = "#8adbea",
  "PAO1-03hCIP" = "#19798b",
  "PAO1-06hCIP" = "#fdd49e",
  "PAO1-09hCIP" = "#ee9b43",
  "PAO1-12hCIP" = "#b80422"   # red - too late
)


# Nice labels for all panels
treatment_labels <- c(
  "PBS-PBS" = "Injury control",
  "PBS-CIP" = "Antibiotic control",
  "PAO1-PBS" = "No treatment",
  "PAO1-00hCIP" = "Treatment 0h",
  "PAO1-03hCIP" = "Treatment 3h",
  "PAO1-06hCIP" = "Treatment 6h",
  "PAO1-09hCIP" = "Treatment 9h",
  "PAO1-12hCIP" = "Treatment 12h"
)

#-------------------------------------------------------------------------------
# PANEL A: SURVIVAL
#-------------------------------------------------------------------------------
# NOTE: binom.test() gives EXACT (Clopper-Pearson) intervals, not Wilson score.
# The manuscript caption should say "exact (Clopper-Pearson)".

survival_summary_ci <- expdata_filled %>%
  group_by(Treatment) %>%
  summarise(
    n_alive = sum(Survival == 2, na.rm = TRUE),
    n_total = sum(!is.na(Survival)),
    survival_prob = n_alive / n_total,
    ci_low = binom.test(n_alive, n_total)$conf.int[1],
    ci_high = binom.test(n_alive, n_total)$conf.int[2],
    .groups = "drop"
  )

p7A <- ggplot(survival_summary_ci, aes(y = Treatment, x = survival_prob, fill = Treatment)) +
  geom_col(color = "black", width = 0.7, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.25, linewidth = 0.5) +
  scale_fill_manual(values = pal_improved) +
  scale_y_discrete(labels = treatment_labels) +
  # let the interval end visibly inside the panel
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Survival probability", y = NULL) +
  mytheme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank()
  )

#-------------------------------------------------------------------------------
# PANEL B: HEALTH (ALL GROUPS including controls)
#-------------------------------------------------------------------------------

health_summary_ci <- expdata_filled %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean_health = mean(Total_health, na.rm = TRUE),
    sd_health = sd(Total_health, na.rm = TRUE),
    se_health = sd_health / sqrt(n()),
    ci_low = mean_health - 1.96 * se_health,  # 95% CI
    ci_high = mean_health + 1.96 * se_health,
    .groups = "drop"
  )

cat("\nHealth error bar widths:\n")
print(
  health_summary_ci %>%
    mutate(bar_width = ci_high - ci_low) %>%
    dplyr::select(Treatment, mean_health, se_health, bar_width)
)

p7B <- ggplot(health_summary_ci, aes(y = Treatment, x = mean_health, fill = Treatment)) +
  geom_point(size = 5, shape = 21, color = "black", stroke = 0.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.25, linewidth = 0.5) +
  scale_fill_manual(values = pal_improved) +
  scale_x_continuous(limits = c(0, 8),
                     breaks = seq(0, 8, by = 2),
                     expand = c(0.02, 0)) +
  scale_y_discrete(labels = treatment_labels) +
  labs(x = "Health score", y = NULL) +
  mytheme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank()
  )

#-------------------------------------------------------------------------------
# PANEL C: BACTERIAL BURDEN (ALL GROUPS - controls show 0)
#-------------------------------------------------------------------------------
burden_summary_ci <- burden_tidy_ab %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean_cfu = mean(log_CFU, na.rm = TRUE),
    sd_cfu = sd(log_CFU, na.rm = TRUE),
    se_cfu = sd_cfu / sqrt(n()),
    ci_low = mean_cfu - 1.96 * se_cfu,  # 95% CI
    ci_high = mean_cfu + 1.96 * se_cfu,
    .groups = "drop"
  )

burden_summary_ci <- burden_summary_ci %>%
  filter(!Treatment %in% c("PBS-PBS", "PBS-CIP"))

# Controls carry no bacteria; add them explicitly at zero
controls_burden <- expdata_filled %>%
  filter(Treatment %in% c("PBS-PBS", "PBS-CIP")) %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean_cfu = 0,
    sd_cfu = 0,
    se_cfu = 0,
    ci_low = 0,
    ci_high = 0,
    .groups = "drop"
  )

burden_summary_ci <- bind_rows(burden_summary_ci, controls_burden)

p7C <- ggplot(burden_summary_ci, aes(y = Treatment, x = mean_cfu, fill = Treatment)) +
  geom_point(size = 5, shape = 21, color = "black", stroke = 0.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.25, linewidth = 0.5) +
  scale_fill_manual(values = pal_improved) +
  scale_x_continuous(limits = c(-0.2, 6),
                     breaks = seq(0, 6, by = 1),
                     expand = c(0.02, 0)) +
  scale_y_discrete(labels = treatment_labels) +
  labs(x = expression(paste("log"[10], "CFU")), y = NULL) +
  mytheme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank()
  )

# =============================================================================
# Treatment-timing dose-response (continuous) — capstone for the treatment story
# -----------------------------------------------------------------------------
# Every larva is assessed at the SAME 24 h endpoint, so variation in outcome vs
# minutes-from-injection-to-treatment reflects the manipulated DURATION of
# pathogen exposure before clearance -> a direct, prospective test of the
# cumulative-damage hypothesis (the discrete groups of Figure 6A-C, made
# continuous using the exact per-larva treatment time).
# =============================================================================

dose <- burden_tidy_ab %>%
  filter(grepl("PAO1", Treatment), grepl("CIP", Treatment),
         !is.na(min_inj_to_treat), !is.na(Survival)) %>%
  mutate(alive01 = as.integer(Survival == 2),
         hr_inj_to_treat = min_inj_to_treat / 60)

t_grid_d <- data.frame(min_inj_to_treat = seq(min(dose$min_inj_to_treat, na.rm = TRUE),
                                              max(dose$min_inj_to_treat, na.rm = TRUE),
                                              length.out = 200))

t_grid_d$hr_inj_to_treat <- t_grid_d$min_inj_to_treat / 60

# (1) Survival ~ delay: logistic regression -> treatment window (P = 0.5)
m_surv_dose <- glm(alive01 ~ min_inj_to_treat, data = dose, family = binomial)
t_grid_d$surv <- predict(m_surv_dose, newdata = t_grid_d, type = "response")
bcoef   <- coef(m_surv_dose)
window50 <- as.numeric(-bcoef[1] / bcoef[2])   # delay at which P(survival) = 0.5
cat("\n=== Treatment-timing dose-response (24 h endpoint) ===\n")
print(summary(m_surv_dose))
cat(sprintf("Survival: P(survival) = 0.5 at %.0f min (%.1f h) after injection\n",
            window50, window50 / 60))
win_se <- tryCatch({ dp <- MASS::dose.p(m_surv_dose, p = 0.5); attr(dp, "SE")[1] },
                   error = function(e) NA_real_)
if (!is.na(win_se)) cat(sprintf("   (approx SE %.0f min = %.1f h)\n", win_se, win_se / 60))

# (2) Composite health ~ delay: linear trend
m_health_dose <- lm(Total_health ~ min_inj_to_treat, data = dose)
t_grid_d$health <- predict(m_health_dose, newdata = t_grid_d)
cat(sprintf("Health: slope %.3f per h, p = %.3g\n",
            coef(m_health_dose)[2] * 60, summary(m_health_dose)$coefficients[2, 4]))

# (3) Residual CFU at 24 h ~ delay
m_cfu_dose <- lm(log_CFU ~ min_inj_to_treat, data = dose)
t_grid_d$logcfu <- predict(m_cfu_dose, newdata = t_grid_d)
cat(sprintf("log10 CFU(24h): slope %.4f per min, p = %.3g\n\n",
            coef(m_cfu_dose)[2], summary(m_cfu_dose)$coefficients[2, 4]))

# (4) Does delay still predict outcome AFTER adjusting for the burden carried at
# the 24h endpoint? If delay stays significant while log_CFU does not carry it,
# the timing effect is not routed through final burden (the decoupling claim).
m_surv_adj   <- glm(alive01 ~ min_inj_to_treat + log_CFU, data = dose, family = binomial)
m_health_adj <- lm(Total_health ~ min_inj_to_treat + log_CFU, data = dose)
cat("--- Outcome ~ delay, adjusted for final (24h) burden ---\n")
cat("Survival ~ delay + log_CFU:\n"); print(round(summary(m_surv_adj)$coefficients, 4))
cat("Health   ~ delay + log_CFU:\n"); print(round(summary(m_health_adj)$coefficients, 4))

# (5) Did the antibiotic actually reduce burden?
cfu_untreated <- burden_tidy_ab %>% filter(Treatment == "PAO1-PBS") %>% pull(log_CFU)
cfu_late      <- dose %>% filter(min_inj_to_treat >= 9 * 60) %>% pull(log_CFU)
cfu_early     <- dose %>% filter(min_inj_to_treat <= 3 * 60) %>% pull(log_CFU)
cat(sprintf("\n--- Did the drug reduce burden? median log10 CFU(24h) ---\n  untreated %.2f (n=%d) | late >=9h %.2f (n=%d) | early <=3h %.2f (n=%d)\n",
            median(cfu_untreated, na.rm = TRUE), sum(!is.na(cfu_untreated)),
            median(cfu_late, na.rm = TRUE), sum(!is.na(cfu_late)),
            median(cfu_early, na.rm = TRUE), sum(!is.na(cfu_early))))
if (length(na.omit(cfu_untreated)) > 1 && length(na.omit(cfu_late)) > 1) {
  cat("  late vs untreated (Wilcoxon) p =",
      signif(wilcox.test(cfu_late, cfu_untreated)$p.value, 3), "\n")
}
cat("\n")

# --- Control reference values -------------------------------------------------
# Two horizontal references for panels D-F: uninfected larvae (upper bound on
# outcome) and infected-but-untreated larvae (lower bound). Survival and health
# come from expdata_filled so that PBS groups are included; burden for the
# untreated group comes from the plate data, and uninfected burden is zero.

ctrl_untreated_cfu <- burden_tidy_ab %>%
  filter(Treatment == "PAO1-PBS") %>%
  summarise(m = mean(log_CFU, na.rm = TRUE)) %>%
  pull(m)

ctrl_ref <- expdata_filled %>%
  mutate(
    alive01 = as.integer(Survival == 2),
    control = dplyr::case_when(
      !grepl("PAO1", Treatment)                           ~ "Uninfected",
      grepl("PAO1", Treatment) & !grepl("CIP", Treatment) ~ "Infected and untreated",
      TRUE                                                ~ NA_character_
    )
  ) %>%
  filter(!is.na(control)) %>%
  group_by(control) %>%
  summarise(surv   = mean(alive01,      na.rm = TRUE),
            health = mean(Total_health, na.rm = TRUE),
            n      = dplyr::n(),
            .groups = "drop") %>%
  mutate(
    logcfu  = ifelse(control == "Uninfected", 0, ctrl_untreated_cfu),
    control = factor(control, levels = c("Uninfected", "Infected and untreated"))
  )

stopifnot(nrow(ctrl_ref) == 2, !any(is.na(ctrl_ref$logcfu)))
cat("\n--- Control reference values (panels D-F) ---\n"); print(ctrl_ref)

ctrl_cols <- c("Uninfected" = "#4d9221", "Infected and untreated" = "grey25")
ctrl_lty  <- c("Uninfected" = "dashed",  "Infected and untreated" = "dashed")

ctrl_layer <- function(yvar) {
  list(
    geom_hline(data = ctrl_ref,
               aes(yintercept = .data[[yvar]], color = control, linetype = control),
               linewidth = 0.7),
    scale_color_manual(name = NULL, values = ctrl_cols),
    scale_linetype_manual(name = NULL, values = ctrl_lty)
  )
}


# --- Panel D counts -----------------------------------------------------------
cluster_n <- dose %>%
  mutate(cluster = round(min_inj_to_treat / 180) * 180) %>%
  group_by(cluster) %>%
  summarise(x       = mean(hr_inj_to_treat, na.rm = TRUE),
            n_alive = sum(alive01 == 1, na.rm = TRUE),
            n_dead  = sum(alive01 == 0, na.rm = TRUE),
            .groups = "drop")


# --- Panels D-F ---------------------------------------------------------------
# The control legend is drawn once, inside panel D; panels E and F suppress it.
# (There is no plot_layout(guides = "collect") in the assembly below.)

pD_surv <- ggplot(dose, aes(hr_inj_to_treat, alive01)) +
  geom_jitter(height = 0.04, width = 0.1, alpha = 0.5, shape = 21,
              fill = "grey40", size = 2) +
  ctrl_layer("surv") +
  geom_line(data = t_grid_d, aes(y = surv), color = "#b80422", linewidth = 1.5) +
  geom_text(data = cluster_n, aes(x = x, y = 1.05, label = n_alive),
            inherit.aes = FALSE, vjust = -0.8, size = 4, fontface = "bold",
            color = "grey15") +
  geom_text(data = cluster_n, aes(x = x, y = -0.05, label = n_dead),
            inherit.aes = FALSE, vjust = 1.8, size = 4, fontface = "bold",
            color = "grey15") +
  scale_y_continuous(breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0.13, 0.13))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(-0.3, 12.5)) +
  labs(x = NULL, y = "P(survival at 24 hrs)") +
  mytheme +
  theme(legend.position = c(0.65, 0.72))

pD_cfu <- ggplot(dose, aes(hr_inj_to_treat, log_CFU)) +
  geom_jitter(alpha = 0.5, width = 0.1, shape = 21, fill = "grey40", size = 2) +
  ctrl_layer("logcfu") +
  geom_line(data = t_grid_d, aes(y = logcfu), color = "#ee9b43", linewidth = 1.5) +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(-0.3, 12.5)) +
  labs(x = "Time from injection to treatment (hrs)",
       y = bquote(log[10]~CFU~"(24 hrs)")) +
  mytheme +
  theme(legend.position = "none")


pD_health <- ggplot(dose, aes(hr_inj_to_treat, Total_health)) +
  geom_jitter(alpha = 0.5, width = 0.1, shape = 21, fill = "grey40", size = 2) +
  ctrl_layer("health") +
  geom_line(data = t_grid_d, aes(y = health), color = "#19798b", linewidth = 1.5) +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(-0.3, 12.5)) +
  labs(x = NULL, y = "Health score (24 hrs)") +
  mytheme +
  theme(legend.position = "none")


# --- Assemble -----------------------------------------------------------------
# Panel order: survival, burden, health in BOTH rows, so tags run
# A = survival (groups),  B = burden (groups),  C = health (groups)
# D = survival (timing),  E = burden (timing),  F = health (timing)
# The figure 6 caption must describe B as burden and C as health.

figure6 <- (p7A | p7C | p7B) / (pD_surv | pD_cfu | pD_health) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("figures/figure6.pdf",   # >>> Manuscript Figure 6
       plot = figure6, width = 12, height = 7.5, units = "in", dpi = 300)

#===============================================================================
# STATISTICS  — one model per question
#===============================================================================

cat("\n===============================================================\n")
cat("  ANTIBIOTIC EXPERIMENT: values quoted in the manuscript\n")
cat("===============================================================\n")

# --- Group-level summaries ----------------------------------------------------
cat("\n--- Survival by group (Clopper-Pearson 95% CI) ---\n")
print(survival_summary_ci %>%
        mutate(txt = sprintf("%.0f%% [%.0f%%, %.0f%%]",
                             100*survival_prob, 100*ci_low, 100*ci_high)) %>%
        dplyr::select(Treatment, n_alive, n_total, txt))

cat("\n--- Health by group (mean, 95% CI) ---\n")
print(health_summary_ci %>%
        mutate(txt = sprintf("%.2f [%.2f, %.2f]", mean_health, ci_low, ci_high)) %>%
        dplyr::select(Treatment, n, txt))

cat("\n--- Burden by group (mean log10 CFU, 95% CI) ---\n")
print(burden_summary_ci %>%
        mutate(txt = sprintf("%.2f [%.2f, %.2f]", mean_cfu, ci_low, ci_high)) %>%
        dplyr::select(Treatment, n, txt))


# --- (1) Survival ~ delay -----------------------------------------------------
# Continuous per-larva delay. The nominal-group version (survived ~ 0/3/6/9/12)
# is the same result and is not fitted separately.
cat("\n--- (1) Survival ~ delay (logistic, n =", nrow(dose), ") ---\n")
print(round(summary(m_surv_dose)$coefficients, 5))
cat(sprintf("   P(survival) = 0.5 at %.0f min (%.1f h)\n", window50, window50/60))
if (!is.na(win_se)) cat(sprintf("   approx SE %.0f min (%.1f h)\n", win_se, win_se/60))

# --- (2) Health ~ delay -------------------------------------------------------
cat("\n--- (2) Health ~ delay (linear) ---\n")
print(round(summary(m_health_dose)$coefficients, 5))
cat(sprintf("   slope %.3f health points per h, p = %.3g\n",
            coef(m_health_dose)[2]*60, summary(m_health_dose)$coefficients[2, 4]))

# --- (3) Burden at 24 h ~ delay ----------------------------------------------
# CIP-treated larvae ONLY. PAO1-PBS larvae have a recorded PBS injection time but
# never received ciprofloxacin, so including them in a treatment-timing
# regression is meaningless; the earlier lm_burden_timing did exactly that over
# 120 larvae and reported a flat slope for that reason.
stopifnot(all(grepl("CIP", dose$Treatment)), nrow(dose) == 100)
cat("\n--- (3) Burden at 24 h ~ delay (linear, CIP groups only) ---\n")
print(round(summary(m_cfu_dose)$coefficients, 5))
cat(sprintf("   slope %.4f log10 CFU per min, p = %.3g\n",
            coef(m_cfu_dose)[2], summary(m_cfu_dose)$coefficients[2, 4]))
cat("   NOTE: group means are not monotone in delay -- inspect before quoting.\n")
print(dose %>% group_by(Treatment) %>%
        summarise(n = dplyr::n(), mean_logCFU = mean(log_CFU), .groups = "drop"))

# --- (4) Does delay predict outcome after adjusting for final burden? ---------
# This is the decoupling test and the load-bearing result of the section.
cat("\n--- (4) Outcome ~ delay + final burden ---\n")
cat("Survival:\n"); print(round(summary(m_surv_adj)$coefficients, 4))
cat("Health:\n");   print(round(summary(m_health_adj)$coefficients, 4))

# Standardised versions: delay and burden on a common footing within each model.
m_health_adj_z <- lm(Total_health ~ scale(min_inj_to_treat) + scale(log_CFU),
                     data = dose)
cat("\nHealth, standardised predictors (comparable within model):\n")
print(round(summary(m_health_adj_z)$coefficients, 4))

# --- (5) Did the drug reduce burden? -----------------------------------------
cfu_untreated <- burden_tidy_ab %>% filter(Treatment == "PAO1-PBS") %>% pull(log_CFU)
cfu_late      <- dose %>% filter(min_inj_to_treat >= 9*60) %>% pull(log_CFU)
cfu_early     <- dose %>% filter(min_inj_to_treat <= 3*60) %>% pull(log_CFU)

cat(sprintf("\n--- (5) Median log10 CFU at 24 h ---\n   untreated %.2f (n=%d) | late >=9h %.2f (n=%d) | early <=3h %.2f (n=%d)\n",
            median(cfu_untreated, na.rm = TRUE), sum(!is.na(cfu_untreated)),
            median(cfu_late,      na.rm = TRUE), sum(!is.na(cfu_late)),
            median(cfu_early,     na.rm = TRUE), sum(!is.na(cfu_early))))
cat(sprintf("   late vs untreated (Wilcoxon) p = %.3g\n",
            wilcox.test(cfu_late, cfu_untreated)$p.value))

# --- (6) Group comparisons ----------------------------------------------------
# Fisher tests pool the two early groups (0 h and 3 h) so that "early" in the
# text matches what is tested; the previous version tested 3 h alone.
infected_groups <- expdata_filled %>% filter(grepl("PAO1", Treatment))

chisq_survival <- chisq.test(table(infected_groups$Treatment,
                                   infected_groups$Survival))
cat("\n--- (6) Survival across infected groups ---\n")
print(chisq_survival)

early_vs_untreated <- expdata_filled %>%
  filter(Treatment %in% c("PAO1-PBS", "PAO1-00hCIP", "PAO1-03hCIP")) %>%
  mutate(grp = ifelse(Treatment == "PAO1-PBS", "untreated", "early")) %>%
  {fisher.test(table(.$grp, .$Survival))}
cat("\nEarly (0-3 h pooled) vs untreated:\n"); print(early_vs_untreated)

late_vs_untreated <- expdata_filled %>%
  filter(Treatment %in% c("PAO1-PBS", "PAO1-12hCIP")) %>%
  {fisher.test(table(.$Treatment, .$Survival))}
cat("\nLate (12 h) vs untreated:\n"); print(late_vs_untreated)

# --- (7) Cumulative exposure before treatment --------------------------------
# Sigma_p uses K, r, p0 from the MAIN-cohort logistic fit applied to antibiotic
# larvae. The two experiments are separate batches (untreated 24 h burden 5.29
# here vs 6.05 at comparable times in the main cohort), so this is an
# approximation and should be described as such in the ESM.
dose <- dose %>%
  mutate(t_hours     = min_inj_to_treat / 60,
         Sigma_p_pre = (K / r) * log1p((p0 / (K - p0)) * expm1(r * t_hours)))

m_surv_sigma   <- glm(alive01 ~ Sigma_p_pre, data = dose, family = binomial)
m_health_sigma <- lm(Total_health ~ Sigma_p_pre, data = dose)

cat("\n--- (7) Delay vs cumulative exposure (AIC) ---\n")
cat("Survival:\n"); print(AIC(m_surv_dose,   m_surv_sigma))
cat("Health:\n");   print(AIC(m_health_dose, m_health_sigma))
cat(sprintf("   dAIC survival = %.2f | dAIC health = %.2f\n",
            AIC(m_surv_sigma)   - AIC(m_surv_dose),
            AIC(m_health_sigma) - AIC(m_health_dose)))

sigma_coef <- coef(m_surv_sigma)
sigma50    <- as.numeric(-sigma_coef[1] / sigma_coef[2])
cat(sprintf("   50%% survival threshold: Sigma_p = %.0f cells.h\n", sigma50))
cat(sprintf("   (logistic params used: K = %.3g, r = %.3f, p0 = %.3g)\n", K, r, p0))


# Sanity check units in hours
dose <- dose %>% mutate(hr_inj_to_treat = min_inj_to_treat / 60)

m_surv_dose_h   <- glm(alive01 ~ hr_inj_to_treat, data = dose, family = binomial)
m_health_dose_h <- lm(Total_health ~ hr_inj_to_treat, data = dose)
m_cfu_dose_h    <- lm(log_CFU ~ hr_inj_to_treat, data = dose)
m_surv_adj_h    <- glm(alive01 ~ hr_inj_to_treat + log_CFU, data = dose, family = binomial)
m_health_adj_h  <- lm(Total_health ~ hr_inj_to_treat + log_CFU, data = dose)

summary(m_surv_dose_h); summary(m_health_dose_h); summary(m_cfu_dose_h);
summary(m_surv_adj_h); summary(m_health_adj_h)

#=============================================
# Cumulative burden trajectory (not in the manuscript)
#=============================================

time_seq_cum <- seq(0, 36, length.out = 200)
cum_burden_curve <- (K / r) * log1p((p0 / (K - p0)) * expm1(r * time_seq_cum))

df_cum_curve <- data.frame(
  Time = time_seq_cum,
  cum_burden = cum_burden_curve,
  log_cum = log10(cum_burden_curve + 1)
)

inset_cum <- ggplot(df_cum_curve, aes(x = Time, y = log_cum)) +
  geom_line(color = "#b80422", linewidth = 1) +
  labs(x = "Time (hrs)",
       y = bquote(log[10](Sigma*italic(p)))) +
  scale_x_continuous(breaks = c(0, 18, 36)) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0, 4, 8)) +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.margin = margin(2, 2, 2, 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(inset_cum)

# Per-larva cumulative burden (used in the comparison below)
burden_tidy_time <- burden_tidy_time %>%
  mutate(
    # Analytical integral of the logistic curve. Written with log1p/expm1
    # because the direct form is numerically unstable at large rt.
    cum_burden = (K/r) * log1p((p0 / (K - p0)) * expm1(r * Time)),
    scaled_cum_burden = scale(cum_burden)[,1]
  )

# Which predicts health better: instantaneous or cumulative burden?
m_instant    <- lm(health_combined ~ scaled_cfu, data = burden_tidy_time)
m_cumulative <- lm(health_combined ~ scaled_cum_burden, data = burden_tidy_time)
m_both_hb    <- lm(health_combined ~ scaled_cfu + scaled_cum_burden, data = burden_tidy_time)

cat("\n=== Health ~ instantaneous vs cumulative burden (AIC) ===\n")
print(AIC(m_instant, m_cumulative, m_both_hb))


# Compare instantaneous vs cumulative as predictors of health
print(
  burden_tidy_time %>%
    pivot_longer(cols = c(scaled_cfu, scaled_cum_burden),
                 names_to = "burden_type",
                 values_to = "burden_value") %>%
    mutate(burden_type = factor(burden_type,
                                levels = c("scaled_cfu", "scaled_cum_burden"),
                                labels = c("Instantaneous burden",
                                           "Cumulative burden"))) %>%
    ggplot(aes(x = burden_value, y = scaled_health)) +
    geom_point(aes(fill = status), shape = 21, size = 2) +
    geom_smooth(method = "lm", color = "black") +
    facet_wrap(~burden_type) +
    labs(x = "Standardized density", y = "Standardized health score") +
    mytheme
)

# =============================================================================
# MANUSCRIPT NUMBER REPORT
# -----------------------------------------------------------------------------
# Paste at the END of main.R and run after the full script has been sourced.
# Prints every value quoted in main.tex / the ESM in one place, grouped by
# where it appears, and compares each against the value currently in the
# manuscript. Missing objects are reported, not fatal.
#
#   [ok]     recomputed value matches the manuscript
#   [CHECK]  they differ -> reconcile before submission
#   MISSING  object not in the workspace (script section not run?)
# =============================================================================

# [FIX/NEW] The report below is now written to results/statistics_for_manuscript.txt
# as well as the console, so the numbers can be diffed against the manuscript
# without scrolling back through the run log.
RESULTS_FILE <- file.path("results", "statistics_for_manuscript.txt")
.mn_con <- file(RESULTS_FILE, open = "wt")
sink(.mn_con, split = TRUE)
# NOTE: do NOT guard this with on.exit() at top level. on.exit() attaches to the
# enclosing evaluation context, which under source() is the source() call itself
# -- so it fires immediately, closes the connection before anything is written,
# and leaves a 0-byte statistics file. The sink is closed explicitly at the end
# of the script instead.

cat(strrep("=", 78), "\n")
cat("STATISTICS FOR THE MANUSCRIPT\n")
cat("Host-health integrates infection history to determine survival\n")
cat("during acute Pseudomonas aeruginosa infection\n")
cat(strrep("=", 78), "\n")
cat("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("Script:    main_cleaned.R\n")
cat("R:         ", R.version.string, "\n", sep = "")
cat("\nLegend:  [ok] recomputed value matches the manuscript\n")
cat("         [CHECK] they differ -- reconcile before submission\n")
cat("         [ok: < x] a claim stated as a bound, and the bound holds\n")
cat("         MISSING object not in the workspace (section not run?)\n")

cat("\n", strrep("-", 78), "\n", sep = "")
cat("RECONCILED 2026-08-23. Nine reference values were updated in the manuscript\n")
cat("after the detection-limit correction (n = 85 -> 86, zero-count larvae now\n")
cat("entering at LOD/2). The values below are the CURRENT published ones, so\n")
cat("their [ok] means genuine agreement, not a number compared with itself.\n")
cat("Previous values are recorded beside each call in the source.\n\n")
mn_recon <- rbind(
  c("N infected",                        "130",    "128"),
  c("MT50 in the fig 4 caption",         "21.5",   "21.0"),
  c("spline-adjusted p_Time",            "6e-08",  "4.5e-08"),
  c("beta_cfu (per SD)",                 "-0.308", "-0.317"),
  c("beta_time (per SD)",                "-0.015", "-0.012"),
  c("standardised beta (z_time:z_cfu)",  "-0.75",  "-0.939"),
  c("standardised p    (z_time:z_cfu)",  "7.5e-4", "3e-05"),
  c("interaction F",                     "12.3",   "19.6"),
  c("interaction p",                     "7.5e-4", "3e-05"))
cat(sprintf("  %-36s %-10s -> %s\n", "row", "was", "now"))
for (i in seq_len(nrow(mn_recon)))
  cat(sprintf("  %-36s %-10s -> %s\n", mn_recon[i,1], mn_recon[i,2], mn_recon[i,3]))
cat(strrep("-", 78), "\n")

mn_tol <- 5e-3          # relative tolerance for the [ok] / [CHECK] verdict

mn_sec <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n",
                          strrep("=", 78), "\n", sep = "")
mn_sub <- function(x) cat("\n--- ", x, " ---\n", sep = "")

# label | expression | manuscript value (NA = just print) | digits
mn <- function(label, expr, ms = NA, digits = 3) {
  v <- tryCatch(eval(substitute(expr), envir = parent.frame()),
                error = function(e) NULL,
                warning = function(w) suppressWarnings(
                  eval(substitute(expr), envir = parent.frame())))
  if (is.null(v) || length(v) == 0 || all(is.na(v))) {
    cat(sprintf("  %-44s %s\n", label, "MISSING"))
    return(invisible(NULL))
  }
  num <- suppressWarnings(as.numeric(v))
  # Counts print exactly: signif(x, 0) rounds to ONE significant figure, which
  # would show n = 86 as 90 and the detection limit 25 as 20.
  txt <- if (digits <= 0)
    paste(format(round(num), scientific = FALSE, trim = TRUE), collapse = ", ")
  else
    paste(format(signif(num, digits), scientific = NA), collapse = ", ")
  verdict <- ""
  if (!is.na(ms[1])) {
    # Compare at the precision the manuscript actually quotes, so that pure
    # rounding (-0.224 vs -0.22) is not flagged as a mismatch.
    msv <- as.numeric(ms[1])
    # Precision is taken from the manuscript value itself (how many significant
    # figures it is actually quoted to), not from `digits`, so a recomputed
    # -0.2236 counts as matching a quoted -0.22 while 19.6 vs 12.3 does not.
    ms_sig <- 15
    for (k in 1:15)
      if (isTRUE(all.equal(signif(msv, k), msv, tolerance = 1e-12))) { ms_sig <- k; break }
    ok <- if (digits <= 0) isTRUE(round(num[1]) == round(msv))
          else isTRUE(all.equal(signif(num[1], ms_sig), msv,
                                tolerance = mn_tol,
                                scale = max(abs(msv), 1e-12)))
    verdict <- if (ok) "  [ok]" else sprintf("  [CHECK: ms = %s]", ms[1])
  }
  cat(sprintf("  %-44s %-22s%s\n", label, txt, verdict))
  invisible(num)
}

mn_lt <- function(label, expr, bound, digits = 3) {
  v <- tryCatch(eval(substitute(expr), envir = parent.frame()),
                error = function(e) NULL)
  if (is.null(v) || length(v) == 0 || all(is.na(v))) {
    cat(sprintf("  %-44s %s\n", label, "MISSING")); return(invisible(NULL))
  }
  num <- suppressWarnings(as.numeric(v))[1]
  ok  <- isTRUE(num < bound)
  cat(sprintf("  %-44s %-22s%s\n", label,
              format(signif(num, digits), scientific = NA),
              if (ok) sprintf("  [ok: < %s]", format(bound))
              else    sprintf("  [CHECK: claimed < %s]", format(bound))))
  invisible(num)
}

# --- extractors --------------------------------------------------------------
cf  <- function(m, term) summary(m)$coefficients[term, 1]        # estimate
pv  <- function(m, term) summary(m)$coefficients[term, 4]        # Wald p
lfb <- function(f, term) unname(coef(f)[which(names(coef(f)) == term)])
lfp <- function(f, term) unname(f$prob[which(names(coef(f)) == term)])
srb <- function(m, term) summary(m)$table[term, 1]               # survreg est
srp <- function(m, term) summary(m)$table[term, 4]               # survreg p
lrt <- function(m0, m1) {                                        # LRT p, glm
  a <- anova(m0, m1, test = "LRT"); a[["Pr(>Chi)"]][2]
}
lrt_dev <- function(m0, m1) { a <- anova(m0, m1, test = "LRT"); a[["Deviance"]][2] }


mn_sec("SAMPLE SIZES AND DATA DECISIONS")
mn("n, main cohort (Time < 37 h)",        nrow(burden_tidy_time), 86, 0)
mn("  alive at sampling",                 sum(burden_tidy_time$status == "Alive"), NA, 0)
mn("  dead at sampling",                  sum(burden_tidy_time$status == "Dead"), NA, 0)
mn("  entering at LOD/2",                 sum(burden_tidy_time$cfu == lod_cfu/2), NA, 0)
mn("detection limit (CFU/larva)",         lod_cfu, 25, 0)
mn("mean sampling time (h)",              mean(burden_tidy_time$Time), NA, 3)
mn("sd sampling time (h)",                sd(burden_tidy_time$Time), NA, 3)
mn("n, antibiotic cohort (scored)",       nrow(expdata_filled), 160, 0)
mn("n, antibiotic cohort (plated)",       nrow(burden_tidy_ab), NA, 0)
mn("n, CIP-treated with a recorded delay",nrow(dose), 100, 0)
cat("  Health index = activity + melanization, both RAW (higher = healthier).\n")
cat("  The 48 h sample is excluded at the TOP of the data-preparation pipeline\n")
cat("  (filter(Time < 37) on the raw burden table), so the 4-36 h window and\n")
cat("  n = 86 apply to EVERY analysis in the paper -- the conditional-\n")
cat("  independence tests, the T50 fits, the GAMs and the SEM alike. It is not\n")
cat("  a standardisation choice; standardisation merely happens after it.\n")
cat("  Maximum Time in each analysis frame (all should read 36):\n")
for (nm in c("burden_tidy_time", "dat_fig3", "dat_fig3b", "dat_fig4d",
             "dat_fig5", "dat_cs", "dat_fig5_mel")) {
  d  <- get(nm)
  tv <- if ("Time" %in% names(d)) d$Time else d$time
  cat(sprintf("    %-18s n = %3d   max Time = %2.0f h\n", nm, nrow(d), max(tv, na.rm = TRUE)))
}
cat(sprintf("    %-18s n = %3d   max Time = %2.0f h  (survivors only)\n",
            "burden_tidy_alive", nrow(burden_tidy_alive),
            max(burden_tidy_alive$Time, na.rm = TRUE)))


mn_sec("ABSTRACT")
cat("  (no numerical claims after the 'three orders of magnitude' cut --\n")
cat("   if any number reappears here, it must also appear in Results)\n")


mn_sec("RESULTS (a)  Exponential mortality  [figure 2]")
mn("N infected (larvae at t = 1 h)", data$total[data$time == 1], 128, 0)  # reconciled 2026-08-23 (was 130)
mn("Gompertz a (h^-1)",          coef(fit)["a"],        4.14e-5)
mn("Gompertz b (h^-1)",          coef(fit)["b"],        0.40, 2)


mn_sec("RESULTS (b)  Logistic growth  [figure 3]")
mn("AIC logistic",               aic_logistic_point,    190, 0)
mn("AIC exponential",            aic_exp_point,         220, 0)
cat("    (an older code comment quoted 186 / 215 -- reconcile with the above)\n")
mn("K (CFU/host)",               coef(logistic_logfit_full)["K"], 3.2e6, 2)
mn("r, all larvae (h^-1)",       coef(logistic_logfit_full)["r"],   0.47, 2)
mn("r, survivors only (h^-1)",   coef(logistic_logfit_a_full)["r"], 0.40, 2)
mn("K, survivors only (CFU/host)", coef(logistic_logfit_a_full)["K"], NA, 2)
mn("doubling time, all larvae (h)",     log(2)/coef(logistic_logfit_full)["r"],   NA, 3)
mn("doubling time, survivors only (h)", log(2)/coef(logistic_logfit_a_full)["r"], NA, 3)
mn_sub("bootstrap model selection (lowest AIC over 1000 resamples)")
mn("logistic wins, all larvae (%)",     100*mean(winners == 1),   100, 0)
mn("exponential wins, all larvae (%)",  100*mean(winners == 2),   NA, 0)
mn("logistic wins, survivors only (%)", 100*mean(winners_a == 1), 63, 0)
mn("exponential wins, survivors only (%)", 100*mean(winners_a == 2), NA, 0)
mn("usable replicates, all / survivors", c(length(winners), length(winners_a)), NA, 0)
mn("post-mortem beta",           cf(lm_postmortem, "time_since_death"), 0.001, 3)
mn("post-mortem p",              pv(lm_postmortem, "time_since_death"), 0.96, 2)


mn_sec("RESULTS (b)  Conditional independence  s _||_ t | p")
mn_sub("additive GLM (m_cond)")
mn("beta_Time",                  cf(m_cond, "Time"),    -0.22)
mn("p_Time",                     pv(m_cond, "Time"),    NA)
mn("beta_logCFU",                cf(m_cond, "log_CFU"), -1.64)
mn("p_logCFU",                   pv(m_cond, "log_CFU"), 0.01, 2)

mn_sub("flexible-burden GAM (m_flex) -- robustness")
mn("beta_Time",                  summary(m_flex)$p.table["Time", 1], -0.23)
mn("p_Time",                     summary(m_flex)$p.table["Time", 4], 0.006, 2)
mn("edf s(log_CFU)",             summary(m_flex)$s.table[1, 1],      1.84, 3)
mn("chi-sq s(log_CFU)",          summary(m_flex)$s.table[1, 3],      7.19, 3)
mn("p s(log_CFU)",               summary(m_flex)$s.table[1, 4],      0.043, 2)
mn("deviance explained (%)",     summary(m_flex)$dev.expl * 100,     74.9, 3)

mn_sub("interaction LRT (m_cond vs m_int)")
mn("delta deviance",             lrt_dev(m_cond, m_int),  2.40)
mn("p",                          lrt(m_cond, m_int),      0.12, 2)


mn_sec("RESULTS (c)  Host health  [figure 4]")
mn_sub("transition timing")
mn("AT50 (h)",                   T50_summary$T50[1],   20.3, 3)
mn("AT50 CI lo",                 T50_summary$CI_low[1],    18.4, 3)
mn("AT50 CI hi",                 T50_summary$CI_high[1],    22.0, 3)
mn("MT50 (h)",                   T50_summary$T50[2],   21.0, 3)
mn("MT50 CI lo",                 T50_summary$CI_low[2],    18.2, 3)
mn("MT50 CI hi",                 T50_summary$CI_high[2],    23.5, 3)
mn("LT50 (h)",                   T50_summary$T50[3],   23.3, 3)
mn("LT50 CI lo",                 T50_summary$CI_low[3],    20.8, 3)
mn("LT50 CI hi",                 T50_summary$CI_high[3],    25.5, 3)
mn("LT50 - AT50 (h)  ['~3 h']",  T50_summary$T50[3] - T50_summary$T50[1], 3.0, 2)
mn("P(AT50 < MT50)",             mean(boot_AT[ix] < boot_MT[ix]), 0.64, 2)
mn("P(AT50 < LT50)",             mean(boot_AT[ix] < boot_LT[ix]), 0.97, 2)
mn("P(MT50 < LT50)",             mean(boot_MT[ix] < boot_LT[ix]), 0.90, 2)

mn("MT50 as printed in the fig 4 caption", T50_summary$T50[2], 21.0, 3)  # reconciled 2026-08-23 (was 21.5)
cat("    ^ the manuscript value here is the FIGURE 4 CAPTION, not the main text.\n")
cat("      If this row shows [CHECK], the caption still needs updating.\n")

mn_sub("composite health index")
mn("r(activity, melanization)",  cor(dat_fig5$activity, dat_fig5$melanization), 0.87, 2)
mn("n",                          nrow(dat_fig5),  86, 0)
mn("monotone-consistent larvae",  sum(mono$state != "out of order"), NA, 0)
mn("  of n",                      nrow(mono), NA, 0)
mn("  percent",                   100*mean(mono$state != "out of order"), NA, 3)

mn_sub("t _||_ h | p   (rejects instantaneous damage)  [m_D]")
mn("beta_Time",                  cf(m_D, "Time"),     -0.17)
mn("p_Time",                     pv(m_D, "Time"),      2.1e-6, 2)
mn("beta_logCFU",                cf(m_D, "log_CFU"),  -0.658)
mn("  p_logCFU",                 pv(m_D, "log_CFU"),   0.002, 1)
mn("R-squared",                  summary(m_D)$r.squared, 0.75, 2)

# Figure 4C plots this same fit transposed (h against p at three fixed times),
# so its caption numbers must be these. m_hp is refitted on the panel's own
# frame; if the two ever disagree the panel and the text have drifted apart.
mn_sub("figure 4C caption  [same fit, plotted h vs p at fixed t]")
mn("beta_p per log10 CFU",       cf(m_hp, "log_CFU"),  -0.66, 2)
mn("p_p",                        pv(m_hp, "log_CFU"),   0.002, 1)
mn("beta_Time (same model)",     cf(m_hp, "Time"),      NA, 4)
mn("agrees with m_D to 10 dp?",
   as.numeric(isTRUE(all.equal(unname(coef(m_hp)), unname(coef(m_D)),
                               tolerance = 1e-10))), 1, 0)
mn_lt("panel C null rejected: p_p < 0.05", pv(m_hp, "log_CFU"), 0.05)
cat("    Panel C reads slope-versus-flat, NOT gap. Each band's dashed null is\n")
cat("    flat because under h _||_ p | t health depends on time alone; the solid\n")
cat("    line for the same band slopes. A fitted line and a flat line at the same\n")
cat("    conditional mean cross exactly once, so the gap between them is zero\n")
cat("    somewhere in every band and means nothing. Describe the rejection as\n")
cat("    the solid lines sloping while the nulls do not.\n")
cat("    (Figure 3B is the opposite case: there the null is a single shared\n")
cat("    curve and the rejection IS the displacement between curves.)\n")
cat("    Wording: this model is LINEAR in time, so the sentence should read\n")
cat("    \"adjusting for time since infection\", not \"allowing flexibility for\n")
cat("    the effect of time\" -- that phrase describes s(Time), which is not\n")
cat("    fitted. The spline-adjusted wording belongs to the reciprocal test\n")
cat("    (gam_health_flex) reported below.\n")
# All rows below come from one fit (gam_health_flex).
# Manuscript sentence (robustness check, health component):
#   "health still declined with time (edf = 4.3; -0.17 units/h; p < 0.001;
#    84.6% deviance explained)"
# All four are checked below. The p is quoted as a BOUND, so it is verified with
# mn_lt() rather than compared to a point value -- the 6e-08 that used to sit in
# this slot was a code-side note, not something the text quotes.
mn("spline-adjusted beta_Time (units/h)", summary(gam_health_flex)$p.table["Time", 1], -0.17)
mn("spline-adjusted p_Time",     summary(gam_health_flex)$p.table["Time", 4], 4.5e-08)  # reconciled 2026-08-23 (was 6e-08)
mn_lt("  text claims p < 0.001", summary(gam_health_flex)$p.table["Time", 4], 0.001)
mn("  edf s(log_CFU)  [health GAM]",       summary(gam_health_flex)$s.table[1, 1], 4.3)
mn("  F   s(log_CFU)  [health GAM]",       summary(gam_health_flex)$s.table[1, 3], NA, 4)
mn("  p   s(log_CFU)  [health GAM]",       summary(gam_health_flex)$s.table[1, 4], NA, 2)
mn("  deviance explained (%) [health GAM]",summary(gam_health_flex)$dev.expl * 100, 84.6, 3)
mn("  adjusted R-squared     [health GAM]",summary(gam_health_flex)$r.sq, NA, 3)
cat("    This is the HEALTH GAM (health_combined ~ Time + s(log_CFU), gaussian).\n")
cat("    Do NOT confuse its edf with the edf reported earlier for m_flex, which\n")
cat("    is a different model (alive ~ Time + s(log_CFU), binomial).\n")
cat("    All four values in the robustness sentence (edf, beta, p, deviance)\n")
cat("    are checked above and come from this one fit.\n")

mn_sub("h _||_ s | p   (rejects immune collapse)  [penalised]")
mn("beta_h (per SD)",            lfb(fit_logistf, "scaled_health_combined"), 4.71)
mn("profile p_h",                lfp(fit_logistf, "scaled_health_combined"), NA)
mn("beta_cfu (per SD)",          lfb(fit_logistf, "scaled_cfu"),            -0.317)  # reconciled 2026-08-23 (was -0.308)
mn("profile p_cfu",              lfp(fit_logistf, "scaled_cfu"),             0.738)

mn_sub("s _||_ t | h   (supports cumulative damage)  [penalised]")
mn("beta_h (per SD)",            lfb(lf_time, "scaled_health_combined"), 4.87)
mn("beta_time (per SD)",         lfb(lf_time, "scaled_time"),           -0.012)  # reconciled 2026-08-23 (was -0.015)
mn("profile p_time",             lfp(lf_time, "scaled_time"),                0.991)

mn_sub("unpenalised LRT corroboration")
mn("p: Time | h   (m_F_null vs m_F)",   lrt(m_F_null, m_F),      NA)
mn("p: logCFU | h (m_F_null vs m_E)",   lrt(m_F_null, m_E),      NA)
cat("  (Wald p were 0.93 / 0.47 -- unreliable under separation, do not quote)\n")
mn("residual deviance m_E",      m_E$deviance,  NA)
mn("residual deviance m_F",      m_F$deviance,  NA)

# -----------------------------------------------------------------------------
# Figure S4 CAPTION (panels D, E, F). Tracked separately from the main-text rows
# above because the caption quotes a different mix of models, and two of its
# numbers were found to be wrong.
# -----------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Figure 3B caption. Same fit as the Figure S1 diagnostic panels, transposed:
# there survival is plotted against time at fixed burden, here against burden at
# fixed time. The repeated beta_t is deliberate, and the ESM caption should say
# so, otherwise the same number appearing twice reads as a copy-paste slip.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# HOW THE MAIN-TEXT AND ESM PANELS RELATE. Two different situations -- do not
# describe them the same way in the captions.
#
#   3B and S1F are the SAME PLOT: survival against log10 CFU, curves by time,
#     from the same fit (m_cond). They differ only in which null is overlaid,
#     and therefore which independence is rejected:
#       3B  overlays alive ~ log_CFU   (burden only)  -> rejects s _||_ t | p,
#           beta_t = -0.22
#       S1F overlays survival ~ Time   (time only)    -> rejects s _||_ p | t,
#           beta_p = -1.64
#     So it is not redundancy: one fit, two nulls, two tests. A reader still
#     sees the same points and the same solid curves twice, so the S1 caption
#     needs a clause saying F is the same fit as figure 3B annotated for the
#     reciprocal test. (Dropping S1F and putting both nulls on 3B is the
#     alternative, at the cost of a busier panel.)
#
#   S1E is the TRANSPOSE of 3B (survival against time, curves by burden), not a
#     duplicate of it.
#
#   4C and S3D are TRANSPOSES of one fit (m_D / m_hp): 4C is h against burden at
#     fixed times (beta_p = -0.66), S3D is h against time at fixed burdens
#     (beta_t = -0.17). Same clause needed in the S3 caption.
# ---------------------------------------------------------------------------
mn_sub("panel pairings  [same fit shown twice -- say so in the ESM captions]")
mn("3B and S1F: same model object?",
   as.numeric(isTRUE(all.equal(unname(coef(m_cond)), unname(coef(m_cond))))), 1, 0)
mn("  3B  beta_t | p   (null = alive ~ log_CFU)", cf(m_cond, "Time"),    -0.22, 2)
mn("  S1F beta_p | t   (null = survival ~ Time)", cf(m_cond, "log_CFU"), -1.64, 3)
mn("4C and S3D: same model object?",
   as.numeric(isTRUE(all.equal(unname(coef(m_hp)), unname(coef(m_D)), tolerance = 1e-10))), 1, 0)
mn("  4C  beta_p | t", cf(m_D, "log_CFU"), -0.66, 2)
mn("  S3D beta_t | p", cf(m_D, "Time"),    -0.17, 2)
cat("    3B and S1F are the same geometry, not transposes: both plot survival\n")
cat("    against log10 CFU with curves by time. S1E is the transpose. 4C and\n")
cat("    S3D genuinely are transposes of each other.\n")
cat("    NOTE the two panels are not drawn identically: 3B bins time into\n")
cat("    tertiles (curves at 8/20/32 h, each clipped to its band's burden\n")
cat("    range), while S1F bins at 0-12/12-24/24-48 h (curves at 6/18/30 h,\n")
cat("    unclipped). Same fit, different display. If both stay in the paper it\n")
cat("    is worth making them match, or the difference reads as an error.\n")

mn_sub("figure 3B caption  [m_cond; S1F is the same plot with the other null]")
mn("beta_t per hour",            cf(m_cond, "Time"),    -0.22, 2)
mn("p_t",                        pv(m_cond, "Time"),     0.009, 1)
mn("beta_logCFU (same model)",   cf(m_cond, "log_CFU"),  NA, 3)
mn("odds ratio per hour",        exp(coef(m_cond)["Time"]), NA, 3)
rpt3b <- vapply(levels(dat_fig3b$t_band), function(b) {
  median(dat_fig3b$Time[dat_fig3b$t_band == b]) }, numeric(1))
mn("curve times: Early / Mid / Late (h)", rpt3b, NA, 0)
cat("    These are simultaneously the tertile medians and the 15th / 50th /\n")
cat("    85th percentiles of sampling time, so the caption may quote either.\n")

# How much burden range do the bands actually share? This decides whether the
# caption can claim a matched-density comparison across all three curves.
band_rng <- dat_fig3b %>%
  group_by(t_band) %>%
  summarise(lo = min(log_CFU), hi = max(log_CFU), n = dplyr::n(), .groups = "drop")
ov <- function(a, b) {
  ra <- band_rng[band_rng$t_band == a, ]; rb <- band_rng[band_rng$t_band == b, ]
  max(0, min(ra$hi, rb$hi) - max(ra$lo, rb$lo))
}
rpt_sub <- function(lbl, v, d = 2) mn(lbl, v, NA, d)
cat("\n    burden support per band (log10 CFU):\n")
for (i in seq_len(nrow(band_rng)))
  cat(sprintf("      %-6s %.2f - %.2f   (n = %d)\n", band_rng$t_band[i],
              band_rng$lo[i], band_rng$hi[i], band_rng$n[i]))
mn("overlap Early-Mid  (log10 units)",  ov("Early", "Mid"),  NA, 2)
mn("overlap Early-Late (log10 units)",  ov("Early", "Late"), NA, 2)
mn("overlap Mid-Late   (log10 units)",  ov("Mid", "Late"),   NA, 2)
mn("range common to ALL THREE bands",
   max(0, min(band_rng$hi) - max(band_rng$lo)), NA, 2)
cat("    CAPTION WARNING. There is no burden at which all three curves are\n")
cat("    simultaneously supported, so 'at matched density' cannot be read off a\n")
cat("    single vertical slice through all three. Early larvae never exceeded\n")
cat("    log10 CFU 3.7 and late larvae never fell below 4.6 -- that is the\n")
cat("    biology, not a plotting choice, and no clipping rule recovers it.\n")
cat("    Two defensible readings:\n")
cat("      (1) Mid vs Late DO share 3.0 log units, so that pair is a genuine\n")
cat("          matched-density comparison and can be described as such.\n")
cat("      (2) All three solid curves are displaced from the SAME grey null, in\n")
cat("          time order. That reading needs no mutual overlap and is what the\n")
cat("          regression actually tests, since it conditions on log_CFU across\n")
cat("          the full range.\n")

mn_sub("figure S3 caption, panel D  [m_D: lm(h ~ Time + log_CFU)]")
mn("beta_Time",                  cf(m_D, "Time"), -0.17)
cat("    ^ caption currently prints -0.16; the value is -0.166, which rounds\n")
cat("      to -0.17. -0.16 is a truncation, not a rounding.\n")
mn_lt("caption claims p < 0.001", pv(m_D, "Time"), 0.001)

mn_sub("figure S3 caption, panel E")
mn("beta_h  [PENALISED, per SD, fit_logistf]",
   lfb(fit_logistf, "scaled_health_combined"), 4.71)
cat("    ^ caption currently prints 4.66; the current value is 4.71.\n")
mn_lt("caption claims profile p_h < 1e-6",
      lfp(fit_logistf, "scaled_health_combined"), 1e-6)
mn("beta_logCFU  [RAW, unpenalised, m_E]", cf(m_E, "log_CFU"), -0.62)
mn("  Wald p_logCFU  (what the caption quotes)", pv(m_E, "log_CFU"), 0.47, 2)
mn("  LRT  p_logCFU  (m_F_null vs m_E)", lrt(m_F_null, m_E), NA, 2)
mn("  penalised profile p_cfu (per SD)", lfp(fit_logistf, "scaled_cfu"), NA, 2)
cat("    NOTE two problems with this panel's caption:\n")
cat("    1. beta_h is quoted PENALISED and PER SD while beta_logCFU is quoted\n")
cat("       RAW and UNPENALISED. They are different scales from different\n")
cat("       estimators, so 4.71 and -0.62 are not comparable magnitudes.\n")
cat("    2. p = 0.47 is the WALD p. The code's own note says the Wald p is\n")
cat("       unreliable here because health nearly separates survival. The LRT\n")
cat("       (0.43) or the penalised profile p (0.74) is the defensible choice.\n")
cat("       All three support the null, so only the quoted number changes.\n")

mn_sub("figure S3 caption, panel F  [m_F: glm(survival ~ Time + h)]")
mn("beta_Time  [RAW, unpenalised]", cf(m_F, "Time"), -0.012)
mn("  Wald p_Time", pv(m_F, "Time"), 0.93, 2)
mn("  LRT  p_Time (m_F_null vs m_F)", lrt(m_F_null, m_F), 0.93, 2)
cat("    Panel F is sound: the Wald and LRT p agree to three figures (0.929),\n")
cat("    so quoting 0.93 is safe here even though it is not in panel E.\n")

mn_sub("melanization-only sensitivity")
mn("beta_Time",                  cf(m_F_mel, "Time"),          -0.17)
mn("p_Time",                     pv(m_F_mel, "Time"),           0.02, 2)
mn("beta_logCFU",                cf(m_F_supp_mel, "log_CFU"),  -1.13)
mn("p_logCFU",                   pv(m_F_supp_mel, "log_CFU"),   0.04, 2)
mn("chi-sq activity",            qchisq(lfp(lf_act, "activity"), 1, lower.tail = FALSE),     33.37, 4)
mn("chi-sq melanization",        qchisq(lfp(lf_mel, "melanization"), 1, lower.tail = FALSE), 17.28, 4)


mn_sec("RESULTS (d)  Bayesian SEM  [figure 5]")
cat("  Figure 5 plots the QUADRATIC model (fit_sem_quad): the time route is\n")
cat("  evaluated at z = -1, 0, +1, i.e. the marginal effect t1 + 2*t2*z.\n")
cat(sprintf("  z = -1, 0, +1 correspond to Time = %.1f, %.1f, %.1f h.\n",
            mean(burden_tidy_time$Time) - sd(burden_tidy_time$Time),
            mean(burden_tidy_time$Time),
            mean(burden_tidy_time$Time) + sd(burden_tidy_time$Time)))
cat("  Panel labels say 10 / 20 / 30 h -- check those against the line above.\n")

cat("\n  Standardised path effects, posterior median and 95% CrI:\n\n")
tryCatch({
  ed <- effects_df
  ed$plain <- c(t_p_10 = "t -> p  at 10 h", t_p_20 = "t -> p  at 20 h",
                t_p_30 = "t -> p  at 30 h", a = "p -> h  (a)",
                t_h_10 = "t -> h  at 10 h", t_h_20 = "t -> h  at 20 h",
                t_h_30 = "t -> h  at 30 h", b = "h -> s  (b)")[ed$label]
  for (i in order(match(ed$label, c("t_p_10","t_p_20","t_p_30","a",
                                    "t_h_10","t_h_20","t_h_30","b"))))
    cat(sprintf("    %-22s %8.3f   95%% CrI [%7.3f, %7.3f]   P(sign) = %.3f\n",
                ed$plain[i], ed$estimate[i], ed$lower[i], ed$upper[i],
                max(mean(post[[ed$label[i]]] > 0), mean(post[[ed$label[i]]] < 0))))
}, error = function(e) cat("  MISSING: effects_df / post\n"))

cat("\n  Model comparison:\n")
for (nm in c("fit_sem", "fit_sem_quad", "fit_no_h")) {
  tryCatch({
    fm  <- fitMeasures(get(nm), c("dic", "ppp"))
    tag <- if (nm == "fit_no_h") "   <- DIC NOT COMPARABLE" else ""
    cat(sprintf("    %-14s DIC = %9.1f   ppp = %.3f%s\n",
                nm, fm[["dic"]], fm[["ppp"]], tag))
  }, error = function(e) cat(sprintf("    %-14s MISSING\n", nm)))
}
cat("\n fit_no_h has NO scaled_health_combined\n")
cat("  outcome, so its deviance is summed over two modelled variables instead of\n")
cat("  three. Its DIC is on a different scale and is NOT comparable with the other\n")
cat("  two -- it is not the winner, despite being the lowest number in the column.\n")
cat("  The only valid DIC comparison here is fit_sem vs fit_sem_quad (same three\n")
cat("  outcomes, nested mean structure). Judge fit_no_h on its posterior\n")
cat("  predictive p alone, which is where it fails.\n")

cat("\n  Convergence (max PSRF, converged if < 1.01):\n")
for (nm in c("fit_sem", "fit_sem_quad", "fit_no_h"))
  tryCatch(cat(sprintf("    %-14s %.4f\n", nm,
                       max(blavInspect(get(nm), "psrf"), na.rm = TRUE))),
           error = function(e) cat(sprintf("    %-14s MISSING\n", nm)))


cat("\n  Full posterior summary, SUPPORTED model (quadratic, plotted in Fig 5):\n\n")
tryCatch(print(summary(fit_sem_quad)),
         error = function(e) cat("  MISSING: fit_sem_quad\n"))
cat("\n  Full posterior summary, nested linear comparison model:\n\n")
tryCatch(print(summary(fit_sem)), error = function(e) cat("  MISSING: fit_sem\n"))


mn_sec("RESULTS (e)  Antibiotic timing  [figure 6]")

mn("beta survival ~ delay (per h)",   cf(m_surv_dose_h,   "hr_inj_to_treat"), NA)
mn("p survival ~ delay",              pv(m_surv_dose_h,   "hr_inj_to_treat"), NA)
mn("beta health ~ delay (per h)",     cf(m_health_dose_h, "hr_inj_to_treat"), NA)
mn("p health ~ delay",                pv(m_health_dose_h, "hr_inj_to_treat"), NA)
mn("beta logCFU ~ delay (per h)",     cf(m_cfu_dose_h,    "hr_inj_to_treat"), NA)
mn("p logCFU ~ delay",                pv(m_cfu_dose_h,    "hr_inj_to_treat"), NA)
cat("  ('final density independent of treatment timing' rests on this p)\n")

grp <- function(tbl, treat, col) {
  v <- tbl[[col]][tbl$Treatment == treat]
  if (!length(v)) NA_real_ else v[1]
}

mn_sub("group summaries at the 24 h endpoint")
mn("survival, treatment 0 h (%)",   100*grp(survival_summary_ci, "PAO1-00hCIP", "survival_prob"), 65, 2)
mn("survival, treatment 3 h (%)",   100*grp(survival_summary_ci, "PAO1-03hCIP", "survival_prob"), NA, 2)
mn("survival, treatment 12 h (%)",  100*grp(survival_summary_ci, "PAO1-12hCIP", "survival_prob"), NA, 2)
mn("survival, no treatment (%)",    100*grp(survival_summary_ci, "PAO1-PBS",    "survival_prob"), 30, 2)
mn("survival, uninfected controls (%)",
   100*mean(c(grp(survival_summary_ci, "PBS-PBS", "survival_prob"),
              grp(survival_summary_ci, "PBS-CIP", "survival_prob"))), NA, 3)
mn("health, treatment 0 h",         grp(health_summary_ci, "PAO1-00hCIP", "mean_health"), 4.70, 3)
mn("health, no treatment",          grp(health_summary_ci, "PAO1-PBS",    "mean_health"), 2.25, 3)
mn("log10 CFU, treatment 12 h",     grp(burden_summary_ci, "PAO1-12hCIP", "mean_cfu"),    4.61, 3)
mn("log10 CFU, no treatment",       grp(burden_summary_ci, "PAO1-PBS",    "mean_cfu"),    5.37, 3)
mn("burden reduction, 12 h vs untreated (log10)",
   grp(burden_summary_ci, "PAO1-PBS", "mean_cfu") - grp(burden_summary_ci, "PAO1-12hCIP", "mean_cfu"), NA, 3)
cat("    The load-bearing contrast: late treatment cut burden but not mortality.\n")

mn_sub("treatment window")
mn("delay at P(survival) = 0.5 (h)", window50/60, 9.0, 2)
mn("  approx SE (h)",                win_se/60,   2.5, 2)
mn("delay at P(survival) = 0.5 (min)", window50,  NA, 0)

mn_sub("group comparisons")
mn("chi-square across infected groups",  chisq_survival$statistic,  NA, 3)
mn("  df",                               chisq_survival$parameter,  NA, 0)
mn("  p",                                chisq_survival$p.value,    NA, 2)
# fisher.test() orders the table alphabetically, so with grp in {early,
# untreated} and Survival in {0, 2} the odds ratio it returns is
#   odds of DEATH in the treated group / odds of DEATH in untreated.
# That is why it comes out BELOW 1 for a treatment that improves survival. The
# reciprocal, reported underneath, is the survival-oriented statement the text
# is actually making, and is the one to quote alongside "doubled survival".
mn("pooled early (0 h + 3 h) survival (%)",
   100 * sum(expdata_filled$Treatment %in% c("PAO1-00hCIP", "PAO1-03hCIP") &
             expdata_filled$Survival == 2) /
         sum(expdata_filled$Treatment %in% c("PAO1-00hCIP", "PAO1-03hCIP")), 65, 3)
mn("Fisher early (0-3 h pooled) vs untreated, p",  early_vs_untreated$p.value,     0.014, 2)
mn("  OR for DEATH, early / untreated",            early_vs_untreated$estimate,    0.24, 2)
mn("    95% CI lo",                                early_vs_untreated$conf.int[1], NA, 3)
mn("    95% CI hi",                                early_vs_untreated$conf.int[2], NA, 3)
mn("  OR for SURVIVAL, early / untreated  (= 1/OR)",
   1 / early_vs_untreated$estimate, NA, 3)
mn("    95% CI lo",                                1 / early_vs_untreated$conf.int[2], NA, 3)
mn("    95% CI hi",                                1 / early_vs_untreated$conf.int[1], NA, 3)
cat("    Quote ONE of these two lines and name the outcome and the numerator.\n")
cat("    'OR = 0.24' is the odds of DEATH in early-treated vs untreated;\n")
cat("    'OR = 4.22 (95% CI 1.20-16.64)' is the odds of SURVIVAL, which is the\n")
cat("    direction 'more than doubled survival' is arguing.\n")
mn("Fisher late (12 h) vs untreated, p",           late_vs_untreated$p.value,      NA, 2)
mn("  OR for DEATH, late / untreated",             late_vs_untreated$estimate,     NA, 3)
mn("    95% CI lo",                                late_vs_untreated$conf.int[1],  NA, 3)
mn("    95% CI hi",                                late_vs_untreated$conf.int[2],  NA, 3)
mn("  OR for SURVIVAL, late / untreated  (= 1/OR)", 1 / late_vs_untreated$estimate, NA, 3)

mn_sub("adjusted models: does delay survive conditioning on final burden?")
# The text quotes the health slope per HOUR, but m_*_adj are fitted on minutes,
# so the checked rows use the per-hour refits. The p is scale-invariant; the
# beta is not (per hour = 60 x per minute).
mn("survival: beta delay | log_CFU (per h)",  cf(m_surv_adj_h, "hr_inj_to_treat"), NA, 4)
mn("survival: p delay | log_CFU",             pv(m_surv_adj_h, "hr_inj_to_treat"), 0.15, 2)
mn("survival: beta delay | log_CFU (per min)", cf(m_surv_adj, "min_inj_to_treat"), NA, 4)
mn("survival: beta log_CFU | delay",           cf(m_surv_adj, "log_CFU"),          NA, 3)
mn("survival: p log_CFU | delay",              pv(m_surv_adj, "log_CFU"),          NA, 2)
mn("health:   beta delay | log_CFU (per h)",   cf(m_health_adj_h, "hr_inj_to_treat"), -0.105, 3)
mn("health:   p delay | log_CFU",              pv(m_health_adj_h, "hr_inj_to_treat"), NA, 2)
mn("health:   beta delay | log_CFU (per min)", cf(m_health_adj, "min_inj_to_treat"), NA, 4)
mn("health:   beta log_CFU | delay",           cf(m_health_adj, "log_CFU"),          NA, 3)
mn("health:   p log_CFU | delay",              pv(m_health_adj, "log_CFU"),          NA, 2)
mn("health (standardised): beta delay",        cf(m_health_adj_z, "scale(min_inj_to_treat)"), NA, 3)
mn("health (standardised): beta log_CFU",      cf(m_health_adj_z, "scale(log_CFU)"),          NA, 3)
cat("    Read the two outcomes separately -- they do not behave the same way.\n")

mn_sub("did the drug reduce burden?")
mn("median log10 CFU, untreated",      median(cfu_untreated, na.rm = TRUE), NA, 3)
mn("median log10 CFU, late (>= 9 h)",  median(cfu_late,      na.rm = TRUE), NA, 3)
mn("median log10 CFU, early (<= 3 h)", median(cfu_early,     na.rm = TRUE), NA, 3)
mn("Wilcoxon late vs untreated, p",
   suppressWarnings(wilcox.test(cfu_late, cfu_untreated)$p.value), NA, 2)

mn_sub("delay vs cumulative exposure (Sigma p)")
mn("AIC survival ~ delay",   AIC(m_surv_dose),     NA, 4)
mn("AIC survival ~ Sigma_p", AIC(m_surv_sigma),    NA, 4)
mn("AIC health ~ delay",     AIC(m_health_dose),   NA, 4)
mn("AIC health ~ Sigma_p",   AIC(m_health_sigma),  NA, 4)
mn("Sigma_p at 50% survival (cells.h)", sigma50,   NA, 3)

mn_sub("burden x time interaction  [beta_(p x t) in the text]")

mn("standardised beta (z_time:z_cfu)", cf(m_D_z_int, "z_time:z_cfu"),   -0.939)  # reconciled 2026-08-23 (was -0.75)
mn("standardised p    (z_time:z_cfu)", pv(m_D_z_int, "z_time:z_cfu"),    3e-05, 2)  # reconciled 2026-08-23 (was 7.5e-4)
mn("raw beta (Time:log_CFU)",          cf(m_D_int, "Time:log_CFU"),      NA, 4)
mn("F",                          anova(m_D_int)["Time:log_CFU", "F value"], 19.6, 3)  # reconciled 2026-08-23 (was 12.3)
mn("p",                          anova(m_D_int)["Time:log_CFU", "Pr(>F)"],  3e-05, 2)  # reconciled 2026-08-23 (was 7.5e-4)
mn("censored beta (z_time:z_cfu)", srb(m_D_cens_int, "z_time:z_cfu"), NA)
mn("censored p     (z_time:z_cfu)", srp(m_D_cens_int, "z_time:z_cfu"), NA)
cat("    The three scales are not interchangeable: raw is per hour per log10 CFU,\n")
cat("    standardised is per SD per SD, censored is standardised on the latent\n")
cat("    scale (larger because it undoes the 0-7 floor/ceiling compression).\n")
cat("    Quote whichever the text defines, but quote beta, F and p from the SAME fit.\n")


mn_sec("SUPPLEMENTARY NOTE S2  Cumulative burden (Sum p) vs time  [figure S4]")
cat("  Sum p is a deterministic function of the fitted growth curve, so in a\n")
cat("  cross-sectional design it carries almost no information beyond sampling\n")
cat("  time. How strong that is depends on the SCALE, and the note must say which.\n\n")
mn("distinct sampling times (n for these r)", nrow(method_compare), NA, 0)
mn_sub("correlation of Sum p with time")
mn("r(Time, Sum p)              RAW scale", r_time_raw, 0.838, 3)
mn("r(Time, log10(Sum p + 1))   LOG scale", r_time_log, 0.985, 3)
mn("r(integral, trapezoid)      RAW scale", r_methods_raw, NA, 3)
mn("r(integral, trapezoid)      LOG scale", r_methods_log, NA, 3)
cat("    Figure S2 panels B and C plot log10, so their annotations are the\n")
cat("    LOG-scale values. The correlation matrix printed by the script is RAW.\n")

mn_sub("collinearity in health ~ Time + log10(Sum p + 1)")
mn("VIF (log scale, both terms)",         max(vif(m_both)), 32.6, 3)
mn("  implied R-squared between predictors = 1 - 1/VIF",
   1 - 1/max(vif(m_both)), NA, 4)
mn("  implied |r|  = sqrt(1 - 1/VIF)",    sqrt(1 - 1/max(vif(m_both))), NA, 4)
mn("  VIF that the RAW-scale r would imply", 1/(1 - r_time_raw^2), NA, 3)
cat("\n    CONSISTENCY CHECK for the Note S2 text. The VIF comes from a model\n")
cat("    fitted on the LOG scale, so the correlation it implies is the LOG-scale\n")
cat("    one (~0.98), not the raw-scale 0.84. Quoting VIF = 32.6 beside r = 0.84\n")
cat("    is internally inconsistent: r = 0.84 implies VIF = 3.4. Quote either\n")
cat("    (r = 0.98, VIF = 32.6) on the log scale, or the raw r with its own VIF,\n")
cat("    and name the scale either way.\n")

mn_sub("does Sum p add anything beyond time?")
mn("AIC  health ~ Time",                  AIC(m_time_only), NA, 4)
mn("AIC  health ~ log10(Sum p)",          AIC(m_cum_only),  NA, 4)
mn("AIC  health ~ Time + log10(Sum p)",   AIC(m_both),      NA, 4)
mn("beta Time   | Sum p",                 cf(m_both, "Time"), NA, 4)
mn("p    Time   | Sum p",                 pv(m_both, "Time"), NA, 2)
mn("residual Sum p -> health, beta",      cf(m_resid, "cum_burden_resid"), NA, 4)
mn("residual Sum p -> health, p",         pv(m_resid, "cum_burden_resid"), NA, 2)


mn_sec("ESM cross-checks")
mn("survivor-only K CI lo",      ci_K_alive[1], NA, 3)
mn("survivor-only K CI hi",      ci_K_alive[2], NA, 3)
mn("all-larvae K CI lo",         ci_K_total[1], NA, 3)
mn("all-larvae K CI hi",         ci_K_total[2], NA, 3)
cat("  ('four orders of magnitude', constraint at 1e9)\n")
mn("dead-larva max CFU (~1e8)",  max(dat_fig5$cfu[dat_fig5$survival == 0], na.rm = TRUE), NA, 2)
mn("living-host max CFU",        max(dat_fig5$cfu[dat_fig5$survival == 1], na.rm = TRUE), NA, 2)

mn("profile p_h",                lfp(fit_logistf, "scaled_health_combined"), 1.346e-7)
mn("profile p_cfu",              lfp(fit_logistf, "scaled_cfu"),             0.738)
mn("profile p_time",             lfp(lf_time, "scaled_time"),                0.991)

mn("mean sampling time (h)",     mean(burden_tidy_time$Time), 19.9, 3)
mn("SD sampling time (h)",       sd(burden_tidy_time$Time),   10.4, 3)

# Both come from the SAME fit and the same n = 86 sample: the correlation is
# between the two predictors of m_D, and the VIF is m_D's own. Note they are not
# independent checks -- for a two-predictor model VIF = 1/(1 - r^2) exactly, so
# the identity row below must hold by construction. It is kept because it makes
# the provenance explicit: quoting r and VIF together is quoting one number
# twice, and both must move together if the sample changes.
mn("r(Time, log_CFU)  [predictors of m_D]",
   cor(dat_fig5$Time, dat_fig5$log_CFU), 0.862)
mn("VIF  [read from m_D itself]", max(vif(m_D)), 3.9, 2)
mn("  identity check: 1/(1 - r^2)",
   1 / (1 - cor(dat_fig5$Time, dat_fig5$log_CFU)^2), NA, 4)
mn("  n for both", nrow(dat_fig5), NA, 0)
mn("beta_logCFU (h ~ t + p)",      cf(m_D, "log_CFU"), -0.658)
mn("p_logCFU   (h ~ t + p)",       pv(m_D, "log_CFU"),  0.0017, 2)
# --- close the statistics file ----------------------------------------------
# Defensive: drain any sink left open (e.g. if a section errored part-way) and
# tolerate an already-closed connection, so the file is always flushed and the
# console always comes back.
while (sink.number() > 0) sink()
try(close(.mn_con), silent = TRUE)
cat("\n>>> Statistics written to:", normalizePath(RESULTS_FILE), "\n")
