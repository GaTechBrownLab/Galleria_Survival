# Analysis code for "Host health, not instantaneous pathogen burden, 
# determines survival during acute Pseudomonas ae-ruginosa infection"
# Author: Canan Karakoç

# Load required libraries
library(minpack.lm)  # For non-linear regression
library(dplyr)
library(tidyr)
library(boot)
library(MASS)
library(mgcv)

# Plotting 
library(ggplot2)  
library(cowplot)
library(patchwork)


# SEM
library(lavaan)
library(blavaan) #bayesian lavaan
library(lavaanPlot)
library(semPlot)
library(brglm2) #biased-reduced regression 
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
library(emmeans)

library(purrr) 

# Theme 
# panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
mytheme <- theme_bw() +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14, color = "black"), axis.title = element_text(size = 16)) +
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    size = 1
  )) +
  theme(strip.text.x = element_text(size = 14), strip.background = element_blank()) +
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

# Working directory, change with yours.
setwd("~/Documents/GitHub/Galleria_Survival")

survival      <- read.table("data/AliveDead.csv", header = T, sep = ",", dec =".")
burden        <- read.table("data/bacterial_burden.csv", header = T, sep = ",", dec =".")
health        <- read.table("data/health_assesment.csv", header = T, sep = ",", dec =".")
death         <- read.table("data/time_to_death.csv", header = T, sep = ",", dec =".")
control_surv  <- read.table("data/control_survival.csv", header = T, sep = ",", dec =".")

# Common time sequence used throughout
time_seq <- seq(0, 36, length.out = 200)

#===============================================================================
# FIGURE 1
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
) %>%
  mutate(
    # Normalize AFTER filtering to display range
    mort_norm = (mortality - min(mortality)) / (max(mortality) - min(mortality))
  )

# Create inset plot
inset_plot <- ggplot(mort_fitted_df, aes(x = time, y = mortality)) +
  geom_line(color = "#b80422", linewidth = 1) +
  scale_x_continuous(breaks = c(0, 18, 36)) +
  scale_y_continuous(breaks = c(0, 35, 70)) +
  labs(x = "Time (h)", y = (bquote('m(t) '(h^-1)))) +
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
  geom_point(data = control, aes(fill = group), size = 2, shape = 21) +
  geom_point(data = plot_data, aes(fill = group), size = 2, shape = 21) +
  geom_line(data = gompertz_fit, aes(y = fitted_survival, linetype = group, color = group), linewidth = 1) +
  annotate("segment", x = 1, xend = 36, y = 1, yend = 1, 
           color = "grey50", linewidth = 0.8, linetype = "dashed") +
  labs(x = "Time (h)", y = "Proportion of survival") +
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
    xmin = -0.05, xmax = 19, ymin = -0.05, ymax = 0.45
  )

#===============================================================================
# SCENARIO PANELS - Using Empirical Parameters
#===============================================================================

# Extract BOTH Gompertz (from survival) and Logistic (from burden) parameters
params_gompertz <- coef(fit)
a <- params_gompertz["a"]
b <- params_gompertz["b"]

K_estimated  <- 1e9
p0_empirical <- 2e3  # Injected concentration
r_empirical  <- b    # assuming survival rate is equal to pathogen growth rate 

# Time sequence
t_seq <- seq(0, 36, length.out = 100)

# Normalize function
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

#-------------------------------------------------------------------------------
# SCENARIO 1: Exponential growth + Linear m(p)
#-------------------------------------------------------------------------------
# In the manuscript used i instead of k to avoid confusion with K
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
) %>%
  mutate(
    time_norm = normalize(time),
    p_norm = normalize(p),
    m_norm = normalize(m_t)
  )

scientific_10 <- function(x) {
  ifelse(x == 0, "0", 
         parse(text = gsub("e\\+?", " %*% 10^", scales::scientific(x))))
}

p1C <- ggplot(df_s1, aes(x = time, y = p)) +
  geom_line(color = "#19798b", linewidth = 1.5) +
  labs(y = "Pathogen p(t)", x = "Time") +
  scale_x_continuous(breaks = c(0, 12, 24, 36)) +
  scale_y_continuous(labels = scientific_10)+
  mytheme +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p1D <- ggplot(df_s1, aes(x = p, y = m_t)) +
  geom_line(color = "#e76f51", linewidth = 1.5) +
  labs(x = "Pathogen, p (CFU)", y = "Mortality m(p)") +
  scale_x_continuous(labels = scientific_10)+
  scale_y_continuous(breaks = c(0, 35, 70))+
  mytheme + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45))

# Combine all panels
scenario1 <- (p1C / p1D) + plot_layout(heights = c(1, 1))

figure1 <- (pA_with_inset | scenario1) +
  plot_layout(widths = c(3, 1))+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 14)) 

figure1

# Save
ggsave("figures/figure1.pdf",  # >>> Manuscript Figure 1 (survival + Gompertz scenarios) 
       plot = figure1, width = 9, height = 6, units = "in", dpi = 300)

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
burden_tidy <- burden %>%
  pivot_longer(cols = rep1:rep4, names_to = "replicates", values_to = "cfu") %>%
  filter(! Countable == "no") %>%
  filter(!(Sample == "S11" & Larvae %in% c("L5", "L6"))) %>%
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
    cfu = ifelse(is.nan(cfu), 0, cfu),  
    survival = case_when(
      survival == 2 ~ 1,  
      survival == 0 ~ 0),  
    status = case_when(
      survival == 1 ~ "Alive",  
      survival == 0 ~ "Dead",    
      TRUE ~ NA_character_  
    )
  ) %>%
  mutate(
  ) %>%
  # ---- Melanization orientation (read this before touching anything) -------
# RAW health_assesment.csv coding: 1 = fully black (worst) ... 4 = pale/
# healthy (best). The RAW scale is therefore HEALTH-oriented (higher = better).
# The line below reorients it to SEVERITY (higher = worse) ONLY so the
# standalone melanization plot (Fig 4B) reads intuitively. After this line,
# `melanization` means 4 - raw  ->  0 = healthy ... 3 = black.
# Every health COMPOSITE below writes (4 - melanization), which UNDOES this
# reorientation and recovers the raw value (4 = healthy). So each composite is
# activity + (raw melanization): higher = healthier, lower = worse health.
mutate(melanization = 4 - melanization,   # now SEVERITY: higher = worse (use in plots)
       log_CFU = log10(cfu+1), 
       scaled_time = scale(Time), 
       scaled_health = scale(activity + (4 - melanization)),  # (4 - melanization) = raw (4 = healthy) -> higher = healthier
       scaled_activity = scale(activity), 
       scaled_melanization = scale(melanization), 
       scaled_cfu = scale(log_CFU),
       scaled_immune_morb = scale(activity + melanization)
) %>%
  filter(log_CFU > 0) %>% 
  ungroup() 

burden_tidy$Time <- as.numeric(burden_tidy$Time) 

#==============================================================================
# FIT MODELS
#==============================================================================

# remove final time point (larvae are dead for a long time)
burden_tidy_time <- burden_tidy %>%
  filter(Time < 37) %>%
  ungroup()

# only alive proportion
burden_tidy_alive <- burden_tidy %>% 
  filter(status == "Alive")

# Summary for starting values
burden_tidy_time_av <- burden_tidy_time %>% 
  group_by(Time) %>% 
  reframe(cfu = mean(cfu, na.rm = TRUE), log_cfu = mean(log_CFU))

burden_tidy_alive_av <- burden_tidy_alive %>% 
  group_by(Time) %>% 
  reframe(cfu = mean(cfu, na.rm = TRUE), log_cfu = mean(log_CFU))

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
coef(logistic_logfit_full)

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
                                start = list(
                                  K = quantile(burden_tidy_alive$cfu, 0.9, na.rm = TRUE),
                                  p0 = quantile(burden_tidy_alive$cfu, 0.1, na.rm = TRUE),
                                  r = 0.2
                                ),
                                lower = c(K = 100, p0 = 0.1, r = 0.01),
                                upper = c(K = 1e9, p0 = 5e3, r = 2),
                                control = list(maxiter = 1000))

# Bootstrap for alive only
pred_time_a_full <- data.frame(Time = seq(min(burden_tidy_alive$Time), 
                                          max(burden_tidy_alive$Time), 
                                          by = 0.1))

boot_preds_alive <- function(data, indices) {
  d <- data[indices, ]
  fit <- tryCatch({
    nlsLM(log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
          data = d,
          start = list(
            K = quantile(d$cfu, 0.9, na.rm = TRUE),
            p0 = quantile(d$cfu, 0.1, na.rm = TRUE),
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

# Compare K values
ci_K_total <- quantile(boot_curve_full$t[, 1], probs = c(0.025, 0.975), na.rm = TRUE)
ci_K_alive <- quantile(boot_curve_a_full$t[, 1], probs = c(0.025, 0.975), na.rm = TRUE)

#===============================================================================
# FIGURE 2
#===============================================================================

# Define colors and linetypes
line_colors <- c("Total population" = "#b80422", "Alive only" = "#19798b")
line_types  <- c("Total population" = "solid", "Alive only" = "dashed")

figure2 <- ggplot(burden_tidy_time, aes(x = Time, y = log_CFU)) +
  geom_jitter(aes(fill = status), size = 3, shape = 21, alpha = 0.6, color = "black") +
  # Total population fit
  geom_ribbon(data = pred_time_full, aes(x = Time, ymin = lower, ymax = upper), 
              fill = "#ee9b43", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = pred_time_full, aes(x = Time, y = fit, 
                                       color = "Total population", 
                                       linetype = "Total population"), 
            size = 1.2, inherit.aes = FALSE) +
  # Alive only fit
  geom_ribbon(data = pred_time_a_full, aes(x = Time, ymin = lower, ymax = upper), 
              fill = "#19798b", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = pred_time_a_full, aes(x = Time, y = fit, 
                                         color = "Alive only", 
                                         linetype = "Alive only"), 
            size = 1.2, inherit.aes = FALSE) +
  xlab("Time (h)") +
  ylab(bquote(log[10](CFU))) +
  scale_fill_manual(
    name = "Status",
    values = c("Alive" = "#19798b", "Dead" = "#ee9b43")
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
  theme(legend.position = c(0.2, 0.8)) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  )

ggsave("figures/figure2.pdf",  # >>> Manuscript Figure 2 (logistic growth) 
       plot = figure2, width = 6, height = 5, units = "in", dpi = 300)


#------------------------------------------------------------------------------
# Comparisons with other models 
#------------------------------------------------------------------------------
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

cat("  Logistic wins:", sum(winners == 1), "(", round(100 * mean(winners == 1)), "%)\n")
cat("  Exponential wins:", sum(winners == 2), "(", round(100 * mean(winners == 2)), "%)\n")
cat("  Linear wins:", sum(winners == 3), "(", round(100 * mean(winners == 3)), "%)\n")

winners_a <- boot_compare_alive$t[, 4]
winners_a <- winners_a[!is.na(winners_a)]

cat("  Logistic wins:", sum(winners_a == 1), "(", round(100 * mean(winners_a == 1)), "%)\n")
cat("  Exponential wins:", sum(winners_a == 2), "(", round(100 * mean(winners_a == 2)), "%)\n")
cat("  Linear wins:", sum(winners_a == 3), "(", round(100 * mean(winners_a == 3)), "%)\n")


#-------------------------------------------------------------------------------
# CFU "time since death"
# Supplementary figure S1
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
burden_tidy_death %>%
  filter(status == "Dead") %>%
  dplyr::select(Time, time_to_death_h, time_since_death, log_CFU) %>%
  summary()

# Plot: CFU vs time since death (dead larvae only)
p_postmortem <- ggplot(burden_tidy_death %>% filter(status == "Dead"), 
                       aes(x = time_since_death, y = log_CFU)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(
    x = "Time since death (h)",
    y = expression(log[10](CFU))
  ) +
  mytheme

print(p_postmortem)

# Save
ggsave("figures/figureS1.pdf",  # >>> Supplementary Figure S1 (post-mortem burden) 
       plot = p_postmortem, width = 7, height = 6, units = "in", dpi = 300)

# Test for post-mortem growth
lm_postmortem <- lm(log_CFU ~ time_since_death, 
                    data = burden_tidy_death %>% filter(status == "Dead"))
summary(lm_postmortem)

# Compare to time since infection for context
p_comparison <- burden_tidy_death %>%
  filter(status == "Dead") %>%
  ggplot(aes(x = Time, y = log_CFU)) +
  geom_point(aes(color = time_since_death), size = 2) +
  scale_color_viridis_c(name = "Hours\npost-death") +
  labs(
    x = "Time since infection (h)",
    y = expression(log[10](CFU))
  ) +
  mytheme

print(p_comparison)

#-------------------------------------------------------------------------------
# m(p) mapping - Excluded from the main text, but included in the 
# supplement to show how m(p) behaves as p approaches K
#-------------------------------------------------------------------------------

params_gompertz <- coef(fit)
a <- params_gompertz["a"]
b <- params_gompertz["b"]

# Real fitted logistic-growth parameters (raw-CFU scale). logistic_logfit_full
# is fitted above; without pulling K/p0/r from it here, this figure falls back
# on whatever values happen to linger in the session — which is why the blow-up
# landed near log10(CFU) = 1 instead of the real carrying capacity.
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

figure_mp <- ggplot(df_m, aes(x = p_log, y = m_p)) +
  geom_line(color = "#e76f51", linewidth = 1.5) +
  labs(x = bquote(log[10](CFU)), 
       y = bquote(m(p)~(h^-1))) +
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
params_logistic <- coef(logistic_logfit_full)
K  <- params_logistic["K"]; p0 <- params_logistic["p0"]; r <- params_logistic["r"]
a  <- coef(fit)["a"];        b  <- coef(fit)["b"]

expo <- 1                       # clean schematic; set expo <- b / r for the exact exponent
m_signed <- function(p) { R <- (p*(K - p0))/(p0*(K - p)); a * sign(R) * abs(R)^expo }

df_pole <- rbind(
  data.frame(pk = seq(0.02, 0.999, length.out = 500), side = "below"),
  data.frame(pk = seq(1.001, 2.5,  length.out = 500), side = "above")
)
df_pole$m <- m_signed(df_pole$pk * K)
ywin <- max(abs(m_signed(0.8 * K)), abs(m_signed(2 * K))) * 1.1

figure_s3 <- ggplot(df_pole, aes(pk, m, group = side)) +
  annotate("rect", xmin = 1, xmax = 2.5, ymin = -ywin, ymax = ywin, fill = "grey50", alpha = 0.08) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_line(colour = "#e76f51", linewidth = 1.3) +
  annotate("text", x = 1, y = ywin, label = "p == K", parse = TRUE, vjust = -0.4, size = 4) +
  annotate("text", x = 1.75, y = 0.82 * ywin, label = "p > K\n(unreachable)", size = 3.4, lineheight = 0.9) +
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("0", "K", "2K")) +
  labs(x = "Pathogen burden (p)", y = bquote(m(p)~(h^-1))) +
  coord_cartesian(ylim = c(-ywin*4, ywin*4), clip = "on") +
  mytheme + theme(panel.grid = element_blank(), plot.margin = margin(14, 10, 6, 6))

ggsave("figures/figureS3.pdf",
       plot = figure_s3, width = 5, height = 4, units = "in", dpi = 300)

#===============================================================================
# Health variables 
#===============================================================================

# Remove the 48h sample
dynamic_period <- burden_tidy %>% filter (Time < 37) 

# ==============================================================================
# Figure 4 Binomial Gompertz + hysteresis
# ------------------------------------------------------------------------------
# Layout:
#   A (top, full width):  Three Gompertz survival/decline curves with T50 markers
#                         (Activity, Melanization, Survival) — fit by binomial
#                         MLE on individual-level current-status data.
#   B (bottom-left):      Activity vs log10(CFU), points coloured by Time.
#   C (bottom-right):     Melanization vs log10(CFU), points coloured by Time.

# ============================================================================

# ----------------------------------------------------------------------------
# Prepare current-status data
# ----------------------------------------------------------------------------
# Definitions (matching existing code conventions):
#   active:    activity score >= 2     (raw 0-3 scale; not flipped)
#   unmelan:   melanization <= 1       (post-flip; raw "no melanization")
#   alive:     status == "Alive"
# Each row = one larva at its sampling time.

dat_cs <- burden_tidy %>%
  filter(Time < 37) %>%
  filter(!is.na(activity), !is.na(melanization), !is.na(status)) %>%
  transmute(
    time   = Time,
    active = as.integer(activity   >= 2),
    unmel  = as.integer(melanization <= 1),
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
  metric  = c("Activity (AT50)", "Melanization (MT50)", "Survival (LT50)"),
  T50     = c(T50_act, T50_mel, T50_surv),
  CI_low  = c(quantile(boot_AT, 0.025), quantile(boot_MT, 0.025), quantile(boot_LT, 0.025)),
  CI_high = c(quantile(boot_AT, 0.975), quantile(boot_MT, 0.975), quantile(boot_LT, 0.975))
)
print(T50_summary)

# Pairwise ordering probabilities (paired bootstraps)
n_paired <- min(length(boot_AT), length(boot_MT), length(boot_LT))
set.seed(456)
ix <- sample.int(n_paired)
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
  make_curve(fit_AT, "Activity (AT50)"),
  make_curve(fit_MT, "Melanization (MT50)"),
  make_curve(fit_LT, "Survival (LT50)")
)

# Observed proportions per timepoint for diagnostic overlay
obs_props <- dat_cs %>%
  pivot_longer(c(active, unmel, alive),
               names_to = "metric_short", values_to = "y") %>%
  group_by(time, metric_short) %>%
  summarise(p = mean(y), n = n(), .groups = "drop") %>%
  mutate(metric = recode(metric_short,
                         active = "Activity (AT50)",
                         unmel  = "Melanization (MT50)",
                         alive  = "Survival (LT50)"))

palette_T50 <- c(
  "Activity (AT50)"     = "#19798b",
  "Melanization (MT50)" = "#ee9b43",
  "Survival (LT50)"     = "#b80422"
)


# ----------------------------------------------------------------------------
# Panel C: T50 ordering
# ----------------------------------------------------------------------------
p4C <- ggplot(curves, aes(x = time, y = p, color = metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(data = obs_props, aes(x = time, y = p, color = metric, size = n),
             alpha = 0.7, show.legend = FALSE) +
  geom_segment(data = T50_summary,
               aes(x = T50, xend = T50, y = 0, yend = 0.5, color = metric),
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.4) +
  scale_color_manual(values = palette_T50, name = NULL) +
  scale_size_continuous(range = c(2, 4)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(limits = c(0, 36), breaks = seq(0, 36, 6)) +
  labs(x = "Time (h)", y = "Proportion") +
  mytheme +
  theme(legend.position = c(0.22, 0.28),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.text = element_text(size = 11))

# ----------------------------------------------------------------------------
# Panels B, C: Hysteresis (CFU vs health, time-coloured)
# ----------------------------------------------------------------------------
# Use viridis 'magma' reversed: early = light, late = dark.
# Black SCAM line gives the marginal CFU-only fit for visual reference.

# Fit GAM for CFU vs Activity relationship
gam_cfu_act <- scam(activity ~ s(log_CFU, k = 8, bs = "mpd"), data = dynamic_period)

# Generate predictions for smooth line
cfu_seq_act <- data.frame(log_CFU = seq(min(dynamic_period$log_CFU),
                                        max(dynamic_period$log_CFU),
                                        length.out = 100))
cfu_seq_act$smoothed_activity <- predict(gam_cfu_act, newdata = cfu_seq_act)

# Fit GAM for CFU vs Melanization relationship
gam_cfu_mel <- scam(melanization ~ s(log_CFU, k = 8, bs = "mpi"), data = dynamic_period)

# Generate predictions for smooth line
cfu_seq_mel <- data.frame(log_CFU = seq(min(dynamic_period$log_CFU),
                                        max(dynamic_period$log_CFU),
                                        length.out = 100))
cfu_seq_mel$smoothed_melanization <- predict(gam_cfu_mel, newdata = cfu_seq_mel)


p4A <- ggplot(dynamic_period, aes(x = log_CFU, y = activity)) +
  geom_jitter(aes(color = Time), size = 2.5, width = 0, height = 0.12, alpha = 0.85) +
  geom_line(data = cfu_seq_act, aes(x = log_CFU, y = smoothed_activity),
            color = "black", linewidth = 0.9, inherit.aes = FALSE) +
  scale_color_viridis_c(option = "magma", direction = -1, end = 0.92,
                        name = "Time (h)") +
  labs(x = expression(log[10](CFU)), y = "Activity") +
  mytheme +
  theme(legend.position = c(0.85, 0.75),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.key.height = unit(0.5, "cm"), 
  )

# Fig 4B: melanization shown in SEVERITY orientation (higher = more melanized = worse).
p4B <- ggplot(dynamic_period, aes(x = log_CFU, y = melanization)) +
  geom_jitter(aes(color = Time), size = 2.5, width = 0, height = 0.12, alpha = 0.85) +
  geom_line(data = cfu_seq_mel, aes(x = log_CFU, y = smoothed_melanization),
            color = "black", linewidth = 0.9, inherit.aes = FALSE) +
  scale_color_viridis_c(option = "magma", direction = -1, end = 0.92,
                        name = "Time (h)") +
  labs(x = expression(log[10](CFU)), y = "Melanization") +
  mytheme +
  theme(legend.position = "none")  # share legend with p4B visually

# ----------------------------------------------------------------------------
# Individual-level binary points (replaces obs_props block)
# ----------------------------------------------------------------------------
metric_offsets <- c("Activity (AT50)"     = -0.5,
                    "Melanization (MT50)" =  0.0,
                    "Survival (LT50)"     =  0.5)

indiv_pts <- dat_cs %>%
  pivot_longer(c(active, unmel, alive),
               names_to = "metric_short", values_to = "y") %>%
  mutate(metric = recode(metric_short,
                         active = "Activity (AT50)",
                         unmel  = "Melanization (MT50)",
                         alive  = "Survival (LT50)"),
         time_plot = time + metric_offsets[metric])

# ----------------------------------------------------------------------------
# Panel C: T50 ordering with individual binary points
# ----------------------------------------------------------------------------

# Observed proportions per timepoint, per metric (for overlay)
obs_props <- dat_cs %>%
  pivot_longer(c(active, unmel, alive),
               names_to = "metric_short", values_to = "y") %>%
  group_by(time, metric_short) %>%
  summarise(p = mean(y), n = n(), .groups = "drop") %>%
  mutate(metric = recode(metric_short,
                         active = "Activity (AT50)",
                         unmel  = "Melanization (MT50)",
                         alive  = "Survival (LT50)"))
obs_props %>% arrange(time, metric) %>% print(n = Inf)
# Panel A
p4C <- ggplot(curves, aes(x = time, y = p, color = metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(data = obs_props,
             aes(x = time, y = p, color = metric),
             alpha = 0.75, size = 4,
             position = position_dodge(width = 0.6)) +
  geom_segment(data = T50_summary,
               aes(x = T50, xend = T50, y = 0, yend = 0.5, color = metric),
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.4) +
  scale_color_manual(values = palette_T50, name = NULL) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(limits = c(0, 36), breaks = seq(0, 36, 6)) +
  labs(x = "Time (h)", y = "Proportion") +
  mytheme +
  theme(legend.position = c(0.28, 0.2),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.text = element_text(size = 11))

# ----------------------------------------------------------------------------
# Compose Figure 4
# ----------------------------------------------------------------------------
figure4 <- (p4A | p4B | p4C) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))+
  theme(legend.title = element_text(size = 12))

ggsave("figures/figure4.pdf",  # >>> Manuscript Figure 4 (activity/melanization + T50)
       plot = figure4, width = 14, height = 4.5, units = "in", dpi = 300)

# ----------------------------------------------------------------------------
# Hysteresis quantification (for the text)
# ----------------------------------------------------------------------------
# At matched CFU, how much does activity drop / melanization rise per hour?
# Use a simple linear model with both predictors so the time slope is
# the hysteresis effect after accounting for instantaneous burden.

m_act_hyst <- lm(activity     ~ log_CFU + Time, data = dynamic_period)
m_mel_hyst <- lm(melanization ~ log_CFU + Time, data = dynamic_period)

cat("\n--- Hysteresis (effect of Time | log_CFU) ---\n")
cat("Activity:    "); print(round(summary(m_act_hyst)$coefficients["Time", ], 4))
cat("Melanization:"); print(round(summary(m_mel_hyst)$coefficients["Time", ], 4))

#==============================================================================
# Overview figure: survival | pathogen growth | composite-health decline
# Reuses p1A (survival, Fig 1) and figure2 (logistic burden growth, Fig 2),
# and adds a panel fitting a logistic DECLINE to composite health over time
# (the mirror of logistic growth). Falls back to a monotone-decreasing smooth
# if the nls fit does not converge.
#==============================================================================
hdat <- burden_tidy_time %>% filter(!is.na(health_combined), !is.na(Time))
t_grid_h <- data.frame(Time = seq(min(hdat$Time, na.rm = TRUE),
                                  max(hdat$Time, na.rm = TRUE), length.out = 200))

# Logistic decline: health falls from an upper plateau U to a lower plateau L,
# with steepness k and midpoint t50 (k > 0 => decreasing in Time).
fit_health_decline <- tryCatch(
  nlsLM(health_combined ~ L + (U - L) / (1 + exp(k * (Time - t50))),
        data  = hdat,
        start = list(U   = max(hdat$health_combined, na.rm = TRUE),
                     L   = min(hdat$health_combined, na.rm = TRUE),
                     k   = 0.2,
                     t50 = median(hdat$Time, na.rm = TRUE)),
        lower = c(U = 0, L = 0, k = 0, t50 = 0),
        control = nls.lm.control(maxiter = 200)),
  error = function(e) NULL)

if (!is.null(fit_health_decline)) {
  t_grid_h$fit <- predict(fit_health_decline, newdata = t_grid_h)
  cat("\n=== Composite health decline (logistic fit) ===\n")
  print(round(coef(fit_health_decline), 3))
} else {
  sc_health <- scam(health_combined ~ s(Time, k = 8, bs = "mpd"), data = hdat)
  t_grid_h$fit <- predict(sc_health, newdata = t_grid_h)
  cat("\nLogistic health-decline fit did not converge; using monotone-decreasing smooth.\n")
}

p_health <- ggplot(hdat, aes(x = Time, y = health_combined)) +
  geom_jitter(aes(fill = status), size = 3, shape = 21, alpha = 0.6,
              color = "black", width = 0.4, height = 0.1) +
  geom_line(data = t_grid_h, aes(x = Time, y = fit),
            color = "black", linewidth = 1.2) +
  scale_fill_manual(name = "Status",
                    values = c("Alive" = "#19798b", "Dead" = "#ee9b43")) +
  labs(x = "Time (h)", y = "Health score") +
  scale_x_continuous(breaks = c(4, 12, 20, 28, 36)) +
  mytheme +
  theme(legend.position = c(0.8, 0.85))


p1A_overview <- p1A + theme(legend.position  = c(0.25, 0.2)) 
figure2_overview <- figure2 + theme(legend.position = c(0.25, 0.8))
figure_overview <- (p1A_overview | figure2_overview | p_health) +
  plot_annotation(tag_levels = "A")

ggsave("figures/figure_overview.pdf",
       plot = figure_overview, width = 15, height = 5, units = "in", dpi = 300)

#===============================================================================
### ANTIBIOTIC TREATMENT ###
#FIGURE 7# 
#===============================================================================

ab_data   <- read.table("data/bacterial_burden_ab.csv", header = T, sep = ",", dec =".")
expdata   <- read.table("data/Galleria_AB_Data_3rd_trial.csv", header = T, sep = ",", dec =".")

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
  )


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

survival_summary <- expdata_filled %>%
  group_by(Sample, Treatment) %>%
  summarise(
    n_alive = sum(Survival == 2, na.rm = TRUE),
    n_total = sum(!is.na(Survival)),   # count rows with valid Survival
    survival_prob = n_alive / n_total
  )

# Summarize mean ± SE of total health by group
health_summary <- expdata_filled %>%
  group_by(Sample, Treatment) %>%
  summarise(
    mean_health = mean(Total_health, na.rm = TRUE),
    se_health   = sd(Total_health, na.rm = TRUE) / sqrt(n())
  )

# tidy plate counts 
burden_tidy_ab <- ab_data %>%
  pivot_longer(cols = rep1:rep5, names_to = "replicates", values_to = "cfu") %>%
  filter(cfu > 1) %>%
  mutate(dilution_factor = (Dilution/Volume_ul)*Buffer_dilution_factor) %>%
  mutate(count = cfu*dilution_factor) %>%
  group_by(Sample, Larvae) %>%
  summarise(cfu = mean(count, na.rm = T)) %>%
  ungroup() %>%
  mutate(Sample = as.character(Sample), Larvae = as.character(Larvae)) %>%
  full_join(expdata_filled %>% mutate(Sample = as.character(Sample), Larvae = as.character(Larvae)), 
            by = c("Sample", "Larvae")) %>%
  mutate(
    
    cfu = ifelse(is.nan(cfu), 0, cfu),  
    status = case_when(
      Survival == 2 ~ "Alive",  
      Survival == 0 ~ "Dead",    
      TRUE ~ NA_character_  
    )
  ) %>%
  mutate(
  ) %>%
  mutate(melanization = 4 - Melanization, 
         log_CFU = log10(cfu+1),
         scaled_health = scale(Total_health), 
         scaled_activity = scale(Activity), 
         scaled_melanization = scale(Melanization), 
         scaled_cfu = scale(log_CFU),
         scaled_immune_morb = scale(Activity + Melanization)
  )

burden_tidy_ab_sum <- burden_tidy_ab %>% 
  filter(!Sample == "G0" & 
           !Sample == "G1") %>% 
  group_by(Treatment) %>%
  summarise(
    mean_cfu = mean(log_CFU, na.rm = TRUE),
    se_cfu   = sd(log_CFU, na.rm = TRUE) / sqrt(n())
  )

# Create a cleaner version with better visual distinction
library(scales)

# Improved color palette emphasizing early vs late
pal_improved <- c(
  "PBS-PBS" = "#d9d9d9",      # lightest grey - vehicle control
  "PBS-CIP" = "#969696",      # medium grey - injection control  
  "PAO1-PBS" = "#252525",     # dark grey/black - no treatment
  "PAO1-00hCIP" = "#1b7837",  # dark green - immediate rescue
  "PAO1-03hCIP" = "#5aae61",  # medium green - early rescue
  "PAO1-06hCIP" = "#fdd49e",  # light orange - partial rescue
  "PAO1-09hCIP" = "#fc8d59",  # medium orange - limited rescue
  "PAO1-12hCIP" = "#d7301f"   # red - too late
)

# Nice labels for all panels
treatment_labels <- c(
  "PBS-PBS" = "Vehicle control",
  "PBS-CIP" = "Injection control",
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
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  scale_y_discrete(labels = treatment_labels) +
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

# Let's check if any error bars are suspiciously small
cat("\nHealth error bar widths:\n")
health_summary_ci %>%
  mutate(bar_width = ci_high - ci_low) %>%
  dplyr::select(Treatment, mean_health, se_health, bar_width) %>%
  print()

p7B <- ggplot(health_summary_ci, aes(y = Treatment, x = mean_health, fill = Treatment)) +
  geom_point(size = 4.5, shape = 21, color = "black", stroke = 0.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                 height = 0.25, linewidth = 0.5) +  
  scale_fill_manual(values = pal_improved) +
  scale_x_continuous(limits = c(0, 9.5), 
                     breaks = seq(0, 9, by = 3),
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
# Include PBS controls with their actual (zero) CFU values
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

# If controls are missing from burden_tidy_ab, add them manually
controls_burden <- expdata_filled %>%
  filter(Treatment %in% c("PBS-PBS", "PBS-CIP")) %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean_cfu = 0,  # Controls have no bacteria
    sd_cfu = 0,
    se_cfu = 0,
    ci_low = 0,
    ci_high = 0,
    .groups = "drop"
  )

# Combine with infected groups if needed
burden_summary_ci <- bind_rows(burden_summary_ci, controls_burden)

p7C <- ggplot(burden_summary_ci, aes(y = Treatment, x = mean_cfu, fill = Treatment)) +
  geom_point(size = 4.5, shape = 21, color = "black", stroke = 0.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                 height = 0.25, linewidth = 0.5) +  # Increased from 0.5
  scale_fill_manual(values = pal_improved) +
  scale_x_continuous(limits = c(-0.2, 4),  # Start slightly below 0 to show controls
                     breaks = seq(0, 4, by = 1),
                     expand = c(0.02, 0)) +
  scale_y_discrete(labels = treatment_labels) +
  labs(x = expression(paste("Pathogen burden (log"[10], " CFU)")), y = NULL) +
  mytheme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank()
  )

#-------------------------------------------------------------------------------
# COMBINE
#-------------------------------------------------------------------------------
figure7 <- (p7A | p7B | p7C) +
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(
    tag_levels = 'A'
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.02, 1.015)
  )

print(figure7)

ggsave("figures/figure7.pdf",  # >>> Manuscript Figure 7 (antibiotic intervention)
       plot = figure7, width = 12, height = 5, units = "in", dpi = 300)


# =============================================================================
# Treatment-timing dose-response (continuous) — capstone for the treatment story
# -----------------------------------------------------------------------------
# Every larva is assessed at the SAME 24 h endpoint, so variation in outcome vs
# minutes-from-injection-to-treatment reflects the manipulated DURATION of
# pathogen exposure before clearance -> a direct, prospective test of the
# cumulative-damage hypothesis (the discrete groups of Fig 7, made continuous
# using the exact per-larva treatment time).
# =============================================================================

# Infected larvae that received ciprofloxacin at a recorded delay.
# (Adjust the filter if you'd rather match the Sample subset G3-G7.)
dose <- burden_tidy_ab %>%
  filter(grepl("PAO1", Treatment), grepl("CIP", Treatment),
         !is.na(min_inj_to_treat), !is.na(Survival)) %>%
  mutate(alive01 = as.integer(Survival == 2))

t_grid_d <- data.frame(min_inj_to_treat = seq(min(dose$min_inj_to_treat, na.rm = TRUE),
                                              max(dose$min_inj_to_treat, na.rm = TRUE),
                                              length.out = 200))

# (1) Survival ~ delay: logistic regression -> treatment window (P = 0.5)
m_surv_dose <- glm(alive01 ~ min_inj_to_treat, data = dose, family = binomial)
t_grid_d$surv <- predict(m_surv_dose, newdata = t_grid_d, type = "response")
bcoef   <- coef(m_surv_dose)
window50 <- as.numeric(-bcoef[1] / bcoef[2])   # delay at which P(survival) = 0.5
cat("\n=== Treatment-timing dose-response (24 h endpoint) ===\n")
cat(sprintf("Survival: P(survival) = 0.5 at %.0f min (%.1f h) after injection\n",
            window50, window50 / 60))
win_se <- tryCatch({ dp <- MASS::dose.p(m_surv_dose, p = 0.5); attr(dp, "SE")[1] },
                   error = function(e) NA_real_)
if (!is.na(win_se)) cat(sprintf("   (approx SE %.0f min)\n", win_se))

# (2) Composite health ~ delay: linear trend (your Total Health figure)
m_health_dose <- lm(Total_health ~ min_inj_to_treat, data = dose)
t_grid_d$health <- predict(m_health_dose, newdata = t_grid_d)
cat(sprintf("Health: slope %.3f per h, p = %.3g\n",
            coef(m_health_dose)[2] * 60, summary(m_health_dose)$coefficients[2, 4]))

# (3) Residual CFU at 24 h ~ delay: later treatment -> more uncleared burden
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

# (5) Sam's key robustness check (comment #230): did the antibiotic actually
# reduce burden? If LATE-treated burden is not below UNTREATED, late failure
# could be pharmacological (inoculum/establishment effect) rather than accrued
# damage. Adjust "PAO1-PBS" if your untreated infected group is named otherwise.
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

pD_surv <- ggplot(dose, aes(min_inj_to_treat, alive01)) +
  geom_jitter(height = 0.04, width = 6, alpha = 0.5, shape = 21, fill = "grey40") +
  geom_line(data = t_grid_d, aes(y = surv), color = "#b80422", linewidth = 1.2) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "P(survival at 24 h)") +
  mytheme

pD_health <- ggplot(dose, aes(min_inj_to_treat, Total_health)) +
  geom_jitter(alpha = 0.5, width = 6, shape = 21, fill = "grey40") +
  geom_line(data = t_grid_d, aes(y = health), color = "#19798b", linewidth = 1.2) +
  labs(x = "Minutes from injection to treatment", y = "Composite health (24 h)") +
  mytheme

pD_cfu <- ggplot(dose, aes(min_inj_to_treat, log_CFU)) +
  geom_jitter(alpha = 0.5, width = 6, shape = 21, fill = "grey40") +
  geom_line(data = t_grid_d, aes(y = logcfu), color = "#ee9b43", linewidth = 1.2) +
  labs(x = NULL, y = bquote(log[10]~CFU~"(24 h)")) +
  mytheme

figure8 <- (pD_surv | pD_health | pD_cfu) +
  plot_annotation(tag_levels = "A")

ggsave("figures/figure8.pdf",  # >>> Manuscript Figure 8 (continuous treatment-timing dose-response)
       plot = figure8, width = 12, height = 4, units = "in", dpi = 300)


# Survival stats
survival_summary_ci %>%
  mutate(
    ci_text = sprintf("%.1f%% [%.1f%%, %.1f%%]", 
                      survival_prob * 100, 
                      ci_low * 100, 
                      ci_high * 100)
  ) %>%
  dplyr::select(Treatment, ci_text) %>%
  print()


health_summary_ci %>%
  mutate(
    ci_text = sprintf("%.2f [%.2f, %.2f]", 
                      mean_health, ci_low, ci_high)
  ) %>%
  dplyr::select(Treatment, ci_text) %>%
  print()

burden_summary_ci %>%
  mutate(
    ci_text = sprintf("%.2f [%.2f, %.2f]", 
                      mean_cfu, ci_low, ci_high)
  ) %>%
  dplyr::select(Treatment, ci_text) %>%
  print()

# Key comparisons
cat("Early (0-3h) vs No treatment:\n")
early_surv <- survival_summary_ci %>% 
  filter(Treatment %in% c("PAO1-00hCIP", "PAO1-03hCIP")) %>%
  summarise(mean_surv = mean(survival_prob))
no_treat_surv <- survival_summary_ci %>% 
  filter(Treatment == "PAO1-PBS") %>%
  pull(survival_prob)
cat(sprintf("  Survival: %.1f%% vs %.1f%% (%.1f percentage point increase)\n",
            early_surv$mean_surv * 100, no_treat_surv * 100,
            (early_surv$mean_surv - no_treat_surv) * 100))

late_surv <- survival_summary_ci %>% 
  filter(Treatment == "PAO1-12hCIP") %>%
  pull(survival_prob)
cat(sprintf("\nLate (12h) vs No treatment:\n"))
cat(sprintf("  Survival: %.1f%% vs %.1f%% (no significant difference)\n",
            late_surv * 100, no_treat_surv * 100))

# Burden comparison
early_burden <- burden_summary_ci %>%
  filter(Treatment %in% c("PAO1-00hCIP", "PAO1-03hCIP")) %>%
  summarise(mean_cfu = mean(mean_cfu))
late_burden <- burden_summary_ci %>%
  filter(Treatment == "PAO1-12hCIP") %>%
  pull(mean_cfu)
no_treat_burden <- burden_summary_ci %>%
  filter(Treatment == "PAO1-PBS") %>%
  pull(mean_cfu)

cat(sprintf("\nBacterial burden:\n"))
cat(sprintf("  Early treatment: log10 %.2f CFU\n", early_burden$mean_cfu))
cat(sprintf("  Late treatment: log10 %.2f CFU\n", late_burden))
cat(sprintf("  No treatment: log10 %.2f CFU\n", no_treat_burden))
cat(sprintf("  Late treatment reduced burden by %.2f log10 units vs untreated\n",
            no_treat_burden - late_burden))
cat(sprintf("  But survival remained poor (%.1f%% vs %.1f%%)\n",
            late_surv * 100, no_treat_surv * 100))

#-------------------------------------------------------------------------------
# STATISTICS FOR INTERVENTION EXPERIMENT
#-------------------------------------------------------------------------------

# Prepare data - exclude PBS controls for infected comparisons
infected_groups <- expdata_filled %>%
  filter(grepl("PAO1", Treatment))

#-------------------------------------------------------------------------------
# SURVIVAL ANALYSIS
#-------------------------------------------------------------------------------

# Create binary survival outcome
survival_data <- expdata_filled %>%
  mutate(
    survived = ifelse(Survival == 2, 1, 0),
    treatment_time = case_when(
      Treatment == "PAO1-00hCIP" ~ 0,
      Treatment == "PAO1-03hCIP" ~ 3,
      Treatment == "PAO1-06hCIP" ~ 6,
      Treatment == "PAO1-09hCIP" ~ 9,
      Treatment == "PAO1-12hCIP" ~ 12,
      Treatment == "PAO1-PBS" ~ NA_real_
    )
  )

timing_data <- survival_data |> 
  dplyr::filter(!is.na(treatment_time))  # infected groups with a defined time

glm_survival_timing <- glm(survived ~ treatment_time,
                           data = timing_data,
                           family = binomial())
summary(glm_survival_timing)


# Chi-square test across all infected groups
survival_table <- table(infected_groups$Treatment, infected_groups$Survival)
chisq_survival <- chisq.test(survival_table)

# Key comparisons for text
early_vs_untreated <- fisher.test(table(
  expdata_filled %>% 
    filter(Treatment %in% c("PAO1-PBS", "PAO1-03hCIP")) %>%
    dplyr::select(Treatment, Survival)
))

late_vs_untreated <- fisher.test(table(
  expdata_filled %>% 
    filter(Treatment %in% c("PAO1-PBS", "PAO1-12hCIP")) %>%
    dplyr::select(Treatment, Survival)
))

cat("\nEarly treatment (3h) vs No treatment:\n")
print(early_vs_untreated)
cat("\nLate treatment (12h) vs No treatment:\n")
print(late_vs_untreated)

#-------------------------------------------------------------------------------
# HEALTH SCORE ANALYSIS
#-------------------------------------------------------------------------------

lm_health_timing <- lm(Total_health ~ treatment_time, data = timing_data)
summary(lm_health_timing)

#-------------------------------------------------------------------------------
# BACTERIAL BURDEN ANALYSIS
#-------------------------------------------------------------------------------

burden_timing <- burden_tidy_ab |>
  dplyr::filter(!is.na(min_inj_to_treat), !is.na(log_CFU))

burden_timing <- burden_timing %>%
  mutate(treatment_hours = min_inj_to_treat / 60)

lm_burden_timing <- lm(log_CFU ~ treatment_hours, data = burden_timing)
summary(lm_burden_timing)

# =============================================================================
# Causal Analysis Figures for Galleria Paper
# =============================================================================

# =============================================================================
# Three-node DAG panels function
# =============================================================================

# Custom function to draw DAG panel
draw_dag_panel <- function(edges, title, subtitle) {
  
  
  # Node positions (matching uploaded image layout)
  nodes <- data.frame(
    name = c("t", "p", "s"),
    x = c(0, 1, 2),
    y = c(0, 1, 0),
    label = c("italic(t)", "italic(p)", "italic(s)")
  )
  
  # Create edge dataframe based on which edges are present
  edge_df <- data.frame()
  
  if ("t_p" %in% edges) {
    edge_df <- rbind(edge_df, data.frame(
      x = 0, y = 0, xend = 1, yend = 1, curve = 0
    ))
  }
  if ("p_s" %in% edges) {
    edge_df <- rbind(edge_df, data.frame(
      x = 1, y = 1, xend = 2, yend = 0, curve = 0
    ))
  }
  if ("t_s" %in% edges) {
    edge_df <- rbind(edge_df, data.frame(
      x = 0, y = 0, xend = 2, yend = 0, curve = 0
    ))
  }
  
  # Box dimensions
  box_w <- 0.35
  box_h <- 0.25
  
  # Create node rectangles
  node_rects <- nodes %>%
    mutate(
      xmin = x - box_w/2,
      xmax = x + box_w/2,
      ymin = y - box_h/2,
      ymax = y + box_h/2
    )
  
  # Adjust arrow endpoints to stop at box edges
  adjust_arrow <- function(x1, y1, x2, y2, box_w, box_h) {
    # Direction vector
    dx <- x2 - x1
    dy <- y2 - y1
    len <- sqrt(dx^2 + dy^2)
    
    # Unit vector
    ux <- dx / len
    uy <- dy / len
    
    # Adjust start point (move away from center of start node)
    # Find intersection with box edge
    if (abs(ux) > 0.01) {
      t_x <- (box_w/2) / abs(ux)
    } else {
      t_x <- Inf
    }
    if (abs(uy) > 0.01) {
      t_y <- (box_h/2) / abs(uy)
    } else {
      t_y <- Inf
    }
    t_start <- min(t_x, t_y)
    
    x1_adj <- x1 + ux * t_start
    y1_adj <- y1 + uy * t_start
    
    # Adjust end point (stop before center of end node)
    x2_adj <- x2 - ux * t_start
    y2_adj <- y2 - uy * t_start
    
    return(c(x1_adj, y1_adj, x2_adj, y2_adj))
  }
  
  # Adjust all edges
  if (nrow(edge_df) > 0) {
    edge_df_adj <- edge_df %>%
      rowwise() %>%
      mutate(
        adj = list(adjust_arrow(x, y, xend, yend, box_w, box_h)),
        x_adj = adj[1],
        y_adj = adj[2],
        xend_adj = adj[3],
        yend_adj = adj[4]
      ) %>%
      ungroup()
  }
  
  # Build plot
  p <- ggplot() +
    # Draw boxes
    geom_rect(data = node_rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", color = "#19798b", linewidth = 0.7) +
    # Draw arrows
    {if (nrow(edge_df) > 0) 
      geom_segment(data = edge_df_adj,
                   aes(x = x_adj, y = y_adj, xend = xend_adj, yend = yend_adj),
                   arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
                   color = "#19798b", linewidth = 0.7)
    } +
    # Draw node labels
    geom_text(data = nodes, aes(x = x, y = y, label = label),
              parse = TRUE, size = 7, color = "black") +
    # Add subtitle (conditional independence statement)
    annotate("text", x = 1, y = -0.5, label = subtitle, 
             size = 7, color = "black", parse = T) +
    # Styling
    coord_fixed(ratio = 1, xlim = c(-0.5, 2.5), ylim = c(-0.7, 1.4)) +
    theme_void() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(p)
}

#===============================================================================
# Four-node DAG panels function
#===============================================================================

draw_dag_panel_4node <- function(edges, title, subtitles, dashed_edges = NULL) {
  
  nodes <- data.frame(
    name = c("t", "p", "h", "s"),
    x = c(0, 1, 1, 2),
    y = c(0.5, 1.2, -0.2, 0.5),
    label = c("italic(t)", "italic(p)", "italic(h)", "italic(s)")
  )
  
  edge_list <- list(
    "t_p" = c(0, 0.5, 1, 1.2),
    "t_h" = c(0, 0.5, 1, -0.2),
    "t_s" = c(0, 0.5, 2, 0.5),
    "p_h" = c(1, 1.2, 1, -0.2),
    "p_s" = c(1, 1.2, 2, 0.5),
    "h_s" = c(1, -0.2, 2, 0.5),
    "h_p" = c(1, -0.2, 1, 1.2)
  )
  
  box_w <- 0.35; box_h <- 0.25
  
  node_rects <- nodes %>%
    mutate(xmin = x - box_w/2, xmax = x + box_w/2,
           ymin = y - box_h/2, ymax = y + box_h/2)
  
  adjust_arrow <- function(x1, y1, x2, y2, box_w, box_h) {
    dx <- x2 - x1; dy <- y2 - y1; len <- sqrt(dx^2 + dy^2)
    if (len == 0) return(c(x1, y1, x2, y2))
    ux <- dx / len; uy <- dy / len
    t_x <- if (abs(ux) > 0.01) (box_w/2)/abs(ux) else Inf
    t_y <- if (abs(uy) > 0.01) (box_h/2)/abs(uy) else Inf
    t_start <- min(t_x, t_y)
    c(x1 + ux*t_start, y1 + uy*t_start, x2 - ux*t_start, y2 - uy*t_start)
  }
  
  # Helper: build adjusted edge df from a vector of edge codes
  build_edges <- function(edge_codes) {
    df <- data.frame()
    for (e in edge_codes) {
      if (e %in% names(edge_list)) {
        co <- edge_list[[e]]
        df <- rbind(df, data.frame(x = co[1], y = co[2], xend = co[3], yend = co[4]))
      }
    }
    if (nrow(df) > 0) {
      df <- df %>% rowwise() %>%
        mutate(adj = list(adjust_arrow(x, y, xend, yend, box_w, box_h)),
               x_adj = adj[1], y_adj = adj[2], xend_adj = adj[3], yend_adj = adj[4]) %>%
        ungroup()
    }
    df
  }
  
  edge_df_adj   <- build_edges(edges)
  dashed_df_adj <- build_edges(dashed_edges)
  
  p <- ggplot() +
    geom_rect(data = node_rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", color = "#19798b", linewidth = 0.5) +
    # solid edges
    {if (nrow(edge_df_adj) > 0)
      geom_segment(data = edge_df_adj,
                   aes(x = x_adj, y = y_adj, xend = xend_adj, yend = yend_adj),
                   arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
                   color = "#19798b", linewidth = 0.5)} +
    # grey dashed "also-rejected variant" edges
    {if (nrow(dashed_df_adj) > 0)
      geom_segment(data = dashed_df_adj,
                   aes(x = x_adj, y = y_adj, xend = xend_adj, yend = yend_adj),
                   arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
                   color = "grey60", linetype = "dashed", linewidth = 0.5)} +
    geom_text(data = nodes, aes(x = x, y = y, label = label),
              parse = TRUE, size = 5, color = "black") +
    coord_fixed(ratio = 1, xlim = c(-0.5, 2.5), ylim = c(-1.1, 1.7)) +
    theme_void() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "italic"))
  
  for (i in seq_along(subtitles)) {
    p <- p + annotate("text", x = 1, y = -0.65 - (i-1)*0.25,
                      label = subtitles[i], size = 4, color = "black", parse = TRUE)
  }
  p
}


# Show correlation between Activity and Melanization.
# Correlate them in the SAME (health) orientation they enter the composite:
# activity (higher = healthier) and (4 - melanization) = raw melanization
# (higher = healthier). Using the stored severity `melanization` here would only
# flip the SIGN of r (magnitude and p-value are identical) and report a negative
# r that contradicts the "both measure health, so sum them" justification.
melan_health <- 4 - burden_tidy_time$melanization      # raw melanization, 4 = healthy
cor_AM <- cor(burden_tidy_time$activity, melan_health,
              use = "complete.obs")
cat("Correlation between Activity and Melanization (both health-oriented):", round(cor_AM, 3), "\n")

# Test correlation significance
cor_test <- cor.test(burden_tidy_time$activity, melan_health)
cat("Correlation p-value:", cor_test$p.value, "\n\n") 

# Create combined health metric

burden_tidy_time <- burden_tidy_time %>%
  mutate(
    # Combined health = activity + (4 - melanization).
    # `melanization` is the severity form (4 - raw), so (4 - melanization)
    # recovers the raw health-oriented value (4 = healthy). Higher = healthier.
    health_combined = activity + (4 - melanization),
    scaled_health_combined = scale(health_combined)[,1]
  )

#==============================================================================
# Causal model formalization (dagitty)
#------------------------------------------------------------------------------
# Encodes the candidate DAGs underlying Figures 3 and 5 and prints the testable
# conditional independencies derived by d-separation. These are exactly the
# implications tested by regression downstream:
#   Fig 3: m_cond            (t -> p -> s implies t _||_ s | p ; etc.)
#   Fig 5: m_D, m_E, m_F, m_F_supp
# Reference: Textor et al. 2016 (dagitty). dagitty is loaded above.
#==============================================================================

# --- Three-node candidate models (Figure 3: t, p, s) ---
dagi3_mediation <- dagitty("dag { t -> p -> s }")            # 3A: pure pathogen mediation
dagi3_no_med    <- dagitty("dag { t -> p ; t -> s }")        # 3B: no pathogen mediation
dagi3_multipath <- dagitty("dag { t -> p -> s ; t -> s }")   # 3C: multi-path (supported)

cat("\n=== Implied conditional independencies: 3-node models (Figure 3) ===\n")
cat("\n3A  t -> p -> s  (pure pathogen mediation):\n")
print(impliedConditionalIndependencies(dagi3_mediation))   # expect: t _||_ s | p
cat("\n3B  t -> p ; t -> s  (no pathogen mediation):\n")
print(impliedConditionalIndependencies(dagi3_no_med))      # expect: p _||_ s | t
cat("\n3C  t -> p -> s ; t -> s  (multi-path, supported):\n")
print(impliedConditionalIndependencies(dagi3_multipath))   # saturated: no implications

# --- Four-node candidate models (Figure 5: t, p, h, s) ---
dagi4_pmed_health <- dagitty("dag { t -> p -> h -> s }")                  # 5A: pathogen-mediated health
dagi4_collapse    <- dagitty("dag { t -> h -> p -> s }")                  # 5B: immune collapse
dagi4_bottleneck  <- dagitty("dag { t -> p ; t -> h ; p -> h ; h -> s }") # 5C: health bottleneck (supported)

cat("\n=== Implied conditional independencies: 4-node models (Figure 5) ===\n")
cat("\n5A  t -> p -> h -> s  (pathogen-mediated health):\n")
print(impliedConditionalIndependencies(dagi4_pmed_health))  # includes t _||_ h | p  (Fig 5D)
cat("\n5B  t -> h -> p -> s  (immune collapse):\n")
print(impliedConditionalIndependencies(dagi4_collapse))     # includes h _||_ s | p  (Fig 5E)
cat("\n5C  t -> p ; t -> h ; p -> h ; h -> s  (supported):\n")
print(impliedConditionalIndependencies(dagi4_bottleneck))   # s _||_ t | h and s _||_ p | h  (Fig 5F)

#==============================================================================
# FIGURE 3 (NEW): 3-node DAGs + diagnostic  plots
# Panels A-C: DAGs for alternative causal models (t,p,S)
# Panels D-E: Diagnostic plots falsifying models A and B
#==============================================================================

# --- DAG panels ---

dag_A <- draw_dag_panel(
  edges = c("t_p", "p_s"),
  title = "",
  subtitle = "italic(t)~symbol('\\136')~italic(s)~'|'~italic(p)"
)

dag_B <- draw_dag_panel(
  edges = c("t_p", "t_s"),
  title = "",
  subtitle = "italic(p)~symbol('\\136')~italic(s)~'|'~italic(t)"
)

dag_C <- draw_dag_panel(
  edges = c("t_p", "p_s", "t_s"),
  title = "",
  subtitle = ""
)

# This is the final figure 3

# ----------------------------------------------------------------------------
# Data prep
# ----------------------------------------------------------------------------
dat_fig3 <- burden_tidy %>%
  filter(Time < 37, !is.na(status), !is.na(log_CFU)) %>%
  mutate(
    alive    = as.integer(status == "Alive"),
    cfu_bin  = cut(log_CFU,
                   breaks = quantile(log_CFU, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                   labels = c("Low", "Medium", "High"),
                   include.lowest = TRUE),
    time_bin = cut(Time,
                   breaks = c(-0.01, 12, 24, 48),
                   labels = c("0–12h", "12–24h", "24–48h"))
  )


# ----------------------------------------------------------------------------
# Single additive logistic GLM — matches the test you report (β_Time, β_logCFU)
# ----------------------------------------------------------------------------
m_cond <- glm(alive ~ Time + log_CFU, data = dat_fig3, family = binomial)
summary(m_cond)  # β_Time and β_logCFU are your reported test statistics

# Illustrative levels of the conditioning variable
cfu_levels  <- quantile(dat_fig3$log_CFU, c(0.15, 0.5, 0.85), na.rm = TRUE)
time_levels <- c(6, 18, 30)
names(cfu_levels)  <- c("Low", "Medium", "High")
names(time_levels) <- c("0–12h", "12–24h", "24–48h")


# ─── Null prediction for Panel D: what we'd see if DAG A (t ⊥ S | p) held ───
# Under DAG A, survival depends only on log_CFU. Fit that null model.
m_null_D <- glm(survival ~ log_CFU, data = burden_tidy_time, family = binomial)

null_D <- purrr::map_dfr(names(cfu_levels), function(lbl) {
  nd <- data.frame(Time = c(0, 36), log_CFU = cfu_levels[[lbl]])
  nd$fit <- plogis(predict(m_null_D, newdata = nd))
  nd$cfu_bin <- factor(lbl, levels = c("Low", "Medium", "High"))
  nd
})


# ─── Null prediction for Panel E: what we'd see if DAG B (p ⊥ S | t) held ───
m_null_E <- glm(survival ~ Time, data = burden_tidy_time, family = binomial)

null_E <- purrr::map_dfr(names(time_levels), function(lbl) {
  nd <- data.frame(log_CFU = range(burden_tidy_time$log_CFU, na.rm = TRUE),
                   Time = time_levels[[lbl]])
  nd$fit <- plogis(predict(m_null_E, newdata = nd))
  nd$bin <- factor(lbl, levels = names(time_levels))   # <-- self-consistent
  nd
})


# ----------------------------------------------------------------------------
# Panel D: P(alive) vs Time at three illustrative CFU levels
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
              width = 0.3, height = 0.035, alpha = 0.5, size = 1.7) +
  #geom_ribbon(data = pred_D,
  #            aes(x = Time, ymin = lwr, ymax = upr, fill = bin),
  #            alpha = 0.12) +
  geom_line(data = null_D, aes(x = Time, y = fit, color = cfu_bin),
            linetype = "dashed", linewidth = 0.7, alpha = 0.6)+
  geom_line(data = pred_D,
            aes(x = Time, y = fit, color = bin),
            linewidth = 1.1) +
  scale_color_manual(values = bin_palette,
                     name = expression(log[10](CFU))) +
  scale_fill_manual(values = bin_palette, guide = "none") +
  scale_y_continuous(breaks = c(0, 1), labels = c("Dead", "Alive"),
                     limits = c(-0.12, 1.12)) +
  scale_x_continuous(limits = c(0, 36), breaks = seq(0, 36, 6)) +
  labs(x = "Time (h)", y = "Survival") +
  mytheme +
  theme(legend.position = c(0.25, 0.25),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ----------------------------------------------------------------------------
# Panel E: P(alive) vs log_CFU at three illustrative time levels
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
  nd$bin <- factor(lbl, levels = c("0–12h", "12–24h", "24–48h"))
  nd
})

time_palette <- c("0–12h" = "#19798b", "12–24h" = "#ee9b43", "24–48h" = "#b80422")

p3E <- ggplot() +
  geom_jitter(data = dat_fig3,
              aes(x = log_CFU, y = alive, color = time_bin),
              width = 0.08, height = 0.035, alpha = 0.5, size = 1.7) +
  #geom_ribbon(data = pred_E,
  #            aes(x = log_CFU, ymin = lwr, ymax = upr, fill = bin),
  #            alpha = 0.12) +
  geom_line(data = null_E, aes(x = log_CFU, y = fit, color = bin),
            linetype = "dashed", linewidth = 0.7, alpha = 0.6)+
  geom_line(data = pred_E,
            aes(x = log_CFU, y = fit, color = bin),
            linewidth = 1.1) +
  scale_color_manual(values = time_palette, name = "Time") +
  scale_fill_manual(values = time_palette, guide = "none") +
  scale_y_continuous(breaks = c(0, 1), labels = c("Dead", "Alive"),
                     limits = c(-0.12, 1.12)) +
  labs(x = expression(log[10](CFU)), y = "Survival") +
  mytheme+
  theme(legend.position = c(0.25, 0.35),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Shared display limits for the 3-node DAG panels. These match the coord_fixed()
# limits set inside draw_dag_panel() so the conditional-independence subtitle
# (placed at y = -0.5) is not clipped when the panels are composed by patchwork.
dag_xlim <- c(-0.5, 2.5)
dag_ylim <- c(-0.7, 1.4)

dag_A <- dag_A + coord_cartesian(xlim = dag_xlim, ylim = dag_ylim, clip = "off")
dag_B <- dag_B + coord_cartesian(xlim = dag_xlim, ylim = dag_ylim, clip = "off")
dag_C <- dag_C + coord_cartesian(xlim = dag_xlim, ylim = dag_ylim, clip = "off")


figure3 <- (dag_A | dag_B | dag_C) / (p3D | p3E) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("figures/figure3.pdf",  # >>> Manuscript Figure 3 (3-node DAGs + diagnostics)
       plot = figure3, width = 12, height = 7, units = "in", dpi = 300)


# ============================================================================
# FIGURE 5: Health-mediated causal analysis
# A,B,C: DAGs (h mediates p; p mediates h; h as unique bottleneck)
# D: rejects A via t ⊥ h | p
# E: rejects B via h ⊥ S | p
# F: supports C via S ⊥ t | h (and report S ⊥ p | h in caption)
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Health composite. `melanization` is the severity form (4 - raw), so
#    (4 - melanization) recovers raw (4 = healthy); higher h = healthier.
# ----------------------------------------------------------------------------
dat_fig5 <- burden_tidy_time %>%
  filter(!is.na(activity), !is.na(melanization),
         !is.na(survival), !is.na(log_CFU)) %>%
  mutate(h = activity + (4 - melanization))   # higher h = healthier


dag5_A <- draw_dag_panel_4node(
  edges = c("t_p", "p_h", "h_s"),
  dashed_edges = "p_s",           # p → s also-rejected variant
  title = "", subtitles = "italic(t)~symbol('\\136')~italic(h)~'|'~italic(p)"
)

dag5_B <- draw_dag_panel_4node(
  edges = c("t_h", "h_p", "p_s"),
  dashed_edges = "t_p",           # t → p also-rejected variant
  title = "", subtitles = "italic(h)~symbol('\\136')~italic(S)~'|'~italic(p)"
)

dag5_C <- draw_dag_panel_4node(
  edges = c("t_p", "t_h", "p_h", "h_s"),
  dashed_edges = NULL,            # survivor, no variant
  title = "",
  subtitles = c("italic(s)~symbol('\\136')~italic(t)~'|'~italic(h)",
                "italic(s)~symbol('\\136')~italic(p)~'|'~italic(h)")
)

# ----------------------------------------------------------------------------
# Shared conditioning levels & palettes
# ----------------------------------------------------------------------------
cfu_levels <- quantile(dat_fig5$log_CFU, c(0.15, 0.50, 0.85), na.rm = TRUE)
h_levels   <- quantile(dat_fig5$h,       c(0.15, 0.50, 0.85), na.rm = TRUE)
names(cfu_levels) <- c("Low", "Medium", "High")
names(h_levels)   <- c("Low h", "Medium h", "High h")

pal_cfu <- c("Low" = "#19798b", "Medium" = "#ee9b43", "High" = "#b80422")
pal_h   <- c("Low h" = "#b80422", "Medium h" = "#ee9b43", "High h" = "#19798b") # inverted: high h = good = teal

dat_binned <- dat_fig5 %>%
  mutate(cfu_bin = cut(log_CFU,
                       breaks = quantile(log_CFU, c(0,1/3,2/3,1), na.rm = TRUE),
                       labels = names(cfu_levels), include.lowest = TRUE),
         h_bin = cut(h,
                     breaks = c(-0.01, 2.5, 5.5, 7.01),
                     labels = c("Low h", "Medium h", "High h")))

# ----------------------------------------------------------------------------
# binned data (with pre-computed jitter)
# ----------------------------------------------------------------------------
dat_binned <- dat_fig5 %>%
  mutate(
    cfu_bin = cut(log_CFU,
                  breaks = quantile(log_CFU, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                  labels = c("Low", "Medium", "High"), include.lowest = TRUE),
    h_bin   = cut(h,
                  breaks = c(-0.01, 2.5, 5.5, 7.01),
                  labels = c("Low h", "Medium h", "High h")),
    survival_jitter = survival + runif(n(), -0.04, 0.04)   # pre-computed
  )

cfu_levels <- quantile(dat_fig5$log_CFU, c(0.15, 0.50, 0.85), na.rm = TRUE)
h_levels   <- c("Low h" = 1, "Medium h" = 4, "High h" = 7)
names(cfu_levels) <- c("Low", "Medium", "High")

pal_cfu <- c("Low" = "#19798b", "Medium" = "#ee9b43", "High" = "#b80422")
pal_h   <- c("Low h" = "#b80422", "Medium h" = "#ee9b43", "High h" = "#19798b")


# ----------------------------------------------------------------------------
# 3. Panel D — rejects A: h vs Time at 3 CFU levels
# ----------------------------------------------------------------------------
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
  scale_color_manual(values = pal_cfu, name = expression(log[10](CFU))) +
  scale_x_continuous(limits = c(0, 36), breaks = seq(0, 36, 6)) +
  labs(x = "Time (h)", y = "Health (h)") + mytheme

# ----------------------------------------------------------------------------
# 4. Panel E — rejects B: S vs h at 3 CFU levels
# ----------------------------------------------------------------------------
m_E      <- glm(survival ~ h + log_CFU, data = dat_fig5, family = binomial)
m_E_null <- glm(survival ~ log_CFU,     data = dat_fig5, family = binomial)

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
  scale_color_manual(values = pal_cfu, name = expression(log[10](CFU))) +
  scale_y_continuous(breaks = c(0,1), labels = c("Dead","Alive"), limits = c(-0.12, 1.12)) +
  labs(x = "Health (h)", y = "Survival") + mytheme

# ----------------------------------------------------------------------------
# 5. Panel F — supports C: S vs Time at 3 h levels
# ----------------------------------------------------------------------------
m_F      <- glm(survival ~ Time + h, data = dat_fig5, family = binomial)
m_F_null <- glm(survival ~ h,        data = dat_fig5, family = binomial)
m_F_supp <- glm(survival ~ log_CFU + h, data = dat_fig5, family = binomial)

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
  scale_color_manual(values = pal_h, name = "Health (h)") +
  scale_y_continuous(breaks = c(0,1), labels = c("Dead","Alive"), limits = c(-0.12, 1.12)) +
  scale_x_continuous(limits = c(0, 36), breaks = seq(0, 36, 6)) +
  labs(x = "Time (h)", y = "Survival") + mytheme +
  theme(legend.position = c(0.25, 0.35))

# ----------------------------------------------------------------------------
# 6. Print test stats for caption
# ----------------------------------------------------------------------------
cat("\n=== Figure 5 test statistics ===\n")
cat("D — t ⊥ h | p (DAG A):\n");   print(summary(m_D)$coefficients["Time", ])
cat("\nE — h ⊥ S | p (DAG B):\n"); print(summary(m_E)$coefficients["h", ])
cat("\nF — S ⊥ t | h (DAG C):\n"); print(summary(m_F)$coefficients["Time", ])
cat("F — S ⊥ p | h (DAG C):\n");   print(summary(m_F_supp)$coefficients["log_CFU", ])

# ----------------------------------------------------------------------------
# 7. Assemble
# ----------------------------------------------------------------------------
design <- "
ABC
DEF
"
figure5 <- dag5_A + dag5_B + dag5_C + p5D + p5E + p5F +
  plot_layout(design = design, heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

figure5

ggsave("figures/figure5.pdf",  # >>> Manuscript Figure 5 (4-node health DAGs + diagnostics)
       plot = figure5, width = 14.5, height = 6.5, dpi = 300)

# -----------------------------------------------------------------------------
# Sensitivity analysis: health composite = melanization only (h_mel = 4 - melan.)
# Reproduces the values reported in the text for the melanization-only health
# index (s ⊥ t | h_mel ; s ⊥ p | h_mel). Self-contained — uses *_mel object
# names so it does NOT overwrite the main Figure 5 objects (dat_fig5, m_F, ...).
# -----------------------------------------------------------------------------
dat_fig5_mel <- burden_tidy_time %>%
  filter(!is.na(melanization), !is.na(survival), !is.na(log_CFU)) %>%
  mutate(h_mel = (4 - melanization))   # (4 - melanization) = raw (4 = healthy); higher = healthier

# DAG C support tests with melanization-only health
m_F_mel      <- glm(survival ~ Time    + h_mel, data = dat_fig5_mel, family = binomial)
m_F_supp_mel <- glm(survival ~ log_CFU + h_mel, data = dat_fig5_mel, family = binomial)

cat("\n=== Figure 5 sensitivity: melanization-only health ===\n")
cat("s \u22A5 t | h_mel  (Time | h_mel):\n");    print(summary(m_F_mel)$coefficients["Time", ])
cat("s \u22A5 p | h_mel  (log_CFU | h_mel):\n"); print(summary(m_F_supp_mel)$coefficients["log_CFU", ])

# =============================================================================
# Updated SEM with Health
# =============================================================================

# Prepare time variables
burden_tidy_time$time_linear  <- burden_tidy_time$scaled_time
burden_tidy_time$time_squared <- burden_tidy_time$scaled_time^2

# SEM Model 5: t → p → h → S with t → h (supported model)
model_sem <- '
  # Pathogen growth
  scaled_cfu ~ t1 * time_linear + t2 * time_squared
  
  # Health depends on both pathogen AND time (hysteresis)
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

summary(fit_sem)
print(fitMeasures(fit_sem, c("dic", "ppp")))

# Check convergence
print(blavInspect(fit_sem, "psrf"))

# =============================================================================
# Compare with alternative models
# =============================================================================

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

print(fitMeasures(fit_sem, c("dic", "ppp")))
print(fitMeasures(fit_no_h, c("dic", "ppp")))

dic_diff <- fitMeasures(fit_no_h, "dic") - fitMeasures(fit_sem, "dic")
cat("\nΔDIC (no_h - with_h):", round(dic_diff, 1), "\n")

# =============================================================================
# Extract and Plot Effects
# =============================================================================

sem_summary <- summary(fit_sem)
effects_df_raw <- as.data.frame(sem_summary)

semPaths(fit_sem,
         what = "est",
         layout = "tree2",
         edge.label.cex = 1.2,
         edge.width = 2,
         sizeMan = 10,
         edge.color = "grey20",
         fade = T,
         curvePivot = TRUE,
         label.prop = 1,
         intercepts = FALSE,
         nCharNodes = 0,
         residuals = TRUE,
         nodeLabels = c("Pathogen\nload", "Host\nhealth", "Survival", "Time", "time^2"))

# Convert to dataframe - row names become first column with empty name

sem_summary <- summary(fit_sem)
effects_df_raw <- as.data.frame(sem_summary)

# Extract effects from new SEM output
effects_df <- data.frame(
  path = c(
    "italic(t) %->% italic(p)~(t[1])",
    "italic(t)^2 %->% italic(p)~(t[2])",
    "italic(p) %->% italic(h)~(a)",
    "italic(t) %->% italic(h)~(h[1])",
    "italic(h) %->% italic(s)~(b)",
    "Indirect:~italic(p) %->% italic(h) %->% italic(s)",
    "Indirect:~italic(t) %->% italic(h) %->% italic(s)"
  ),
  estimate = c(
    effects_df_raw["X", "Estimate"],
    effects_df_raw["X.1", "Estimate"],
    effects_df_raw["X.2", "Estimate"],
    effects_df_raw["X.3", "Estimate"],
    effects_df_raw["X.4", "Estimate"],
    effects_df_raw["X.11", "Estimate"],
    effects_df_raw["X.12", "Estimate"]
  ),
  lower = c(
    effects_df_raw["X", "pi.lower"],
    effects_df_raw["X.1", "pi.lower"],
    effects_df_raw["X.2", "pi.lower"],
    effects_df_raw["X.3", "pi.lower"],
    effects_df_raw["X.4", "pi.lower"],
    effects_df_raw["X.11", "pi.lower"],
    effects_df_raw["X.12", "pi.lower"]
  ),
  upper = c(
    effects_df_raw["X", "pi.upper"],
    effects_df_raw["X.1", "pi.upper"],
    effects_df_raw["X.2", "pi.upper"],
    effects_df_raw["X.3", "pi.upper"],
    effects_df_raw["X.4", "pi.upper"],
    effects_df_raw["X.11", "pi.upper"],
    effects_df_raw["X.12", "pi.upper"]
  )
) %>%
  mutate(path = factor(path, levels = rev(path)))

# Convert to numeric
effects_df$estimate <- as.numeric(effects_df$estimate)
effects_df$lower    <- as.numeric(effects_df$lower)
effects_df$upper    <- as.numeric(effects_df$upper)

print(effects_df)

#==============================================================================
# FIGURE 5 (UPDATED): SEM results — supported DAG + effect sizes
# Panel A: Supported 4-node DAG (Model 5: t→p→h→S, t→h)
# Panel B: Effect sizes plot
# semPaths diagram can be generated separately or replaced by the DAG
#==============================================================================

# Panel A: The winning DAG (already exists as dag_5)
dag_sem <- draw_dag_panel_4node(
  edges = c("t_p", "t_h", "p_h", "h_s"),
  title = "Supported model",
  subtitles = character(0)
) +
  theme(
    plot.background  = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 12, face = "italic")
  )

# Panel B: Effect sizes (already exists as effect_plot)
# Recreate to ensure clean state:

effects_df$estimate <- as.numeric(effects_df$estimate)
effects_df$lower    <- as.numeric(effects_df$lower)
effects_df$upper    <- as.numeric(effects_df$upper)

p5B <- ggplot(effects_df, aes(y = path, x = estimate)) +
  geom_col(fill = "grey70", width = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, linewidth = 0.8) +
  labs(y = "", x = "Effect Size") +
  mytheme +
  theme(
    axis.text.y = element_text(hjust = 1, size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )+
  scale_y_discrete(labels = function(x) parse(text = x))

# Compose: bar plot with DAG floating in top-right
#dag_sem <- dag_sem + coord_cartesian(ylim = c(-0.4, 1.6), clip = "off")
figure6 <- p5B +
  inset_element(dag_sem,
                left = 0.62, bottom = 0.02,    # bottom-right corner
                right = 1.00, top = 0.35,
                align_to = "panel",
                on_top = TRUE)

ggsave("figures/figure6.pdf",  # >>> Manuscript Figure 6 (SEM path effects)
       plot = figure6, width = 7.7, height = 7, units = "in", dpi = 300)


#=============================================
# Cumulative burden - For Damage hypothesis
#=============================================

#--------------------------------------------------------------------
# Cumulative burden over time - EXCLUDED from main text
#--------------------------------------------------------------------

time_seq_cum <- seq(0, 36, length.out = 200)
cum_burden_curve <- (K / r) * log1p((p0 / (K - p0)) * expm1(r * time_seq_cum))

df_cum_curve <- data.frame(
  Time = time_seq_cum,
  cum_burden = cum_burden_curve,
  log_cum = log10(cum_burden_curve + 1)
)

inset_cum <- ggplot(df_cum_curve, aes(x = Time, y = log_cum)) +
  geom_line(color = "#b80422", linewidth = 1) +
  labs(x = "Time (h)", 
       y = bquote(log[10](Sigma * p))) +
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

ggsave("figures/cumulative_burden.pdf",
       plot = inset_cum, width = 5, height = 5, units = "in", dpi = 300)


# First, calculate cumulative burden
# Using your fitted logistic parameters from logistic_logfit_full
params_logistic <- coef(logistic_logfit_full)
K <- params_logistic["K"]
p0 <- params_logistic["p0"]
r <- params_logistic["r"]

# For each larva, calculate cumulative burden from t=0 to t=Time
# Integral of logistic: ∫[K/(1 + ((K-p0)/p0)*e^(-rt))]dt
burden_tidy_time <- burden_tidy_time %>%
  mutate(
    # Analytical solution for integral of logistic curve
    #This was numerically unstable, so I rewrote it using log1p and expm1 to improve stability:
    #cum_burden = K * Time - (K/r) * log(1 + ((K-p0)/p0) * exp(-r * Time)) +
    # (K/r) * log(1 + ((K-p0)/p0)),  # subtract initial value
    cum_burden = (K/r) * log1p((p0 / (K - p0)) * expm1(r * Time)),
    scaled_cum_burden = scale(cum_burden)
  )

# which one predicts health better?
# Test models
m_instant <- lm(total_score ~ scaled_cfu, data = burden_tidy_time)
m_cumulative <- lm(total_score ~ scaled_cum_burden, data = burden_tidy_time)
m_both <- lm(total_score ~ scaled_cfu + scaled_cum_burden, data = burden_tidy_time)

# Compare AIC
AIC(m_instant, m_cumulative, m_both)


# Compare instantaneous vs cumulative as predictors of health
p_comparison <- burden_tidy_time %>%
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
  labs(x = "Standardized burden", y = "Standardized health") +
  mytheme


#===============================================================================
# SUPPLEMENTARY ANALYSIS: Cumulative Burden Estimation
# Why cumulative measures are problematic in cross-sectional designs
#===============================================================================

#-------------------------------------------------------------------------------
# METHOD 1: Integral of fitted logistic (current approach)
#-------------------------------------------------------------------------------

# Parameters from fitted logistic model
params_logistic <- coef(logistic_logfit_full)
K  <- params_logistic["K"]
p0 <- params_logistic["p0"]
r  <- params_logistic["r"]

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
# For each timepoint, sum area under curve from t=0 to t=T
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

# The key insight: ALL cumulative measures are near-perfectly correlated with time
cat("\nCorrelation with Time:\n")
cat("  Integral method:    r =", round(cor_matrix["Time", "cum_burden_integral"], 4), "\n")
cat("  Trapezoidal method: r =", round(cor_matrix["Time", "cum_burden_trapezoid"], 4), "\n")
cat("  Discrete method:    r =", round(cor_matrix["Time", "cum_burden_discrete"], 4), "\n")

#-------------------------------------------------------------------------------
# FIGURE S2: Comparison of cumulative burden estimation methods
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
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(x = "Time since infection (h)", 
       y = expression(log[10](Cumulative~burden)),
       color = "Method", linetype = "Method") +
  mytheme +
  theme(legend.position = c(0.3, 0.8))

# Panel B: Integral vs Trapezoidal (method comparison)
method_compare <- burden_tidy_time %>%
 dplyr::select(Time, cum_burden_integral, cum_burden_trapezoid) %>%
  distinct()

pB_supp <- ggplot(method_compare, aes(x = log10(cum_burden_integral + 1), 
                                      y = log10(cum_burden_trapezoid + 1))) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "#E69F00") +
  labs(x = expression(log[10](Integral~of~fitted~model)),
       y = expression(log[10](Trapezoidal~sum~of~data))) +
  annotate("text", x = 6, y = 8, 
           label = paste0("r = ", round(cor(method_compare$cum_burden_integral, 
                                            method_compare$cum_burden_trapezoid), 3)),
           size = 5) +
  mytheme

# Panel C: THE PROBLEM - Cumulative burden vs Time (near-perfect correlation)
pC_supp <- ggplot(method_compare, aes(x = Time, y = log10(cum_burden_integral + 1))) +
  geom_point(size = 3, color = "#E69F00") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Time since infection (h)",
       y = expression(log[10](Cumulative~burden))) +
  annotate("text", x = 10, y = 8, 
           label = paste0("r = ", round(cor(method_compare$Time, 
                                            method_compare$cum_burden_integral), 3)),
           size = 5) +
  annotate("text", x = 10, y = 7, 
           label = "Near-deterministic\nrelationship", size = 4, color = "darkred") +
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
summary(m_resid)

pD_supp <- ggplot(burden_tidy_time, aes(x = cum_burden_resid, y = scaled_health)) +
  geom_point(aes(fill = status), shape = 21, size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_fill_manual(values = c("#19798b", "#ee9b43")) +
  labs(x = "Residual cumulative burden\n(after removing time effect)",
       y = "Standardized health",
       fill = "Status") +
  annotate("text", x = 0, y = 2, 
           label = paste0("β = ", round(coef(m_resid)[2], 3), 
                          ", p = ", round(summary(m_resid)$coefficients[2,4], 3)),
           size = 4) +
  mytheme +
  theme(legend.position = c(0.85, 0.85))

# Combine
figure_cumulative_supp <- (pA_supp | pB_supp) / (pC_supp | pD_supp) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

figure_cumulative_supp

ggsave("figures/figureS2.pdf",  # >>> Supplementary Figure S2 (cumulative-burden methods) 
       plot = figure_cumulative_supp, width = 10, height = 8, dpi = 300)

#-------------------------------------------------------------------------------
# STATISTICAL TEST: Does cumulative burden add information beyond time?
#-------------------------------------------------------------------------------

cat("\n=== Model comparison: Does cumulative burden add info beyond time? ===\n\n")

# Predicting health
m_time_only <- lm(scaled_health ~ Time, data = burden_tidy_time)
m_cum_only  <- lm(scaled_health ~ log10(cum_burden_integral + 1), data = burden_tidy_time)
m_both      <- lm(scaled_health ~ Time + log10(cum_burden_integral + 1), data = burden_tidy_time)

cat("AIC comparison (predicting health):\n")
print(AIC(m_time_only, m_cum_only, m_both))

cat("\nModel with both predictors:\n")
print(summary(m_both))

# Key result: In the combined model, one predictor becomes non-significant
# because they're colinear

# Variance Inflation Factor
library(car)
cat("\nVariance Inflation Factor (VIF) in combined model:\n")
print(vif(m_both))

#-------------------------------------------------------------------------------
# WHY THIS MATTERS FOR SEM
#-------------------------------------------------------------------------------

cat("\n=== Implications for Structural Equation Modeling ===\n\n")

cat("The near-perfect correlation between cumulative burden and time creates problems:\n\n")

cat("1. COLLINEARITY: When both are in a model, estimates become unstable\n")
cat("   VIF > 10 indicates severe collinearity\n\n")

cat("2. IDENTIFIABILITY: SEM cannot distinguish paths like:\n")
cat("   t → Σp → h → S  vs.  t → h → S (with Σp as byproduct)\n")
cat("   because Σp ≈ f(t) in cross-sectional data\n\n")

cat("3. OVERFITTING: Models with Σp achieve low DIC by capturing\n")
cat("   time-correlated noise, leading to poor predictive performance (PPP < 0.01)\n\n")

cat("4. SOLUTION: The antibiotic intervention experiment breaks this collinearity\n")
cat("   by manipulating cumulative exposure while holding assessment time constant.\n")
cat("   Early vs late treatment creates larvae with DIFFERENT cumulative exposures\
")
cat("   assessed at the SAME time (24h), allowing causal identification.\n")

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
  Cons = c("Model-dependent, near-deterministic with t",
           "Sensitive to sampling density",
           "Crude, ignores time intervals")
)

print(summary_table)
