# Load required libraries
library(minpack.lm)  # For non-linear regression
library(ggplot2)   
library(dplyr)
library(tidyr)
library(boot)
library(MASS)
library(mgcv)
library(cowplot)

library(survival)

# SEM
library(lavaan)
library(blavaan) #bayesian lavaan
library(lavaanPlot)
library(semPlot)

#DAG
library(dagitty)
library(ggdag)

library(brglm2) #biased-reduced regression 

#library(gratia) #for estimating inflection points

library(tidybayes)

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

setwd("~/Documents/GitHub/Galleria_Survival")

data    <- read.table("data/AliveDead.csv", header = T, sep = ",", dec =".")
burden  <- read.table("data/bacterial_burden.csv", header = T, sep = ",", dec =".")
health  <- read.table("data/health_assesment.csv", header = T, sep = ",", dec =".")
death   <- read.table("data/time_to_death.csv", header = T, sep = ",", dec =".")
control <- read.table("data/control_survival.csv", header = T, sep = ",", dec =".")

# Calculate survival proportion
data <- data %>%
  mutate(
    total = total - 2,
    alive = alive - 2,
    survival = alive / total)

LT50 <- approx(x = data$survival, y = data$time, xout = 0.5)$y #22h

fit <- survfit(Surv(time, survival) ~ 1, data = data)
summary(fit)

control <- control %>%
  mutate(survival = alive/total)

names(data)[3:4] <- c("Dead", "Alive")  
names(control)[3:4] <- c("Dead", "Alive") 

# Gompertz model function
gompertz_model <- function(t, a, b) {
  exp(-a / b * (exp(b * t) - 1))
}

# Fit the Gompertz model
fit <- nlsLM(
  survival ~ gompertz_model(time, a, b),
  data = data,
  start = list(a = 1, b = 0.1)  # Initial guesses for parameters
)

# Extract fitted parameters
params <- coef(fit)
a <- params["a"]
b <- params["b"]

LT50 <- (1 / b) * log(1 + (log(2) * b / a)) #22

# Add the fitted Gompertz curve to the data
data$fitted_survival <- gompertz_model(data$time, a, b)

# Create a new column to define the manual grouping for legend
data$group <- "Infected"  
control$group <- "Control"  

# Combine both datasets
plot_data <- bind_rows(data, control)

# Define the Gompertz fit as a separate dataset for the legend
gompertz_fit <- data.frame(time = data$time, fitted_survival = data$fitted_survival, group = "Gompertz fit")

# Now plot with the correct legend mapping
surv <- ggplot(plot_data, aes(x = time, y = survival)) +
  geom_point(data = control, aes(y = survival, fill = group), size = 2, shape = 21) +
  annotate("segment", x = 1, xend = 48, y = 1, yend = 1, color = "grey50", size = 1, linetype = "dashed") +
  geom_point(data = data, aes(y = survival, fill = group), size = 2, shape = 21) +
  geom_line(data = gompertz_fit, aes(x = time, y = fitted_survival, linetype = group, color = group), linewidth = 1) +
  labs(
    x = "Time (h)",
    y = "Survival Proportion",
    fill = "Treatment",  # Legend title for points
  ) +
  mytheme +
  scale_fill_manual(
    values = c("Control" = "grey70", "Infected" = "grey20"),
    labels = c("Control", "Infected")
  ) +
  scale_linetype_manual(
    values = c("Gompertz fit" = "solid"),
    labels = c("Gompertz fit")
  ) +
  
  scale_color_manual(
    values = c("Gompertz fit" = "black"),
    labels = c("Gompertz fit")
  ) +
  guides(
    fill = guide_legend(order = 1),  
    linetype = guide_legend(order = 2), 
    color = "none" 
  ) +
  theme(legend.position = c(.7, .7))  

#ggsave("figures/survival.pdf", plot = surv, width = 5.5, height = 5.5, units = "in", dpi = 300)

# Bacterial Burden 
# tidy plate counts 
burden_tidy <- burden %>%
  pivot_longer(cols = rep1:rep4, names_to = "replicates", values_to = "cfu") %>%
  filter(! Countable == "no") %>%
  filter(!(Sample == "S11" & Larvae %in% c("L5", "L6"))) %>%
  mutate(dilution_factor = (Dilution/Volume_ul)*Buffer_dilution_factor) %>%
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
  mutate(melanization = 4 - melanization, 
         log_CFU = log10(cfu+1), 
         scaled_time = scale(Time), 
         scaled_health = scale(total_score), 
         scaled_activity = scale(activity), 
         scaled_melanization = scale(melanization), 
         scaled_cfu = scale(log_CFU),
         scaled_immune_morb = scale(activity + melanization)
         ) 
  
burden_clean <- burden_tidy %>%
  ungroup() %>%
  filter(cfu > 0, cfu < 1e6) %>%
  filter(Time < 48)

burden_tidy_av <- burden_clean %>% 
  group_by(Time) %>% 
  reframe(cfu = mean(cfu, na.rm = T), log_cfu = mean(log_CFU))

linear_model <- lm(cfu ~ Time, data = burden_clean)
summary(linear_model)
AIC(linear_model)

# Fit logistic model without t0 on raw cfu values 
logistic_model_raw <- nlsLM(cfu ~ K / (1 + ((K - p0) / p0) * exp(-r * Time)),
                            data = burden_clean,
                            start = list(K = max(burden_tidy_av$cfu, na.rm = TRUE),
                                         p0 = min(burden_tidy_av$cfu, na.rm = TRUE),
                                         r = 0.1),
                            control = list(maxiter = 1000))

# log fit
logistic_logfit <- nlsLM(log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
                         data = burden_clean,
                         start = list(K = max(burden_tidy_av$cfu, na.rm = TRUE),
                                      p0 = min(burden_tidy_av$cfu, na.rm = TRUE),
                                      r = 0.1),
                         lower = c(K = 1000, p0 = 0, r = 0.01),
                         upper = c(K = 1e6, p0 = 100, r = 2),
                         control = list(maxiter = 1000))

# View model summary
summary(logistic_logfit)
AIC(logistic_logfit)
rss <- sum(residuals(logistic_logfit)^2)
tss <- sum((burden_clean$log_CFU - mean(burden_clean$log_CFU))^2)
pseudo_r2 <- 1 - rss / tss

pred_time <- data.frame(Time = seq(min(burden_clean$Time), max(burden_clean$Time), by = 0.1))

boot_preds <- function(data, indices) {
  d <- data[indices, ]
  fit <- tryCatch({
    nlsLM(
      log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
      data = d,
      start = list(
        K = quantile(d$cfu, 0.9, na.rm = TRUE),
        p0 = quantile(d$cfu, 0.1, na.rm = TRUE),
        r = 0.2
      ),
      lower = c(K = 100, p0 = 0.1, r = 0.01),
      upper = c(K = 1e6, p0 = 5e3, r = 2),
      control = list(maxiter = 1000)
    )
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(rep(NA, 3))
  coef(fit)  # returns K, p0, r
}

set.seed(123)
boot_curve <- boot(data = burden_clean, statistic = boot_preds, R = 1000)

# Extract the K estimates (assuming they are the first column)
boot_K <- boot_curve$t[, 1]
# Remove failed fits (NAs)
boot_K <- boot_K[!is.na(boot_K)]
# Compute 95% confidence interval using percentiles
ci_K <- quantile(boot_K, probs = c(0.025, 0.975), na.rm = TRUE)
print(ci_K)

boot_matrix <- boot_curve$t
valid_matrix <- boot_matrix[apply(boot_matrix, 1, function(row) {
  !any(is.na(row)) && (max(row) - min(row)) > 1  # only keep replicates with real slope
}), ]

# Generate predictions from bootstrap parameter sets
boot_prediction_curves <- apply(valid_matrix, 1, function(params) {
  K <- params[1]; p0 <- params[2]; r <- params[3]
  log10(K / (1 + ((K - p0) / p0) * exp(-r * pred_time$Time)))
})

# boot_prediction_curves is now a matrix: rows = timepoints, cols = bootstrap replicates
boot_prediction_curves <- t(boot_prediction_curves)  # transpose so timepoints are columns

# Calculate CI and median fit at each timepoint
ci_bounds <- apply(boot_prediction_curves, 2, quantile, probs = c(0.025, 0.975))
median_fit <- apply(boot_prediction_curves, 2, median)

# Add to pred_time
pred_time$lower <- ci_bounds[1, ]
pred_time$upper <- ci_bounds[2, ]
pred_time$fit   <- median_fit

# Plot
pred_time$group <- "Logistic fit" 

bact <- ggplot(burden_clean, aes(x = Time, y = log_CFU)) +
  geom_jitter(aes(fill = status), size = 2, shape = 21, alpha = 1, color = "black") +
  geom_ribbon(data = pred_time, aes(x = Time, ymin = lower, ymax = upper), 
               fill = "gray80", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = pred_time, aes(x = Time, y = fit, linetype = group), 
            color = "black", size = 1.2, inherit.aes = FALSE) +
  xlab("Time (h)") +
  ylab("log(Colony Forming Units)") +
  scale_fill_manual(values = c("Alive" = "#19798b", "Dead" = "#ee9b43")) +
  guides(fill = guide_legend(order = 1)) +
  theme_minimal() +
  theme(legend.position = c(.2, .8)) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32, 36)) +
  mytheme

#ggsave("figures/bacterial_burden.pdf", plot = bact, width = 5.5, height = 5.5, units = "in", dpi = 300)

#########################################################################################################
# Health variables 

# Fit ordinal regression
fit <- polr(factor(melanization) ~ Time, data = burden_clean, method = "logistic")

# Null model (intercept only)
fit_null_mel <- polr(factor(melanization) ~ 1, data = burden_clean, method = "logistic")

# Calculate McFadden's pseudo-R²
resid_dev_mel <- deviance(fit)
null_dev_mel <- deviance(fit_null_mel)

pseudo_r2_mel <- 1 - (resid_dev_mel / null_dev_mel)

# Get predicted probabilities for each category (0-4)
probabilities <- predict(fit, type = "probs")

# Compute expected melanization (weighted sum of probability)
expected_melanization <- rowSums(probabilities * matrix(rep(0:4, each = nrow(probabilities)), 
                                                        nrow = nrow(probabilities), byrow = FALSE))

# Reverse predicted melanization scale too
burden_clean$predicted_melanization <- expected_melanization

# Fit GAM for smooth prediction
gam_fit <- gam(predicted_melanization ~ s(Time, k = 5), data = burden_clean)
summary(gam_fit)

# Generate smooth predictions
time_seq <- data.frame(Time = seq(min(burden_clean$Time), max(burden_clean$Time), length.out = 100))
time_seq$smoothed_melanization <- predict(gam_fit, newdata = time_seq)

# Plot
# Create a dataset for the GAM fit with a label
time_seq$group <- "GAM fit" 

melanization_plot <- ggplot(burden_clean, aes(x = Time, y = melanization)) +
  geom_jitter(aes(fill = status), size = 2, shape = 21, alpha = 1, color = "black") +
  geom_line(data = time_seq, aes(x = Time, y = smoothed_melanization, linetype = group, color = group), size = 1) +
  scale_fill_manual(values = c("#19798b", "#ee9b43")) +
  scale_linetype_manual(
    values = c("GAM fit" = "solid"),
    labels = c("GAM fit")
  ) +
  scale_color_manual(
    values = c("GAM fit" = "black"),
    labels = c("GAM fit")
  ) +
  guides(
    fill = guide_legend(order = 1), 
    linetype = guide_legend(order = 2),  
    color = "none"  
  ) + 
  labs(x = "Time (h)", y = "Melanization") +  
  mytheme +
  theme(legend.position = c(.8, .2))  

# Fit ordinal regression
fit_act <- polr(factor(activity) ~ Time, data = burden_clean, method = "logistic")

# Null model (intercept only)
fit_null <- polr(factor(activity) ~ 1, data = burden_clean, method = "logistic")

# Calculate McFadden's pseudo-R²
resid_dev <- deviance(fit_act)
null_dev <- deviance(fit_null)

pseudo_r2 <- 1 - (resid_dev / null_dev)
pseudo_r2

# Get predicted probabilities for each category (0-4)
probabilities_act <- predict(fit_act, type = "probs")

# Compute expected activity (weighted sum of probability)
expected_activity <- rowSums(probabilities_act * matrix(rep(0:3, each = nrow(probabilities_act)), 
                                                        nrow = nrow(probabilities_act), byrow = FALSE))

# Reverse predicted melanization scale too
burden_clean$predicted_activity <- expected_activity

# Fit GAM for smooth prediction
gam_fit_act <- gam(predicted_activity ~ s(Time, k = 5), data = burden_clean)
summary(gam_fit_act)

# Generate smooth predictions
time_seq_act <- data.frame(Time = seq(min(burden_clean$Time), max(burden_clean$Time), length.out = 100))
time_seq_act$smoothed_activity <- predict(gam_fit_act, newdata = time_seq_act)

# Plot
activity_plot <- ggplot(burden_clean, aes(x = Time, y = activity)) +
  geom_jitter(aes(fill = status), size = 2, shape = 21, alpha = 1, color = "black") +
  geom_line(data = time_seq_act, aes(x = Time, y = smoothed_activity), 
            color = "black", size = 1) +  # Smoothed curve
  scale_fill_manual(values = c("#19798b", "#ee9b43")) +
  labs(x = "Time (h)", y = "Activity") +  
  mytheme +
  theme(legend.position = "none")

# Arrange the two small plots together first
small_plots <- plot_grid(surv, activity_plot, melanization_plot, ncol = 3, align = "h", labels = "AUTO")
ggsave("figures/small_plots.pdf", plot = small_plots, width = 13, height = 4, units = "in", dpi = 300)

#################
# Test how interpretations change after selectively fitting the models to living individuals only 

# Subset data to only living individuals
burden_tidy_alive <- burden_clean %>% filter(status == "Alive")
burden_av_alive   <- burden_tidy_alive %>% 
  group_by(Time) %>% 
  reframe(cfu = mean(cfu, na.rm = T), log_cfu = mean(log_CFU))

linear_model_alive <- lm(log10(cfu) ~ log10(Time), data = burden_tidy_alive)
summary(linear_model_alive)
AIC(linear_model_alive)

# log fit
logistic_logfit_a <- nlsLM(log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
                         data = burden_tidy_alive,
                         start = list(K = max(burden_av_alive$log_cfu, na.rm = TRUE),
                                      p0 = min(burden_av_alive$log_cfu, na.rm = TRUE),
                                      r = 0.1),
                         control = list(maxiter = 1000))

# View model summary
summary(logistic_logfit_a)
AIC(logistic_logfit_a)
rss <- sum(residuals(logistic_logfit_a)^2)
tss <- sum((burden_tidy_alive$log_CFU - mean(burden_tidy_alive$log_CFU))^2)
pseudo_r2 <- 1 - rss / tss

pred_time_a <- data.frame(Time = seq(min(burden_tidy_alive$Time), max(burden_tidy_alive$Time), by = 0.1))

boot_preds <- function(data, indices) {
  d <- data[indices, ]
  fit <- tryCatch({
    nlsLM(
      log_CFU ~ log10(K / (1 + ((K - p0) / p0) * exp(-r * Time))),
      data = d,
      start = list(
        K = quantile(d$cfu, 0.9, na.rm = TRUE),
        p0 = quantile(d$cfu, 0.1, na.rm = TRUE),
        r = 0.2
      ),
      lower = c(K = 100, p0 = 0.1, r = 0.01),
      upper = c(K = 1e5, p0 = 5e3, r = 2),
      control = list(maxiter = 1000)
    )
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(rep(NA, 3))
  coef(fit)  # returns K, p0, r
}

set.seed(124)
boot_curve_a <- boot(data = burden_tidy_alive, statistic = boot_preds, R = 1000)

# Extract the K estimates (assuming they are the first column)
boot_K <- boot_curve_a$t[, 1]
# Remove failed fits (NAs)
boot_K <- boot_K[!is.na(boot_K)]
# Compute 95% confidence interval using percentiles
ci_K <- quantile(boot_K, probs = c(0.025, 0.975), na.rm = TRUE)
print(ci_K)

boot_matrix_a <- boot_curve_a$t
valid_matrix_a <- boot_matrix_a[apply(boot_matrix_a, 1, function(row) {
  !any(is.na(row)) && (max(row) - min(row)) > 1  # only keep replicates with real slope
}), ]

# Generate predictions from bootstrap parameter sets
boot_prediction_curves_a <- apply(valid_matrix_a, 1, function(params) {
  K <- params[1]; p0 <- params[2]; r <- params[3]
  log10(K / (1 + ((K - p0) / p0) * exp(-r * pred_time_a$Time)))
})

# boot_prediction_curves is now a matrix: rows = timepoints, cols = bootstrap replicates
boot_prediction_curves_a <- t(boot_prediction_curves_a)  # transpose so timepoints are columns

# Calculate CI and median fit at each timepoint
ci_bounds <- apply(boot_prediction_curves_a, 2, quantile, probs = c(0.025, 0.975))
median_fit <- apply(boot_prediction_curves_a, 2, median)

# Add to pred_time
pred_time_a$lower <- ci_bounds[1, ]
pred_time_a$upper <- ci_bounds[2, ]
pred_time_a$fit   <- median_fit


# Plot
bact2 <- ggplot(burden_clean, aes(x = Time, y = log_CFU)) +
  geom_jitter(aes(fill = status), size = 2, shape = 21, alpha = 1, color = "black") +
  geom_ribbon(data = pred_time, aes(x = Time, ymin = lower, ymax = upper), 
              fill = "#ee9b43", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = pred_time, aes(x = Time, y = fit, linetype = group), 
            color = "#b80422", size = 1.2, inherit.aes = FALSE) +
  geom_ribbon(data = pred_time_a, aes(x = Time, ymin = lower, ymax = upper), 
              fill = "#19798b", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = pred_time_a, aes(x = Time, y = fit), 
            color = "#19798b", linetype = "dashed", size = 1.2, inherit.aes = FALSE) +
  xlab("Time (h)") +
  ylab("log(Colony Forming Units)") +
  scale_fill_manual(values = c("Alive" = "#19798b", "Dead" = "#ee9b43")) +
  guides(fill = guide_legend(order = 1)) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32, 36)) +
  mytheme+
  theme(legend.position = c(.2, .8))

#ggsave("figures/bacterial_burden_aliveDead.pdf", plot = bact2, width = 5.5, height = 5, units = "in", dpi = 300)

# Pathogen growth vs. health measures

cfuMel <- ggplot(burden_clean, aes(y = melanization, x = log10(cfu))) +
  geom_jitter(size = 2, shape =21, color = "black", aes(fill = status)) +
  mytheme +
  #stat_smooth(method = "gam", color = "grey25") + 
  xlab("log(CFU)") +
  ylab("Melanization")+
  scale_fill_manual(values = c("#19798b", "#ee9b43"))+
  theme(legend.position = "none")

cfuAct <- ggplot(burden_clean, aes(y = activity, x = log10(cfu))) +
  geom_jitter(size = 2, shape = 21, color = "black", aes( fill = status)) +
  mytheme +
  #stat_smooth(method = "gam", color = "grey25") + 
  xlab("log(CFU)") +
  ylab("Activity")+
  scale_fill_manual(values = c("#19798b", "#ee9b43"))+
  theme(legend.position = "none")

cfuActMel <- ggplot(burden_clean, aes(x = melanization, y = activity)) +
  geom_jitter(size = 2, shape = 21, color = "black", aes( fill = status)) +
  mytheme +
  #stat_smooth(method = "gam", color = "grey25") + 
  ylab("Activity") +
  xlab("Melanization")+
  scale_fill_manual(values = c("#19798b", "#ee9b43"))+
  theme(legend.position = "none")

small <- plot_grid(cfuAct, cfuMel, ncol = 1, align = "v", labels = c("C", "D"), hjust = -0.5, vjust=1.5)
big <- plot_grid(bact2, small, ncol = 2, align = "v", labels = c("A", ""), rel_widths = c(2, 1.2), hjust = -0.5, vjust=1.5)
#ggsave("figures/bacterial_burden_aliveDead.pdf", plot = big, width = 9, height = 4.5, units = "in", dpi = 300)

######DAG######

dag <- dagitty("dag {
  scaled_time   -> scaled_cfu
  scaled_time   -> scaled_health
  scaled_cfu    -> scaled_health
  scaled_health -> survival
  scaled_time   -> survival
  scaled_cfu    -> survival
  scaled_health <-> scaled_cfu
}")

plot(dag)

# Set coordinates for a clean layout
coordinates(dag) <- list(
  x = c(scaled_time = 0, scaled_cfu = -1, scaled_health = 0, survival = 0),
  y = c(scaled_time = 3, scaled_cfu = 2, scaled_health = 1, survival = 0)
)

ggdag(dag)+
  theme_dag_blank()+
  geom_dag_point(size = 20, shape = 21, fill = "white")+
  geom_dag_text(color = "black")


# Systematical testing 
summary(lm(scaled_cfu    ~ scaled_time, data = burden_clean))
summary(lm(scaled_health ~ scaled_cfu + scaled_time, data = burden_clean))

# bidirection?
summary(lm(scaled_cfu    ~ scaled_health + scaled_time, data = burden_clean))

summary(glm(survival     ~ scaled_time, data = burden_clean, family = binomial(), method = "brglmFit"))
summary(glm(survival     ~ scaled_health, family = binomial, data = burden_clean, method = "brglmFit"))
summary(glm(survival     ~ scaled_health + scaled_time, data = burden_clean, family = binomial(), method = "brglmFit"))
summary(glm(survival     ~ scaled_cfu, data = burden_clean, family = binomial(), method = "brglmFit"))
summary(glm(survival     ~ scaled_cfu + scaled_time, data = burden_clean, family = binomial(), method = "brglmFit"))
summary(glm(survival     ~ scaled_health + scaled_cfu + scaled_time, family = binomial(), data = burden_clean, method = "brglmFit"))

#Nicer Dag diagram 

# Define your DAG with labels and coordinates
dag_final <- dagify(
  survival ~ scaled_health +scaled_cfu + scaled_time,
  scaled_health ~ scaled_cfu + scaled_time,
  scaled_cfu ~ scaled_time,
  exposure = "scaled_cfu",
  outcome = "survival",
  labels = c(
    scaled_cfu = "Pathogen load",
    scaled_health = "Host health",
    scaled_time = "Time",
    survival = "Survival"
  ),
  coords = list(
    x = c(scaled_time = 0, scaled_cfu = 3, scaled_health = 1, survival = 2),
    y = c(scaled_time = 3, scaled_cfu = 2, scaled_health = 1, survival = 0)
  )
)

#%>%tidy_dagitty()


# Create the plot
ggdag_status(dag_final, use_labels = "label", text = FALSE) +
  theme_dag_blank() +
  theme(
    text = element_text(size = 14)
  ) +
  guides(fill = FALSE, color = FALSE)


# SEM 

mediation_model <- '
  scaled_cfu ~ scaled_time
  scaled_health ~ a * scaled_cfu + scaled_time
  survival ~ b * scaled_health + c_prime * scaled_cfu + scaled_time
  indirect := a * b
  total := c_prime + (a * b)
'

burden_clean$time_linear <- burden_clean$scaled_time
burden_clean$time_squared <- burden_clean$scaled_time^2

model_quad <- '
  scaled_cfu ~ t1 * time_linear + t2 * time_squared
  scaled_health ~ a * scaled_cfu + h1 * time_linear 
  survival ~ b * scaled_health + c_prime * scaled_cfu

  # Mediation effect
  indirect := a * b
  total := c_prime + (a * b)
'

fit_bayes_sem1 <- bsem(model_quad, data = burden_clean, burnin = 1000, sample = 5000, n.chains = 3, 
                       save.lvs = TRUE)
summary(fit_bayes_sem1)

lavaanPlot(model = fit_bayes_sem1, 
           labels = list(
             scaled_time = "Time",
             scaled_cfu = "Pathogen load",
             scaled_health= "Host health",
             survival = "Survival"
           ),
           coefs = TRUE,
           node_options = list( fontname = "Helvetica"), 
           edge_options = list(color = "grey"), stars = c("regress"))  

semPaths(fit_bayes_sem1,
         what = "std",
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

# Check for issues
summary(fit_bayes_sem1)
blavInspect(fit_bayes_sem1, "psrf")  # Rhat values
pp_check(fit_bayes_sem1)
resid(fit_bayes_sem1, type = "standardized")

effects_df <- data.frame(
  path = c("8Time -> Pathogen (t1)", 
           "7Time² -> Pathogen (t2)",
           "6Time -> Health (h1)", 
           "5Health -> Survival (b)",
           "4Pathogen -> Health (a)",
           "3Indirect (a × b)", 
           "2Pathogen -> Survival (c')",
           "1Total (a × b + c')"),
  estimate = c(0.905, -0.425, -0.609, 0.515, -0.434, -0.223, 0.037, -0.186),
  lower =    c(0.79, -0.607, -0.866, 0.456, -0.643, -0.333,  -0.021,  --0.300),
  upper =    c(1.050, -0.246, -0.341, 0.572, -0.228, -0.114, 0.095, -0.072)
)

effect_plot <- ggplot(effects_df, aes(y = path, x = estimate)) +
  geom_col(fill = "grey", width = 0.6) +
  geom_vline(xintercept = 0, linetype ="dashed")+
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "", x = "Standardized Effect Size") +
  mytheme

ggsave("figures/effect_plot.pdf", plot = effect_plot, width = 6.5, height = 5.5, units = "in", dpi = 300)

#############################################################################################
# Mapping

normalize <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

a     = 4.139139e-05 
b     = 0.4005271 
p     = 10^pred_time$fit
p_log = pred_time$fit
p0    <- 4.449848e-01   # since original p0 is in log10
K     <- 5.087321e+04   # K already in CFU
r     <- 4.523319e-01

time <- seq(1, 36, length.out = length(p))  # match length
m_t = a*exp(b*time)
m_p <- a * ((p * (p0 - K)) / (p0 * (p - K)))^(b / r)

df <- tibble(
  time = time,
  m_t = m_t,
  p = p,
  m_p = m_p,
  time_norm = normalize(time), 
  m_t_norm = normalize(m_t), 
  p_norm = normalize(p_log), 
  m_p_norm = normalize(m_p)
)

# Plot all 3 panels
p1 <- ggplot(df, aes(x = time_norm, y = p_norm)) +
  geom_line(color = "#19798b", size = 1.2) +
  labs(x = "Time (unit)", y = "Scaled p(t)")+
  mytheme

p2 <- ggplot(df, aes(x = time_norm, y = m_t_norm)) +
  geom_line(color = "#ee9b43", size = 1.2) +
  labs(x = "Time (unit)", y = "Scaled m(t)")+
  mytheme

p3 <- ggplot(df, aes(x = p_norm, y = m_p_norm)) +
  geom_line(color = "#e76f51", size = 1.2) +
  labs(x = "Scaled p(t)", y = "Scaled m(p)")+
  mytheme

plot_grid(p1, p2, p3, labels = c("A", "B", "C"), nrow = 1)

# m_t 

# Data to overlay 
data_m_t <- data %>%
  arrange(time) %>%
  filter(time < 37) %>%
  mutate(
    m_t = log(lag(survival) / survival),
    m_t = ifelse(m_t < 0, NA, m_t),
    time_m = (lag(time) + time) / 2,
    m_t_norm = normalize(m_t), 
    time_norm = normalize(time_m)
  )

#p_t
p_t_data <- tibble(p_t = normalize(burden_clean$log_CFU), 
                   time = burden_clean$Time)

m_t_data  <- tibble(m_t = data_m_t$m_t_norm, 
                    time_norm = data_m_t$time_norm)

m_p_data <- tibble(m_p = df$m_p_norm, 
                   time_norm = df$time_norm, 
                   p_t = df$p_norm)

# Summarize pathogen burden at each time point (mean log10 CFU)
burden_summary <- burden_clean %>%
  group_by(Time) %>%
  reframe(log_CFU_mean = mean(log_CFU, na.rm = TRUE), .groups = "drop", 
          activity = mean(activity, na.rm = TRUE), 
          melanization = mean(melanization, na.rm = TRUE), 
          time = mean(Time, na.rm = TRUE))

p4 <- ggplot(df, aes(x = time_norm, y = p_norm)) +
  geom_line(color = "#19798b", size = 1.2) +
  geom_jitter(data = burden_summary, aes(x = normalize(Time), y = normalize(log_CFU_mean)),
             inherit.aes = FALSE, fill = "#19798b",size = 3, shape = 21, alpha = 0.7) +
  labs(x = "Time (unit)", y = "Scaled p(t)") +
  mytheme


p5 <- ggplot(df, aes(x = time_norm, y = m_t_norm)) +
  geom_line(color = "#ee9b43", size = 1.2) +
  geom_point(data = m_t_data, aes(x = time_norm, y = m_t),
             inherit.aes = FALSE, fill = "#ee9b43", size = 3, shape = 21, alpha = 0.7) +
  
  labs(x = "Time (unit)", y = "Scaled m(t)") +
  mytheme

# Interpolate both datasets over a common time vector

# Common time vector (same as for m(t))
time_common <- seq(1, 36, by = 1)

# Interpolate burden means across full time range
p_interp <- approx(x = burden_summary$Time, y = burden_summary$log_CFU_mean, xout = time_common)$y

# Interpolate mortality estimates from m(t)
m_interp <- approx(x = data_m_t$time_m, y = data_m_t$m_t, xout = time_common)$y

# Normalize both
normalize <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

p_norm <- normalize(p_interp)
m_norm <- normalize(m_interp)

# Create overlay dataframe
map_overlay <- tibble(
  time = time_common,
  p = p_norm,
  m = m_norm
)

# Plot with model-derived m(p) mapping
p6 <- ggplot(map_overlay, aes(x = p, y = m)) +
  geom_point(shape = 21, fill = "grey", size = 3, alpha) +
  geom_line(data = df, aes(x = p_norm, y = m_p_norm), color = "#b80422", size = 1.2, inherit.aes = FALSE) +
  labs(x = "Scaled p(t)", y = "Scaled m(t)") +
  mytheme

plot_grid(p4, p5, p6, labels = c("A", "B", "C"), nrow = 1)


# Interpolate melanization and activity
mel_interp <- approx(x = burden_summary$Time, y = burden_summary$melanization, xout = time_common)$y
act_interp <- approx(x = burden_summary$Time, y = burden_summary$activity, xout = time_common)$y

mel_norm <- normalize(mel_interp)
act_norm <- normalize(act_interp)

map_overlay <- tibble(
  time = time_common,
  p = p_norm,
  m = m_norm,
  activity = act_interp,
  melanization = mel_interp,
  raw_time = time_common
)

# Color points by normalized time
ggplot(map_overlay, aes(x = p, y = m, fill = raw_time)) +
  geom_point(shape = 21, size = 3, color = "black") +
  geom_line(data = df, aes(x = p_norm, y = m_p_norm), 
            inherit.aes = FALSE, color = "#b80422", size = 1.2) +
  scale_fill_viridis_c(option = "B") +
  labs(x = "Scaled p(t)", y = "Scaled m(t)", fill = "Time (h)") +
  mytheme


# Color points by normalized time
mp1 <- ggplot(map_overlay, aes(x = p, y = m, fill = raw_time)) +
  geom_point(shape = 21, size = 3, color = "black") +
  geom_line(data = df, aes(x = p_norm, y = m_p_norm), 
            inherit.aes = FALSE, color = "#b80422", size = 1.2) +
  scale_fill_viridis_c(option = "B", direction = -1) +
  labs(x = "Scaled p(t)", y = "Scaled m(t)", fill = "Time (h)") +
  mytheme

mp2 <- ggplot(map_overlay, aes(x = p, y = m, fill = activity)) +
  geom_point(shape = 21, size = 3, color = "black") +
  geom_line(data = df, aes(x = p_norm, y = m_p_norm), 
            inherit.aes = FALSE, color = "#b80422", size = 1.2) +
  scale_fill_viridis_c(option = "B") +
  labs(x = "Scaled p(t)", y = "Scaled m(t)", fill = "Time (h)") +
  mytheme

mp3 <- ggplot(map_overlay, aes(x = p, y = m, fill = melanization)) +
  geom_point(shape = 21, size = 3, color = "black") +
  geom_line(data = df, aes(x = p_norm, y = m_p_norm), 
            inherit.aes = FALSE, color = "#b80422", size = 1.2) +
  scale_fill_viridis_c(option = "B", direction = -1) +
  labs(x = "Scaled p(t)", y = "Scaled m(t)", fill = "Time (h)") +
  mytheme


plot_grid(mp1, mp2, mp3, labels = c("A", "B", "C"), nrow = 1)
