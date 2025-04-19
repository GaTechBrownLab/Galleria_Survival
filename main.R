# Load required libraries
library(minpack.lm)  # For non-linear regression
library(ggplot2)   
library(dplyr)
library(tidyr)
library(boot)
library(MASS)
library(mgcv)
library(cowplot)

# SEM
library(lavaan)
library(blavaan) #bayesian lavaan
library(lavaanPlot)
library(semPlot)

#DAG
library(dagitty)
library(ggdag)

library(brglm2) #biased-reduced regression 

library(gratia) #for estimating inflection points

library(tidybayes)

# Better fonts 
install.packages("extrafont")
library(extrafont)
extrafont::font_import()   # only once — can take a while!
extrafont::loadfonts()


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
         mutate(total = total - 2, 
         dead = dead - 2, 
         alive = alive - 2, 
         survival = alive/total)

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

LT50 <- (1 / b) * log(1 + (log(2) * b / a))

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
  geom_line(data = gompertz_fit, aes(x = time, y = fitted_survival, linetype = group, color = group), size = 1) +
  labs(
    x = "Time (h)",
    y = "Survival Proportion",
    fill = "Treatment",  # Legend title for points
  ) +
  mytheme +
  scale_fill_manual(
    values = c("Control" = "#C993A2", "Infected" = "#94aec2"),
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

linear_model <- lm(log_CFU ~ Time, 
                      data = burden_tidy)
summary(linear_model)
AIC(linear_model)

# Define the logistic function
logistic_model <- nls(log_CFU ~ K / (1 + exp(-r * (Time - t0))), 
                      data = burden_tidy,
                      start = list(K = max(burden_tidy$log_CFU), 
                                   r = 0.1, 
                                   t0 = mean(burden_tidy$Time)),
                      control = list(maxiter = 500))

# View model summary
summary(logistic_model)
AIC(logistic_model)
rss <- sum(residuals(logistic_model)^2)
tss <- sum((burden_tidy$log_CFU - mean(burden_tidy$log_CFU))^2)
pseudo_r2 <- 1 - rss / tss

# Generate more time points for a smooth curve
time_seq_all <- data.frame(Time = seq(min(burden_tidy$Time), max(burden_tidy$Time), by = 0.1))

# Predict values at these time points
time_seq_all$log_CFU_fitted <- predict(logistic_model, newdata = time_seq_all)

# Define function for bootstrapping
boot_nls <- function(data, indices) {
  df_boot <- data[indices, ]  # Resample with replacement
  fit <- try(nls(log_CFU ~ K / (1 + exp(-r * (Time - t0))),
                 data = df_boot,
                 start = list(K = max(burden_tidy$log_CFU), 
                              r = 0.1, 
                              t0 = mean(burden_tidy$Time))), 
             silent = TRUE)
  
  if (inherits(fit, "try-error")) return(rep(NA, 3))  # Skip failed fits
  return(coef(fit))  # Return model parameters
}

# Run bootstrapping (adjust number of iterations if needed)
set.seed(42)
boot_results <- boot(burden_tidy, boot_nls, R = 1000)

# Store predicted values for each bootstrap iteration
boot_predictions <- matrix(NA, nrow = length(time_seq_all$Time), ncol = nrow(boot_results$t))

# Generate predictions for each bootstrap sample
for (i in 1:nrow(boot_results$t)) {
  K_i <- boot_results$t[i, 1]  
  r_i <- boot_results$t[i, 2]  
  t0_i <- boot_results$t[i, 3]  
  boot_predictions[, i] <- K_i / (1 + exp(-r_i * (time_seq_all$Time - t0_i)))
}

# Compute the final fitted curve using median
time_seq_all$log_CFU_fitted_median <- apply(boot_predictions, 1, median, na.rm = TRUE)

# Compute 95% confidence intervals at each time point
time_seq_all$log_CFU_lower <- apply(boot_predictions, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
time_seq_all$log_CFU_upper <- apply(boot_predictions, 1, function(x) quantile(x, 0.975, na.rm = TRUE))

# Plot with corrected confidence intervals
# Create a label in time_seq for the logistic fit
time_seq_all$group <- "Logistic fit"

bact <- ggplot(burden_tidy, aes(x = Time, y = log_CFU)) +
  geom_jitter(aes(fill = status), size = 2, shape = 21, alpha = 1, color = "black") +
  mytheme +
  geom_ribbon(data = time_seq_all, aes(x = Time, ymin = log_CFU_lower, ymax = log_CFU_upper), 
              inherit.aes = FALSE, fill = "grey", alpha = 0.3) + 
  geom_line(data = time_seq_all, aes(x = Time, y = log_CFU_fitted, linetype = group, color = group), 
            inherit.aes = FALSE, size = 1) +  
  xlab("Time (h)") +
  ylab("log(Colony Forming Units)") +
  scale_fill_manual(values = c("#19798b", "#ee9b43")) +
  scale_linetype_manual(
    values = c("Logistic fit" = "solid"),
    labels = c("Logistic fit")
  ) +
  scale_color_manual(
    values = c("Logistic fit" = "black"),
    labels = c("Logistic fit")
  ) +
  guides(
    fill = guide_legend(order = 1),  
    linetype = guide_legend(order = 2),
    color = "none"
  ) +  
  theme(legend.position = c(.2, .8)) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32, 36, 48))

ggsave("figures/bacterial_burden.pdf", plot = bact, width = 5.5, height = 5.5, units = "in", dpi = 300)

# Fit ordinal regression
fit <- polr(factor(melanization) ~ Time, data = burden_tidy, method = "logistic")

# Null model (intercept only)
fit_null_mel <- polr(factor(melanization) ~ 1, data = burden_tidy, method = "logistic")

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
burden_tidy$predicted_melanization <- expected_melanization

# Fit GAM for smooth prediction
gam_fit <- gam(predicted_melanization ~ s(Time, k = 5), data = burden_tidy)
summary(gam_fit)

# Generate smooth predictions
time_seq <- data.frame(Time = seq(min(burden_tidy$Time), max(burden_tidy$Time), length.out = 100))
time_seq$smoothed_melanization <- predict(gam_fit, newdata = time_seq)

# Predict derivatives
fd_mel <- derivatives(gam_fit)

# Find inflection point as time of steepest slope:
inflection_time_mel <- fd_mel$Time[which.max(abs(fd_mel$.derivative))]

# Plot
# Create a dataset for the GAM fit with a label
time_seq$group <- "GAM fit" 

melanization_plot <- ggplot(burden_tidy, aes(x = Time, y = melanization)) +
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
fit_act <- polr(factor(activity) ~ Time, data = burden_tidy, method = "logistic")

# Null model (intercept only)
fit_null <- polr(factor(activity) ~ 1, data = burden_tidy, method = "logistic")

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
burden_tidy$predicted_activity <- expected_activity

# Fit GAM for smooth prediction
gam_fit_act <- gam(predicted_activity ~ s(Time, k = 5), data = burden_tidy)
summary(gam_fit_act)

# Predict derivatives
fd <- derivatives(gam_fit_act)

# Find inflection point as time of steepest slope:
inflection_time <- fd$Time[which.max(abs(fd$.derivative))]

# Generate smooth predictions
time_seq_act <- data.frame(Time = seq(min(burden_tidy$Time), max(burden_tidy$Time), length.out = 100))
time_seq_act$smoothed_activity <- predict(gam_fit_act, newdata = time_seq_act)

# Plot
activity_plot <- ggplot(burden_tidy, aes(x = Time, y = activity)) +
  geom_jitter(aes(fill = status), size = 2, shape = 21, alpha = 1, color = "black") +
  geom_line(data = time_seq_act, aes(x = Time, y = smoothed_activity), 
            color = "black", size = 1) +  # Smoothed curve
  scale_fill_manual(values = c("#19798b", "#ee9b43")) +
  labs(x = "Time (h)", y = "Activity") +  
  mytheme +
  theme(legend.position = "none")

# Arrange the two small plots together first
small_plots <- plot_grid(activity_plot, melanization_plot, surv, ncol = 3, align = "h", labels = "AUTO")
ggsave("figures/small_plots.pdf", plot = small_plots, width = 13, height = 4, units = "in", dpi = 300)

#################
# Time of death 
death <- ggplot(burden_tidy, aes(x = time_to_death_h, y = log10(cfu))) +
  geom_point(size = 2, shape =21, alpha = 0.5, color = "black", fill = "#ee9b43") +
  mytheme +
  stat_smooth(method = "lm", color = "grey25") + 
  ylab("log(Colony Forming Units)") +
  xlab("Time of death")

ggsave("figures/time_death.pdf", plot = death, width = 5.5, height = 5.5, units = "in", dpi = 300)

#################
# Test how interpretations change after selectively fitting the models to living individuals only 

# Subset data to only living individuals
burden_tidy_alive <- burden_tidy %>% filter(status == "Alive")

# Fit a linear model to only living individuals
linear_model_alive <- lm(log_CFU ~ Time, data = burden_tidy_alive)
summary(linear_model_alive)

# Define the logistic function
logistic_model_alive <- nls(log_CFU ~ K / (1 + exp(-r * (Time - t0))), 
                      data = burden_tidy_alive,
                      start = list(K = max(burden_tidy_alive$log_CFU), 
                                   r = 0.1, 
                                   t0 = mean(burden_tidy_alive$Time)),
                      control = list(maxiter = 500))

# View model summary
summary(logistic_model_alive)

AIC_linear <- AIC(linear_model_alive)
AIC_logistic <- AIC(logistic_model_alive)

rss2 <- sum(residuals(logistic_model_alive)^2)
tss2 <- sum((burden_tidy_alive$log_CFU - mean(burden_tidy_alive$log_CFU))^2)
pseudo_r2_2 <- 1 - rss2 / tss2

# Plot logistic one 

# Generate more time points for a smooth curve
time_seq_alive <- data.frame(Time = seq(min(burden_tidy_alive$Time), max(burden_tidy_alive$Time), by = 0.1))

# Predict values at these time points
time_seq_alive$log_CFU_fitted <- predict(logistic_model_alive, newdata = time_seq_alive)

# Define function for bootstrapping
boot_nls_alive <- function(data, indices) {
  df_boot <- data[indices, ]  # Resample with replacement
  fit <- try(nls(log_CFU ~ K / (1 + exp(-r * (Time - t0))),
                 data = df_boot,
                 start = list(K = max(burden_tidy_alive$log_CFU), 
                              r = 0.1, 
                              t0 = mean(burden_tidy_alive$Time))), 
             silent = TRUE)
  
  if (inherits(fit, "try-error")) return(rep(NA, 3))  # Skip failed fits
  return(coef(fit))  # Return model parameters
}

# Run bootstrapping (adjust number of iterations if needed)
set.seed(43)
boot_results_alive <- boot(burden_tidy_alive, boot_nls_alive, R = 1000)

# Store predicted values for each bootstrap iteration
boot_predictions_alive <- matrix(NA, nrow = length(time_seq_alive$Time), ncol = nrow(boot_results_alive$t))

# Generate predictions for each bootstrap sample
for (i in 1:nrow(boot_results_alive$t)) {
  K_i <- boot_results_alive$t[i, 1]  
  r_i <- boot_results_alive$t[i, 2]  
  t0_i <- boot_results_alive$t[i, 3]  
  boot_predictions_alive[, i] <- K_i / (1 + exp(-r_i * (time_seq_alive$Time - t0_i)))
}

boot_matrix_alive <- boot_results_alive$t
boot_matrix_alive <- boot_matrix_alive[, colSums(is.na(boot_matrix_alive)) == 0]  # Remove failed bootstrap runs

# Compute 95% confidence intervals at each time point
time_seq_alive$log_CFU_fitted_median <- apply(boot_predictions_alive, 1, median, na.rm = TRUE)
time_seq_alive$log_CFU_lower <- apply(boot_predictions_alive, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
time_seq_alive$log_CFU_upper <- apply(boot_predictions_alive, 1, function(x) quantile(x, 0.975, na.rm = TRUE))

# Plot with corrected confidence intervals
# Create a label in time_seq for the logistic fit
time_seq_alive$group <- "Logistic fit"

bact2 <- ggplot(burden_tidy, aes(x = Time, y = log_CFU)) +
  geom_jitter(aes(fill = status), size = 2, shape = 21, alpha = 1, color = "black") +
  mytheme +
  geom_ribbon(data = time_seq_alive, aes(x = Time, ymin = log_CFU_lower, ymax = log_CFU_upper), 
              inherit.aes = FALSE, fill = "#19798b", alpha = 0.3) + 
  geom_line(data = time_seq_alive, aes(x = Time, y = log_CFU_fitted, linetype = group, color = group), 
            inherit.aes = FALSE, size = 1, color = "#19798b", , linetype = "dashed") +  
  geom_ribbon(data = time_seq_all, aes(x = Time, ymin = log_CFU_lower, ymax = log_CFU_upper), 
              inherit.aes = FALSE, fill = "#ee9b43", alpha = 0.3) + 
  geom_line(data = time_seq_all, aes(x = Time, y = log_CFU_fitted, linetype = group, color = group), 
            inherit.aes = FALSE, size = 1, color = "#ee9b43") +  
  xlab("Time (h)") +
  ylab("log(Colony Forming Units)") +
  scale_fill_manual(values = c("#19798b", "#ee9b43")) +
  scale_linetype_manual(
    values = c("Logistic fit" = "solid"),
    labels = c("Logistic fit")
  ) +
  scale_color_manual(
    values = c("Logistic fit" = "black"),
    labels = c("Logistic fit")
  ) +
  guides(
    fill = guide_legend(order = 1),  
    linetype = guide_legend(order = 2),
    color = "none"
  ) +  
  theme(legend.position = c(.2, .8)) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32, 36, 48))

ggsave("figures/bacterial_burden_aliveDead.pdf", plot = bact2, width = 5.5, height = 5, units = "in", dpi = 300)

# Pathogen growth vs. health measures

cfuMel <- ggplot(burden_tidy, aes(x = melanization, y = log10(cfu))) +
  geom_jitter(size = 2, shape =21, color = "black", aes(fill = status)) +
  mytheme +
  stat_smooth(method = "gam", color = "grey25") + 
  ylab("log(CFU)") +
  xlab("Melanization")+
  scale_fill_manual(values = c("#19798b", "#ee9b43"))+
  theme(legend.position = "none")

cfuAct <- ggplot(burden_tidy, aes(x = activity, y = log10(cfu))) +
  geom_jitter(size = 2, shape = 21, color = "black", aes( fill = status)) +
  mytheme +
  #stat_smooth(method = "gam", color = "grey25") + 
  ylab("log(CFU)") +
  xlab("Activity")+
  scale_fill_manual(values = c("#19798b", "#ee9b43"))+
  theme(legend.position = "none")

cfuActMel <- ggplot(burden_tidy, aes(x = melanization, y = activity)) +
  geom_jitter(size = 2, shape = 21, color = "black", aes( fill = status)) +
  mytheme +
  #stat_smooth(method = "gam", color = "grey25") + 
  ylab("Activity") +
  xlab("Melanization")+
  scale_fill_manual(values = c("#19798b", "#ee9b43"))+
  theme(legend.position = "none")

small <- plot_grid(cfuAct, cfuMel, ncol = 1, align = "v", labels = c("C", "D"), hjust = -0.5, vjust=1.5)
big <- plot_grid(bact2, small, ncol = 2, align = "v", labels = c("A", ""), rel_widths = c(2, 1.2), hjust = -0.5, vjust=1.5)
ggsave("figures/bacterial_burden_aliveDead.pdf", plot = big, width = 9, height = 4.5, units = "in", dpi = 300)

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
summary(lm(scaled_cfu    ~ scaled_time, data = burden_tidy))
summary(lm(scaled_health ~ scaled_cfu + scaled_time, data = burden_tidy))

# bidirection?
summary(lm(scaled_cfu    ~ scaled_health + scaled_time, data = burden_tidy))

summary(glm(survival     ~ scaled_time, data = burden_tidy, family = binomial(), method = "brglmFit"))
summary(glm(survival     ~ scaled_health, family = binomial, data = burden_tidy, method = "brglmFit"))
summary(glm(survival     ~ scaled_health + scaled_time, data = burden_tidy, family = binomial(), method = "brglmFit"))
summary(glm(survival     ~ scaled_cfu, data = burden_tidy, family = binomial(), method = "brglmFit"))
summary(glm(survival     ~ scaled_cfu + scaled_time, data = burden_tidy, family = binomial(), method = "brglmFit"))
summary(glm(survival     ~ scaled_health + scaled_cfu + scaled_time, family = binomial(), data = burden_tidy, method = "brglmFit"))


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

burden_tidy$time_linear <- burden_tidy$scaled_time
burden_tidy$time_squared <- burden_tidy$scaled_time^2

model_quad <- '
  scaled_cfu ~ t1 * time_linear + t2 * time_squared
  scaled_health ~ a * scaled_cfu + h1 * time_linear 
  survival ~ b * scaled_health + c_prime * scaled_cfu

  # Mediation effect
  indirect := a * b
  total := c_prime + (a * b)
'

#fit_mediation <- sem(mediation_model, data = burden_tidy, ordered = "survival")
#summary(fit_mediation, standardized = TRUE)

#fit_boot <- sem(mediation_model, data = burden_tidy, se = "boot", bootstrap = 1000)
#parameterEstimates(fit_boot, standardized = TRUE) %>% filter(op == ":=")

fit_bayes_sem1 <- bsem(model_quad, data = burden_tidy, burnin = 1000, sample = 5000, n.chains = 3, 
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
  estimate = c(0.89, -0.50, -0.53, 0.49, -0.42, -0.21, 0.01, -0.20),
  lower = c(0.79, -0.59, -0.67, 0.45, -0.56, -0.28, -0.03, -0.27),
  upper = c(0.99, -0.41, -0.39, 0.54, -0.29, -0.14, 0.06, -0.13)
)

effect_plot <- ggplot(effects_df, aes(y = path, x = estimate)) +
  geom_col(fill = "grey", width = 0.6) +
  geom_vline(xintercept = 0, linetype ="dashed")+
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "", x = "Standardized Effect Size") +
  mytheme

ggsave("figures/effect_plot.pdf", plot = effect_plot, width = 6.5, height = 5.5, units = "in", dpi = 300)
