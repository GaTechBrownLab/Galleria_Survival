library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Time vector
t <- seq(0, 50, by = 0.1)

# Parameters
p0 <- 1
r <- 0.1
m0 <- 0.01
beta <- 0.5
a <- 0.001
b <- 0.2
K <- 10
r_log = 0.3
t0 <- 25

# Scenarios
p_exp <- p0 * exp(r * t)
m_exp_lin <- p_exp

p_lin <- p0 + r * t
m_lin_exp <- m0 * exp(beta * p_lin)

p_const <- rep(5, length(t))
m_exp_time <- m0 * exp(b * t)

p_log <- K / (1 + exp(-r_log * (t - t0)))
m_log_time <- a * exp(b * t)
m_log_of_p <- ifelse(p_log >= K, NA, a * ((p_log * (K - p0)) / (p0 * (K - p_log)))^(b / r_log))

# Gompertz survival curve
N_gompertz <- exp(-(a / b) * (exp(b * t) - 1))

# Combine everything into one dataframe
df_models <- tibble(
  t = t,
  p_exp = p_exp,
  m_exp_lin = m_exp_lin,
  p_lin = p_lin,
  m_lin_exp = m_lin_exp,
  p_const = p_const,
  m_exp_time = m_exp_time,
  p_log = p_log,
  m_log_time = m_log_time,
  m_log_of_p = m_log_of_p,
  N_gompertz = N_gompertz
)

# Panel A: Common observable survival curve
pA <- ggplot(df_models, aes(x = t, y = N_gompertz)) +
  geom_line(size = 2, color = "grey50") +
  labs(x = "Time, t", y = "Survival N(t)") +
  mytheme+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank(),  #remove y axis labels
  )+
  theme(axis.ticks.length = unit(-.25, "cm"))+
  scale_y_continuous(sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())

# Panel G: Direct time-dependent m(t)
pB <- ggplot(df_models, aes(x = t, y = m_exp_time)) +
  geom_line(size = 2, color = "#E64B35FF") +
  labs(x = "Time, t", y = "Per-capita mortality m(t)") +
  mytheme+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank(),  #remove y axis labels
  )+
  theme(axis.ticks.length = unit(-.25, "cm"))+
  scale_y_continuous(sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())


# Panel B: Exponential pathogen
pC <- ggplot(df_models, aes(x = t, y = p_exp)) +
  geom_line(size = 2, color = "#91D1C2FF") +
  labs(x = "Time, t", y = "Pathogen density p(t)") +
  mytheme+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank(),  #remove y axis labels
  )+
  theme(axis.ticks.length = unit(-.25, "cm"))+
  scale_y_continuous(sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())

# Panel C: Linear m(p)
df_mapped1 <- tibble(p = df_models$p_exp, m = df_models$m_exp_lin)
pD <- ggplot(df_mapped1, aes(x = p, y = m)) +
  geom_line(size = 2, color = "#F39B7FFF") +
  labs(x = "Pathogen density, p", y = "Mortality m(p)") +
  mytheme+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank(),  #remove y axis labels
  )+
  theme(axis.ticks.length = unit(-.25, "cm"))+
  scale_y_continuous(sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())

# Panel D: Linear pathogen
pE <- ggplot(df_models, aes(x = t, y = p_lin)) +
  geom_line(size = 2, color = "#91D1C2FF") +
  labs(x = "Time, t", y = "Pathogen density p(t)") +
  mytheme+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank(),  #remove y axis labels
  )+
  theme(axis.ticks.length = unit(-.25, "cm"))+
  scale_y_continuous(sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())

# Panel E: Exponential m(p)
df_mapped2 <- tibble(p = df_models$p_lin, m = df_models$m_lin_exp)
pF <- ggplot(df_mapped2, aes(x = p, y = m)) +
  geom_line(size = 2, color = "#F39B7FFF") +
  labs(x = "Pathogen density, p", y = "Mortality m(p)") +
  mytheme+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank(),  #remove y axis labels
  )+
  theme(axis.ticks.length = unit(-.25, "cm"))+
  scale_y_continuous(sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())

# Panel H: Logistic pathogen
pG <- ggplot(df_models, aes(x = t, y = p_log)) +
  geom_line(size = 2, color = "#91D1C2FF") +
  labs(x = "Time, t", y = "Pathogen density p(t)") +
  mytheme+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank(),  #remove y axis labels
  )+
  theme(axis.ticks.length = unit(-.25, "cm"))+
  scale_y_continuous(sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())

# Panel I: Nonlinear threshold-like mapping m(p)
df_mapped4 <- tibble(p = df_models$p_log, m = df_models$m_log_of_p)
pH <- ggplot(df_mapped4, aes(x = p, y = m)) +
  geom_line(size = 2, color = "#F39B7FFF") +
  labs(x = "Pathogen density, p", y = "Mortality m(p)") +
  mytheme+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank(),  #remove y axis labels
  )+
  theme(axis.ticks.length = unit(-.25, "cm"))+
  scale_y_continuous(sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())

# Assemble logical layout
plot_grid(
  pA,
  pB, pC,
  pD, pE,
  pF, pG,
  pH,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H")
)

