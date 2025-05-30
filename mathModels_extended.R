library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Parameters for time effects
gamma <- 0.05  # multiplicative time effect
c     <- 0.002 # additive time effect

# New mortality models with time effects
m_lin_exp_add_time <- m_lin_exp + c * t
m_lin_exp_mult_time <- m_lin_exp * exp(gamma * t)

# Add to dataframe
df_models <- df_models %>%
  mutate(
    m_lin_exp_add_time = m_lin_exp_add_time,
    m_lin_exp_mult_time = m_lin_exp_mult_time
  )

# Create mapping data
df_mapped_add_time <- tibble(p = df_models$p_lin, m = df_models$m_lin_exp_add_time)
df_mapped_mult_time <- tibble(p = df_models$p_lin, m = df_models$m_lin_exp_mult_time)

# Panel J: Additive time on m(p)
pI <- ggplot(df_mapped_add_time, aes(x = p, y = m)) +
  geom_line(size = 2, color = "#D55E00") +
  labs(x = "Pathogen density, p", y = "Mortality m(p)\n(additive time)") +
  mytheme +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.length = unit(-.25, "cm")) +
  scale_y_continuous(sec.axis = dup_axis()) +
  scale_x_continuous(sec.axis = dup_axis())

# Panel K: Multiplicative time on m(p)
pJ <- ggplot(df_mapped_mult_time, aes(x = p, y = m)) +
  geom_line(size = 2, color = "#0072B2") +
  labs(x = "Pathogen density, p", y = "Mortality m(p)\n(mult. time)") +
  mytheme +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.length = unit(-.25, "cm")) +
  scale_y_continuous(sec.axis = dup_axis()) +
  scale_x_continuous(sec.axis = dup_axis())

# Updated figure
plot_grid(
  pA,
  pB, pC,
  pD, pE,
  pF, pG,
  pH,
  pI, pJ,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
)
