

#Plots

############################# mu ###############################################
library(ggplot2)
library(dplyr)
library(gridExtra)

#dat.new <- readRDS("/Users/alessandrasaracini/Downloads/AdjustResults.rds")

#dattt <- readRDS("/Users/alessandrasaracini/Downloads/AdjustResults.rds")

#d <- readRDS("/Users/alessandrasaracini/Downloads/parDATA_76.rds")
dat.new <- readRDS("/Users/alessandrasaracini/Downloads/AdjustResults3200k_15.rds")

##################### Bias #####################################################

library(ggplot2)
library(dplyr)
library(gridExtra)

wrap_plot_bias <- function(data) {
  ggplot(data) +
    geom_line(aes(x = tau_squared_values, y = tau2_Unadj_bias), color = "#E41A1C", linewidth = 0.6) + #Naive
    geom_line(aes(x = tau_squared_values, y = tau2_DGM_bias), color = "#4DAF4A", linewidth = 0.6) + #DGM
    #geom_line(aes(x = tau_squared_values, y = tau2_LL_bias), color = "#377EB8", linewidth = 0.6) + #w_{0}
    geom_line(aes(x = tau_squared_values, y = tau2_LC_bias), color = "#984EA3", linewidth = 0.6) + #w_{1, gamma=3}
    geom_line(aes(x = tau_squared_values, y = tau2_CL_bias), color = "#FF7F00", linewidth = 0.6) + #w_{2, beta=3}
    geom_line(aes(x = tau_squared_values, y = tau2_LL_bias), color = "#377EB8", linewidth = 0.6) + #w_{0}
    geom_line(aes(x = tau_squared_values, y = tau2_CC1_bias), color = "#F781BF", linewidth = 0.6) + #w_{3, gamma=7, beta=1.5}
    geom_line(aes(x = tau_squared_values, y = tau2_CC2_bias), color = "#A65628", linewidth = 0.6) + #w_{3, gamma=1.5, beta=7}
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = "", 
      y = "", 
      title = bquote("K" == .(data$k_values) ~ "," ~ mu == .(data$mu_values))
    ) +
    ylim(-0.3, 0.1) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}

# Create a separate plot for the legend
legend_plot_bias <- ggplot() +
  geom_line(aes(x = 0, y = 0, color = "#E41A1C")) +
  geom_line(aes(x = 0, y = 0, color = "#4DAF4A")) +
  geom_line(aes(x = 0, y = 0, color = "#377EB8")) +
  geom_line(aes(x = 0, y = 0, color = "#984EA3")) +
  geom_line(aes(x = 0, y = 0, color = "#FF7F00")) +
  geom_line(aes(x = 0, y = 0, color = "#F781BF")) +
  geom_line(aes(x = 0, y = 0, color = "#A65628")) +
  
  scale_color_manual(
    values = c(
      "#E41A1C", 
      "#4DAF4A", 
      "#377EB8", 
      "#984EA3",
      "#FF7F00", 
      "#F781BF", 
      "#A65628"
    ),
    labels = c(
      expression(Naive),
      expression(w[DGM]),
      expression(w[0]),
      expression(w[1] ~ gamma == 3),
      expression(w[2] ~ beta == 3),
      expression(w[3] ~ gamma == 7 ~ "," ~ beta == 1.5),
      expression(w[3] ~ gamma == 1.5 ~ "," ~ beta == 7)
    )
  ) +
  
  
  theme_void() +
  theme(
    legend.text = element_text(size = 6),
    legend.key.size = unit(2, "mm"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(color = guide_legend(title = NULL, label.hjust = 0.5))



#################### Bias setting 1.5 ##########################################

# Filter the data for gamma = 1.5
parDATA_15 <- dat.new %>% filter(gamma == 1.5)

# Generate plots
plots_bias <- parDATA_15 %>%
  group_by(mu_values, k_values) %>%
  do(plot = wrap_plot_bias(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot <- grid.arrange(
  arrangeGrob(grobs = plots_bias, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(tau^2)) - tau^2), rot = 90, gp = grid::gpar(fontsize = 8)),
  bottom = grid::textGrob(expression(tau^2), gp = grid::gpar(fontsize = 8))
)

invisible(final_plot)

###################### Bias setting 0.5 ########################################

parDATA_05 <- dat.new %>% filter(gamma == 0.5)


# Generate plots
plots_bias <- parDATA_05 %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_bias(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot <- grid.arrange(
  arrangeGrob(grobs = plots_bias, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(mu)) - mu), rot=90, gp = grid::gpar(fontsize = 8)), 
  bottom = grid::textGrob(expression(mu), gp = grid::gpar(fontsize = 8)))




invisible(final_plot)



##################### Other performance measures ###############################


######################## MSE ###################################################

#Plots
library(ggplot2)
library(dplyr)
library(gridExtra)

#dat.new <- readRDS("/Users/alessandrasaracini/Downloads/AdjustResults1k.rds")

wrap_plot_mse <- function(data) {
  ggplot(data) +
    geom_line(aes(x = tau_squared_values, y = Unadj_mse), color = "#F8766D", linewidth = 0.6) + #Naive
    geom_line(aes(x = tau_squared_values, y = DGM_mse), color = "#00BA38", linewidth = 0.6) + #DGM
    geom_line(aes(x = tau_squared_values, y = LL_mse), color = "#619CFF", linewidth = 0.6) + #w_{0}
    geom_line(aes(x = tau_squared_values, y = LC_mse), color = "#F564E3", linewidth = 0.6) + #w_{1, gamma=3}
    geom_line(aes(x = tau_squared_values, y = CL_mse), color = "#00C19F", linewidth = 0.6) + #w_{2, beta=3}
    geom_line(aes(x = tau_squared_values, y = CC1_mse), color = "#FF61C3", linewidth = 0.6) + #w_{3, gamma=7, beta=1.5}
    geom_line(aes(x = tau_squared_values, y = CC2_mse), color = "#FFB000", linewidth = 0.6) + #w_{3, gamma=1.5, beta=7}
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = "", 
      y = "", 
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$tau_squared_values / (data$tau_squared_values + 2 / 50), 2)))
    ) +
    ylim(0, 0.4) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}

############################### mse 1.5 ########################################

# Filter the data for gamma = 1.5
parDATA_15 <- dat.new %>% filter(gamma == 1.5)

# Generate plots
plots_mse <- parDATA_15 %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_mse(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot_mse <- grid.arrange(
  arrangeGrob(grobs = plots_mse, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(mu)) - mu), rot = 90, gp = grid::gpar(fontsize = 8)),
  bottom = grid::textGrob(expression(mu), gp = grid::gpar(fontsize = 8))
)

invisible(final_plot_mse)

############################### mse 0.5 ########################################

# Filter the data for gamma = 0.5
parDATA_05 <- dat.new %>% filter(gamma == 0.5)

# Generate plots
plots_mse <- parDATA_05 %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_mse(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot_mse <- grid.arrange(
  arrangeGrob(grobs = plots_mse, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(mu)) - mu), rot = 90, gp = grid::gpar(fontsize = 8)),
  bottom = grid::textGrob(expression(mu), gp = grid::gpar(fontsize = 8))
)

invisible(final_plot_mse)

############################# Coverage #########################################


wrap_plot_cov <- function(data) {
  ggplot(data) +
    geom_line(aes(x = tau_squared_values, y = Unadj_cov), color = "#F8766D", linewidth = 0.6) + #Naive
    geom_line(aes(x = tau_squared_values, y = DGM_cov), color = "#00BA38", linewidth = 0.6) + #DGM
    geom_line(aes(x = tau_squared_values, y = LL_cov), color = "#619CFF", linewidth = 0.6) + #w_{0}
    geom_line(aes(x = tau_squared_values, y = LC_cov), color = "#F564E3", linewidth = 0.6) + #w_{1, gamma=3}
    geom_line(aes(x = tau_squared_values, y = CL_cov), color = "#00C19F", linewidth = 0.6) + #w_{2, beta=3}
    geom_line(aes(x = tau_squared_values, y = CC1_cov), color = "#FF61C3", linewidth = 0.6) + #w_{3, gamma=7, beta=1.5}
    geom_line(aes(x = tau_squared_values, y = CC2_cov), color = "#FFB000", linewidth = 0.6) + #w_{3, gamma=1.5, beta=7}
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = "", 
      y = "", 
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$tau_squared_values / (data$tau_squared_values + 2 / 50), 2)))
    ) +
    ylim(0, 1) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}

################### Coverage 1.5 ###############################################


# Filter the data for gamma = 1.5
parDATA_15 <- dat.new %>% filter(gamma == 1.5)

# Generate plots
plots_cov <- parDATA_15 %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_cov(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot_cov <- grid.arrange(
  arrangeGrob(grobs = plots_cov, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(mu)) - mu), rot = 90, gp = grid::gpar(fontsize = 8)),
  bottom = grid::textGrob(expression(mu), gp = grid::gpar(fontsize = 8))
)

invisible(final_plot_cov)

########################## Coverage 0.5 ########################################

# Filter the data for gamma = 0.5
parDATA_05 <- dat.new %>% filter(gamma == 0.5)

# Generate plots
plots_cov <- parDATA_05 %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_cov(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot_cov <- grid.arrange(
  arrangeGrob(grobs = plots_cov, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(mu)) - mu), rot = 90, gp = grid::gpar(fontsize = 8)),
  bottom = grid::textGrob(expression(mu), gp = grid::gpar(fontsize = 8))
)

invisible(final_plot_cov)

############################# Power ############################################

wrap_plot_pow <- function(data) {
  ggplot(data) +
    geom_line(aes(x = tau_squared_values, y = 1-Unadj_power), color = "#F8766D", linewidth = 0.6) + #Naive
    geom_line(aes(x = tau_squared_values, y = 1-DGM_power), color = "#00BA38", linewidth = 0.6) + #DGM
    geom_line(aes(x = tau_squared_values, y = 1-LL_power), color = "#619CFF", linewidth = 0.6) + #w_{0}
    geom_line(aes(x = tau_squared_values, y = 1-LC_power), color = "#F564E3", linewidth = 0.6) + #w_{1, gamma=3}
    geom_line(aes(x = tau_squared_values, y = 1-CL_power), color = "#00C19F", linewidth = 0.6) + #w_{2, beta=3}
    geom_line(aes(x = tau_squared_values, y = 1-CC1_power), color = "#FF61C3", linewidth = 0.6) + #w_{3, gamma=7, beta=1.5}
    geom_line(aes(x = tau_squared_values, y = 1-CC2_power), color = "#FFB000", linewidth = 0.6) + #w_{3, gamma=1.5, beta=7}
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = "", 
      y = "", 
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$tau_squared_values / (data$tau_squared_values + 2 / 50), 2)))
    ) +
    ylim(0, 1) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}


########################### Power 1.5 ##########################################


# Filter the data for gamma = 1.5
parDATA_15 <- dat.new %>% filter(gamma == 1.5)

# Generate plots
plots_pow <- parDATA_15 %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_pow(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot_pow <- grid.arrange(
  arrangeGrob(grobs = plots_pow, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(mu)) - mu), rot = 90, gp = grid::gpar(fontsize = 8)),
  bottom = grid::textGrob(expression(mu), gp = grid::gpar(fontsize = 8))
)

invisible(final_plot_pow)


####################### Power 0.5 ##############################################

# Filter the data for gamma = 0.5
parDATA_05 <- dat.new %>% filter(gamma == 0.5)

# Generate plots
plots_pow <- parDATA_05 %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_pow(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot_pow <- grid.arrange(
  arrangeGrob(grobs = plots_pow, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(mu)) - mu), rot = 90, gp = grid::gpar(fontsize = 8)),
  bottom = grid::textGrob(expression(mu), gp = grid::gpar(fontsize = 8))
)

invisible(final_plot_pow)

######################## Emirical SE ###########################################

wrap_plot_empSE <- function(data) {
  ggplot(data) +
    geom_line(aes(x = tau_squared_values, y = 1-Unadj_emp_se), color = "#F8766D", linewidth = 0.6) + #Naive
    geom_line(aes(x = tau_squared_values, y = 1-DGM_emp_se), color = "#00BA38", linewidth = 0.6) + #DGM
    geom_line(aes(x = tau_squared_values, y = 1-LL_emp_se), color = "#619CFF", linewidth = 0.6) + #w_{0}
    geom_line(aes(x = tau_squared_values, y = 1-LC_emp_se), color = "#F564E3", linewidth = 0.6) + #w_{1, gamma=3}
    geom_line(aes(x = tau_squared_values, y = 1-CL_emp_se), color = "#00C19F", linewidth = 0.6) + #w_{2, beta=3}
    geom_line(aes(x = tau_squared_values, y = 1-CC1_emp_se), color = "#FF61C3", linewidth = 0.6) + #w_{3, gamma=7, beta=1.5}
    geom_line(aes(x = tau_squared_values, y = 1-CC2_emp_se), color = "#FFB000", linewidth = 0.6) + #w_{3, gamma=1.5, beta=7}
    #geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = "", 
      y = "", 
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$tau_squared_values / (data$tau_squared_values + 2 / 50), 2)))
    ) +
    ylim(0.5, 1) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}

######################### gamma = 1.5 ##########################################


# Filter the data for gamma = 1.5
parDATA_15 <- dat.new %>% filter(gamma == 1.5)

# Generate plots
plots_SE <- parDATA_15 %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_empSE(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot_SE <- grid.arrange(
  arrangeGrob(grobs = plots_SE, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(mu)) - mu), rot = 90, gp = grid::gpar(fontsize = 8)),
  bottom = grid::textGrob(expression(mu), gp = grid::gpar(fontsize = 8))
)

invisible(final_plot_SE)

##################### gamma = 0.5 ##############################################


# Filter the data for gamma = 1.5
parDATA_05 <- dat.new %>% filter(gamma == 0.5)

# Generate plots
plots_SE <- parDATA_05 %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_empSE(.)) %>%
  pull(plot)

# Arrange the plots and legend
final_plot_SE <- grid.arrange(
  arrangeGrob(grobs = plots_SE, nrow = 5, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6),
  left = grid::textGrob(expression(E(hat(mu)) - mu), rot = 90, gp = grid::gpar(fontsize = 8)),
  bottom = grid::textGrob(expression(mu), gp = grid::gpar(fontsize = 8))
)

invisible(final_plot_SE)

