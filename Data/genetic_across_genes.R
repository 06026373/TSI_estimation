library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

df_SGA_gag <- read.csv('gag_windows.csv')
df_SGA_env <- read.csv('env_windows.csv')

# Dataset list
datasets <- list(df_SGA_gag, df_NGS_gag, df_SGA_env, df_NGS_env, df_wide)
names(datasets) <- c("SGA: pdist on gag", "NGS: MAF on gag", "SGA: pdist on env", "NGS: MAF on env", 'NGS: LRTT on env')

# Bootstrap function
boot_mean_ci <- function(x, nboot = 1000) {
  boots <- replicate(nboot, mean(sample(x, replace = TRUE)), simplify = TRUE)
  mean_val <- mean(x)
  ci <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
  data.frame(mean = mean_val, ci_lower = ci[1], ci_upper = ci[2])
}

plot_list <- lapply(seq_along(datasets), function(i) {
  df <- datasets[[i]]
  
  df$TSI <- as.numeric(df$TSI)
  
  df <- df %>%
    mutate(TSI_cat = cut(
      TSI,
      breaks = c(0, 90, 180, 270, 365, 730, 2555),
      labels = c("0-3m", "3-6m", "6-9m", "9-12m", "1-2y", "2-7y"),
      right = TRUE
    ))
  
  df_long <- df %>%
    pivot_longer(cols = starts_with("m"), names_to = "variable", values_to = "value") %>%
    mutate(nombre = as.numeric(sub("m", "", variable)))
  
  results <- df_long %>%
    group_by(TSI_cat, nombre) %>%
    summarise(boot_mean_ci(value), .groups = "drop")
  
  # Change if needed
  xlim <- c(0., 0., 0., 0., 0.01)
  ylim <- c(0.05, 0.02, 0.05, 0.03, 0.16)
  ylabs <- c("pdist", "MAF", "pdist", "MAF", "LRTT")
  
  p <- ggplot(results, aes(x = nombre, y = mean, color = TSI_cat, fill = TSI_cat)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    labs(
      title = names(datasets)[i],
      x = "Sites",
      y = ylabs[i]
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10)) +
    scale_color_viridis_d(name='TSI', na.translate = FALSE) +
    scale_fill_viridis_d(name='TSI', na.translate = FALSE) +
    scale_y_sqrt() +
    coord_cartesian(ylim = c(xlim[i], ylim[i]))
  
  if (i != 2) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "right")
  }
  
  return(p)
})

# Plot
plot_list[[1]] / plot_list[[2]] / plot_list[[3]] / plot_list[[4]] / plot_list[[5]]

