###############################################
# IL-1β Figure Generation – LPS & BCG
# Adapted from SandwichELISA-Automated script
###############################################

# Load packages (install if missing)
packages <- c("ggplot2", "dplyr", "cowplot", "tidyr")
installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]
if (length(to_install) > 0) install.packages(to_install)

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

#####
# Set working directory to script location
base_dir <- "."

#####
# Import data
lps_data <- read.csv(file.path(base_dir, "AIF006_ELISA_IL1b(LPS).csv"))
bcg_data <- read.csv(file.path(base_dir, "AIF006_ELISA_IL1b(BCG).csv"))

#####
# Helper: summarise and plot for each stimulus
make_figure <- function(dat, stimulus_name) {

  # Convert Carprofen and Stim_Conc to factors for proper ordering
  dat$Carprofen <- factor(dat$Carprofen, levels = c(0, 1, 10, 100))
  dat$Stim_Conc <- factor(dat$Stim_Conc, levels = c(0, 1, 10, 100))

  # Summary statistics (mean +/- SEM)
  summary_df <- dat %>%
    group_by(Carprofen, Stim_Conc) %>%
    summarise(
      mean_conc = mean(IL1b_pgml, na.rm = TRUE),
      sd_conc   = sd(IL1b_pgml, na.rm = TRUE),
      n         = sum(!is.na(IL1b_pgml)),
      sem_conc  = sd_conc / sqrt(n),
      .groups   = "drop"
    )

  # Plot title
  plot_title <- bquote("IL-1" * beta ~ "(" * .(stimulus_name) * ")")

  # Create figure: grouped bar plot with individual data points
  p <- ggplot() +
    # Bar plot (mean)
    geom_col(data = summary_df,
             aes(x = Stim_Conc, y = mean_conc, fill = Carprofen),
             position = position_dodge(width = 0.8),
             width = 0.7, alpha = 0.7, color = "black", linewidth = 0.3) +
    # Error bars (SEM)
    geom_errorbar(data = summary_df,
                  aes(x = Stim_Conc,
                      ymin = mean_conc - sem_conc,
                      ymax = mean_conc + sem_conc,
                      group = Carprofen),
                  position = position_dodge(width = 0.8),
                  width = 0.25, linewidth = 0.5) +
    # Individual donor points
    geom_point(data = dat %>% filter(!is.na(IL1b_pgml)),
               aes(x = Stim_Conc, y = IL1b_pgml, group = Carprofen),
               position = position_dodge(width = 0.8),
               size = 1.8, shape = 21, fill = "white", color = "black") +
    # Color scale for Carprofen doses
    scale_fill_manual(
      values = c("0" = "#4DAF4A", "1" = "#377EB8", "10" = "#FF7F00", "100" = "#E41A1C"),
      name = expression("Carprofen (" * mu * "M)")
    ) +
    labs(
      title = plot_title,
      x = paste0(stimulus_name, " concentration (ng/mL)"),
      y = expression("IL-1" * beta ~ "(pg/mL)")
    ) +
    theme_minimal(base_size = 16) +
    theme(
      aspect.ratio = 1,
      axis.title = element_text(size = 18),
      axis.text  = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 18),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks.length = unit(2, "mm"),
      axis.ticks = element_line(color = "black"),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 14)
    )

  return(p)
}

#####
# Generate figures
p_lps <- make_figure(lps_data, "LPS")
p_bcg <- make_figure(bcg_data, "BCG")

print(p_lps)
print(p_bcg)

#####
# Combined panel figure
p_combined <- plot_grid(p_lps, p_bcg, labels = c("A", "B"),
                        ncol = 2, label_size = 20, align = "h")
print(p_combined)

#####
# Export figures
fig_width  <- 15
fig_height <- 12.5

# LPS
ggsave(file.path(base_dir, "AIF006_IL1b_LPS_figure.png"),
       plot = p_lps, width = fig_width, height = fig_height, units = "cm", dpi = 300)
ggsave(file.path(base_dir, "AIF006_IL1b_LPS_figure.pdf"),
       plot = p_lps, width = fig_width, height = fig_height, units = "cm")

# BCG
ggsave(file.path(base_dir, "AIF006_IL1b_BCG_figure.png"),
       plot = p_bcg, width = fig_width, height = fig_height, units = "cm", dpi = 300)
ggsave(file.path(base_dir, "AIF006_IL1b_BCG_figure.pdf"),
       plot = p_bcg, width = fig_width, height = fig_height, units = "cm")

# Combined
ggsave(file.path(base_dir, "AIF006_IL1b_combined_figure.png"),
       plot = p_combined, width = 30, height = 14, units = "cm", dpi = 300)
ggsave(file.path(base_dir, "AIF006_IL1b_combined_figure.pdf"),
       plot = p_combined, width = 30, height = 14, units = "cm")
