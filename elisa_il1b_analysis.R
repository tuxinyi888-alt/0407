# =============================================================
# ELISA IL-1β Dose-Response Analysis
# Experiments: LPS and BCG (treated as separate experiments)
# =============================================================
# Files (same directory as this script):
#   ELISA_IL1b(LPS).xlsx  —  sheet: "IL-1b"
#   ELISA_IL1b(BCG).xlsx  —  sheet: "AIF006_ELISA_IL1b(BCG)"
#
# Data structure (confirmed from actual files):
#   Donor      : D1–D4
#   Carprofen  : 0, 1, 10, 100   (factor levels used as-is)
#   Stim_Conc  : 0, 1, 10, 100   (numeric; log1p x-scale handles zero)
#   IL1b_pgml  : pg/mL; ~27 NAs in LPS, ~51 NAs in BCG
#
# Output: IL1b_LPS.png/.pdf  |  IL1b_BCG.png/.pdf
# =============================================================

library(readxl)
library(dplyr)
library(ggplot2)

# ── Axis unit labels — edit here if units differ ──────────────
unit_lps  <- "ng/mL"    # LPS concentration unit
unit_bcg  <- "µg/mL"    # BCG concentration unit
unit_carp <- "µg/mL"    # Carprofen dose unit

# =============================================================
# 1. Load data
# =============================================================
# Sheet names differ between files — do NOT change these strings

lps_raw <- read_excel("ELISA_IL1b(LPS).xlsx",
                      sheet = "IL-1b",
                      na    = "NA")

bcg_raw <- read_excel("ELISA_IL1b(BCG).xlsx",
                      sheet = "AIF006_ELISA_IL1b(BCG)",
                      na    = "NA")

# =============================================================
# 2. Pre-process
# =============================================================
# • Apply 1 pg/mL detection floor BEFORE any log transform
# • Carprofen → ordered factor (levels from actual file values)
# • Stim_Conc stays numeric (0, 1, 10, 100); log1p x-scale
#   places zero at the far left then log-spaces 1, 10, 100

prep_data <- function(df) {
  df %>%
    mutate(
      IL1b_plot = case_when(
        is.na(IL1b_pgml) ~ NA_real_,   # preserve NAs
        IL1b_pgml < 1    ~ 1,           # floor sub-detection values
        TRUE             ~ IL1b_pgml
      ),
      Carprofen = factor(Carprofen, levels = c(0, 1, 10, 100)),
      Stim_Conc = as.numeric(Stim_Conc)
    )
}

lps <- prep_data(lps_raw)
bcg <- prep_data(bcg_raw)

# =============================================================
# 3. Summary statistics — mean ± 95% CI per group
# =============================================================
# CI is floored at 1 pg/mL (lower bound cannot go below floor).
# Groups with n = 0 are dropped; CI elements require n ≥ 2.

calc_summary <- function(df) {
  df %>%
    group_by(Carprofen, Stim_Conc) %>%
    summarise(
      n    = sum(!is.na(IL1b_plot)),
      mean = mean(IL1b_plot, na.rm = TRUE),
      sem  = sd(IL1b_plot,   na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) %>%
    filter(n >= 1) %>%
    mutate(
      t_val         = ifelse(n > 1, qt(0.975, df = n - 1), NA_real_),
      ci_lo         = pmax(mean - t_val * sem, 1),    # floor lower CI
      ci_hi         = mean + t_val * sem,
      Stim_Conc_plot = ifelse(Stim_Conc == 0, 0.1, Stim_Conc)  # 0 → 0.1 for log10 x
    )
}

lps_sum <- calc_summary(lps)
bcg_sum <- calc_summary(bcg)

# =============================================================
# 4. Colour palette — purple shades, light → dark
# =============================================================
# Carprofen order: 0 (lightest) → 100 (darkest)

purple_pal <- c(
  "0"   = "#C9A0DC",   # soft lavender
  "1"   = "#9B59B6",   # medium purple
  "10"  = "#6C3483",   # dark purple
  "100" = "#2C0054"    # deep violet
)

# =============================================================
# 5. Consistent y-axis settings (used for BOTH figures)
# =============================================================
# Breaks and limits are identical across LPS and BCG plots.
# Upper limit (200) provides headroom above BCG max (~87 pg/mL).

y_breaks <- c(1, 10, 100)
y_labels <- c("1", "10", "100")
y_limits <- c(0.8, 200)        # 0.8 gives visual padding below the 1 pg/mL floor

# =============================================================
# 6. Plot function
# =============================================================
# Arguments:
#   sum_df      — summary data frame from calc_summary()
#   stim_label  — "LPS" or "BCG" (used in title and x-axis)
#   stim_unit   — concentration unit for the x-axis label
#   x_labels    — character vector of 4 tick labels for breaks c(0.1,1,10,100)
#                 LPS: c("L0","L1","L10","L100")
#                 BCG: c("MOI0","MOI1","MOI10","MOI100")
#
# Plot layers:
#   1. CI ribbon  (pale, only for groups with n ≥ 2)
#   2. Mean line  (connects all groups with n ≥ 1)
#   3. Mean point (open circle, white fill, coloured outline)

make_plot <- function(sum_df, stim_label, stim_unit, x_labels) {

  ggplot() +

    # ── CI ribbon: very pale, matching purple, alpha = 0.18 ──
    geom_ribbon(
      data = filter(sum_df, n >= 2),
      aes(x     = Stim_Conc_plot,
          ymin  = ci_lo,
          ymax  = ci_hi,
          fill  = Carprofen,
          group = Carprofen),
      alpha = 0.18, colour = NA, show.legend = FALSE
    ) +

    # ── Mean line ─────────────────────────────────────────────
    geom_line(
      data = sum_df,
      aes(x      = Stim_Conc_plot,
          y      = mean,
          colour = Carprofen,
          group  = Carprofen),
      linewidth = 0.9
    ) +

    # ── Mean point (open circle, white fill) ──────────────────
    geom_point(
      data = sum_df,
      aes(x      = Stim_Conc_plot,
          y      = mean,
          colour = Carprofen,
          group  = Carprofen),
      size = 3.2, shape = 21, fill = "white", stroke = 1.2
    ) +

    # ── Colour scale ──────────────────────────────────────────
    scale_colour_manual(
      values = purple_pal,
      name   = paste0("Carprofen\n(", unit_carp, ")")
    ) +
    scale_fill_manual(
      values = purple_pal,
      guide  = "none"         # ribbon fill not shown in legend
    ) +

    # ── x-axis: true log10; Stim_Conc 0 is plotted at 0.1 ───────
    # All four tick positions are forced via breaks + labels so
    # missing data at any level never drops its tick mark.
    scale_x_log10(
      limits = c(0.1, 100),
      breaks = c(0.1, 1, 10, 100),
      labels = x_labels,
      expand = expansion(mult = 0.06),
      name   = paste0(stim_label, " Concentration (", stim_unit, ")")
    ) +

    # ── y-axis: log10, consistent breaks across both figures ──
    scale_y_log10(
      breaks = y_breaks,
      labels = y_labels,
      limits = y_limits,
      oob    = scales::squish,   # squish any out-of-range CI tips
      name   = "IL-1\u03b2 (pg/mL)"
    ) +

    # ── Titles ────────────────────────────────────────────────
    labs(
      title   = paste0("IL-1\u03b2 Response to ", stim_label),
      caption = paste0(
        "Mean \u00b1 95% CI; n = 4 donors. ",
        "Values \u003c 1\u2009pg/mL set to detection floor (1\u2009pg/mL)."
      )
    ) +

    # ── Legend: show line + open circle per Carprofen group ───
    guides(
      colour = guide_legend(
        override.aes = list(
          shape     = 21,
          fill      = "white",
          linetype  = 1,
          linewidth = 0.9,
          size      = 3
        )
      )
    ) +

    # ── Publication-style theme ───────────────────────────────
    theme_classic(base_size = 13) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.caption       = element_text(hjust = 0, size = 8, colour = "grey50",
                                        margin = margin(t = 6)),
      axis.title         = element_text(size = 12),
      axis.text          = element_text(size = 11, colour = "black"),
      axis.line          = element_line(colour = "black", linewidth = 0.45),
      axis.ticks         = element_line(colour = "black"),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      legend.title       = element_text(face = "bold", size = 11),
      legend.text        = element_text(size = 10),
      legend.key.height  = unit(1.3, "lines"),
      legend.position    = "right",
      plot.margin        = margin(12, 15, 10, 10)
    )
}

# =============================================================
# 7. Build plots
# =============================================================

p_lps <- make_plot(lps_sum, "LPS", unit_lps,
                   x_labels = c("L0", "L1", "L10", "L100"))
p_bcg <- make_plot(bcg_sum, "BCG", unit_bcg,
                   x_labels = c("MOI0", "MOI1", "MOI10", "MOI100"))

# Preview in RStudio Plots pane
print(p_lps)
print(p_bcg)

# =============================================================
# 8. Save — high-res PNG and vector PDF
# =============================================================

# PNG (300 dpi, suitable for reports and journal submissions)
ggsave("IL1b_LPS.png", plot = p_lps,
       width = 7, height = 5, dpi = 300, units = "in")
ggsave("IL1b_BCG.png", plot = p_bcg,
       width = 7, height = 5, dpi = 300, units = "in")

# PDF (vector, best quality for publications)
ggsave("IL1b_LPS.pdf", plot = p_lps,
       width = 7, height = 5, units = "in")
ggsave("IL1b_BCG.pdf", plot = p_bcg,
       width = 7, height = 5, units = "in")

message("\nSaved:")
message("  IL1b_LPS.png  |  IL1b_LPS.pdf")
message("  IL1b_BCG.png  |  IL1b_BCG.pdf")
