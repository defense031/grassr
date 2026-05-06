#!/usr/bin/env Rscript
# run_recovery_figure.R — Clean line-plot version of the q-recovery figure.
# Consumes output/recovery_analysis.rds. Produces a single panel with three
# lines (one per metric), RMSE of q-recovery vs prevalence, aggregated across
# sample sizes. The shape of each line tells the paper's central claim.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ra  <- readRDS("output/recovery_analysis.rds")
agg <- as.data.table(ra$by_prev)

long <- melt(agg,
             id.vars       = "prevalence",
             measure.vars  = c("rmse_kappa", "rmse_PABAK", "rmse_AC1"),
             variable.name = "metric",
             value.name    = "rmse")
long[, metric := factor(metric,
  levels = c("rmse_kappa", "rmse_PABAK", "rmse_AC1"),
  labels = c("kappa", "PABAK", "AC1"))]

p <- ggplot(long, aes(x = prevalence, y = rmse, color = metric, linetype = metric)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2) +
  scale_color_manual(values = c("kappa" = "#E41A1C",
                                "PABAK" = "#377EB8",
                                "AC1"   = "#4DAF4A"),
                     name = NULL) +
  scale_linetype_manual(values = c("kappa" = "solid",
                                   "PABAK" = "solid",
                                   "AC1"   = "solid"),
                        guide = "none") +
  scale_x_continuous(breaks = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = "Prevalence",
    y = expression(paste("RMSE of ", hat(q), " - ", q[true])),
    title = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = c(0.5, 0.95),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 8, 8, 8)
  )

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figures/fig_recovery_lines.png", p,
       width = 7, height = 4, dpi = 300)

cat("Saved output/figures/fig_recovery_lines.png\n")
