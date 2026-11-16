# README figures (ggplot2).

suppressPackageStartupMessages(library(ggplot2))

PAL <- c("Cox" = "#2a78d6", "Royston-Parmar" = "#eb6834", "Random survival forest" = "#1baf7a")
theme_clean <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e4e3df"),
          plot.title = element_text(face = "bold", size = 12), plot.title.position = "plot",
          legend.position = "bottom", legend.title = element_blank(),
          plot.background = element_rect(fill = "#fcfcfb", colour = NA))
}
save <- function(p, root, name, w = 8, h = 4.5) ggsave(file.path(root, "results/figures", name), p,
                                                       width = w, height = h, dpi = 160)

make_figures <- function(root = ".") {
  tb <- function(n) read.csv(file.path(root, "results/tables", n), check.names = FALSE)
  pr <- readRDS(file.path(root, "data/processed/predictions.rds"))

  # 1. Kaplan-Meier: the two cohorts
  km <- function(d, name) {
    f <- survfit(Surv(time, status) ~ 1, data = d)
    data.frame(years = c(0, f$time), surv = c(1, f$surv), cohort = name)
  }
  kmd <- rbind(km(pr$dev, sprintf("Rotterdam (development, n = %s)", format(nrow(pr$dev), big.mark = ","))),
               km(pr$val, sprintf("GBSG (validation, n = %d)", nrow(pr$val))))
  p1 <- ggplot(kmd[kmd$years <= 7.5, ], aes(years, surv, colour = cohort)) +
    geom_step(linewidth = 0.9) + geom_vline(xintercept = 5, linetype = 2, colour = "#52514e") +
    scale_colour_manual(values = c("#2a78d6", "#eb6834")) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Years since surgery", y = "Recurrence-free survival",
         title = "GBSG patients relapse sooner; every woman in that trial was node-positive") +
    theme_clean()
  save(p1, root, "fig1_cohorts_km.png")

  # 2. Calibration curves at 5 years
  cc <- tb("10_calibration_curves.csv")
  ext <- tb("08_external_validation.csv")
  lab <- setNames(sprintf("%s (ICI %.3f, slope %.2f)", ext$model, ext$ici, ext$slope), ext$model)
  cc$label <- factor(lab[cc$model], levels = lab)
  risk_df <- do.call(rbind, lapply(names(pr$risks), function(m) data.frame(model = m, risk = pr$risks[[m]])))
  risk_df$label <- factor(lab[risk_df$model], levels = lab)
  p2 <- ggplot(cc, aes(predicted, observed, colour = label)) +
    geom_abline(linetype = 2, colour = "#52514e") +
    geom_line(linewidth = 1.1) +
    geom_rug(data = risk_df, aes(x = risk, colour = label), inherit.aes = FALSE, sides = "b",
             alpha = 0.25, length = unit(0.02, "npc")) +
    scale_colour_manual(values = setNames(PAL[names(lab)], lab)) +
    coord_equal(xlim = c(0.1, 0.9), ylim = c(0.1, 0.9)) +
    scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) +
    labs(x = "Predicted 5-year risk", y = "Observed 5-year risk (smoothed)",
         title = "External calibration in GBSG at 5 years") +
    theme_clean() + guides(colour = guide_legend(ncol = 1))
  save(p2, root, "fig2_calibration_gbsg.png", w = 6.5, h = 6.5)

  # 3. Decision curves
  dca <- tb("12_decision_curves.csv")
  ref <- unique(dca[, c("threshold", "treat_all")])
  p3 <- ggplot(dca, aes(threshold, net_benefit, colour = model)) +
    geom_hline(yintercept = 0, colour = "#52514e") +
    geom_line(data = ref, aes(threshold, treat_all), inherit.aes = FALSE, colour = "#8a8984",
              linetype = 2, linewidth = 0.8) +
    geom_line(linewidth = 1) +
    annotate("text", x = 0.415, y = 0.13, label = "Treat everyone", colour = "#52514e", hjust = 1, size = 3.3) +
    scale_colour_manual(values = PAL) +
    scale_x_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(-0.05, 0.5)) +
    labs(x = "Risk threshold for treatment", y = "Net benefit",
         title = "Decision curves in GBSG: the models beat treating everyone above about a 32% threshold") +
    theme_clean()
  save(p3, root, "fig3_decision_curves.png")

  # 4. Time-varying effect of progesterone receptor (Royston-Parmar)
  hr <- tb("14_rp_pgr_hazard_ratio_over_time.csv")
  p4 <- ggplot(hr, aes(years, hr)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "#52514e") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#eb6834", alpha = 0.18) +
    geom_line(colour = "#eb6834", linewidth = 1.1) +
    labs(x = "Years since surgery", y = "Hazard ratio (PgR upper vs lower quartile)",
         title = "Progesterone receptor protects early, then the effect fades") +
    theme_clean()
  save(p4, root, "fig4_pgr_time_varying_hr.png", w = 7, h = 4.2)
}
