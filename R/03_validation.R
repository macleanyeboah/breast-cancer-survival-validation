# Develop on Rotterdam, validate in GBSG, write every table to results/.

suppressPackageStartupMessages(library(jsonlite))

run_analysis <- function(root = ".", B = 200) {
  tab <- function(x, name) write.csv(x, file.path(root, "results/tables", name), row.names = FALSE)
  d <- load_data(root)
  dev <- d$dev
  val <- d$val
  tab(case_mix(dev, val), "01_case_mix.csv")

  # ---- development ----
  cox <- fit_cox(dev)
  sel <- select_rp(dev, cox)
  rp <- fit_rp(dev, sel$df, sel$tvc)
  grid <- tune_rsf(dev)
  rsf <- fit_rsf(dev, grid$nodesize[1], grid$mtry[1])
  tab(data.frame(term = names(coef(cox)), hr = exp(coef(cox)), lower = exp(confint(cox))[, 1],
                 upper = exp(confint(cox))[, 2], p = summary(cox)$coefficients[, "Pr(>|z|)"]),
      "02_cox_hazard_ratios.csv")
  zph <- data.frame(term = rownames(sel$zph), sel$zph, row.names = NULL)
  tab(zph, "03_proportional_hazards_tests.csv")
  tab(sel$path, "04_rp_selection_path.csv")
  tab(grid, "05_rsf_tuning_oob.csv")
  vimp <- sort(rsf$importance, decreasing = TRUE)
  tab(data.frame(variable = names(vimp), permutation_importance = vimp), "06_rsf_importance.csv")

  # ---- internal validation ----
  rp_fitter <- function(x) fit_rp(x, sel$df, sel$tvc)
  internal <- rbind(
    cbind(model = "Cox", bootstrap_optimism(dev, fit_cox, risk_cox, B = B)),
    cbind(model = "Royston-Parmar", bootstrap_optimism(dev, rp_fitter, risk_rp, B = B)),
    data.frame(model = "Random survival forest", measure = "c_harrell",
               apparent = c_index(dev$time, dev$status, risk_rsf(rsf, dev), HORIZON)["harrell"],
               optimism = NA, corrected = 1 - tail(rsf$err.rate[!is.na(rsf$err.rate)], 1))
  )
  rownames(internal) <- NULL
  tab(internal, "07_internal_validation.csv")

  # ---- external validation ----
  risks <- list("Cox" = risk_cox(cox, val), "Royston-Parmar" = risk_rp(rp, val),
                "Random survival forest" = risk_rsf(rsf, val))
  ext <- do.call(rbind, lapply(names(risks), function(m)
    data.frame(model = m, t(performance(val$time, val$status, risks[[m]], HORIZON)))))
  tab(ext, "08_external_validation.csv")

  # Bootstrap CIs for the external C-index differences (resampling GBSG patients)
  set.seed(7)
  diffs <- t(replicate(1000, {
    i <- sample(nrow(val), replace = TRUE)
    cs <- sapply(risks, function(r) c_index(val$time[i], val$status[i], r[i], HORIZON)["harrell"])
    c(rp_minus_cox = cs[[2]] - cs[[1]], rsf_minus_cox = cs[[3]] - cs[[1]])
  }))
  cdiff <- data.frame(comparison = colnames(diffs), estimate = colMeans(diffs),
                      lower = apply(diffs, 2, quantile, 0.025), upper = apply(diffs, 2, quantile, 0.975))
  tab(cdiff, "09_external_c_differences.csv")

  # Calibration curves
  grid_p <- seq(0.10, 0.90, by = 0.01)
  curves <- do.call(rbind, lapply(names(risks), function(m)
    cbind(model = m, calibration(val$time, val$status, risks[[m]], HORIZON, grid = grid_p)$curve)))
  tab(curves, "10_calibration_curves.csv")

  # Recalibration-in-the-large (shifts the risk level, keeps the ranking)
  recal <- do.call(rbind, lapply(names(risks), function(m) {
    rc <- recalibrate_in_the_large(risks[[m]], km_risk(val$time, val$status, HORIZON))
    s <- calibration(val$time, val$status, rc$risk, HORIZON)$stats
    data.frame(model = m, shift = rc$a, ici_before = ext$ici[ext$model == m], ici_after = s[["ici"]],
               oe_after = s[["oe"]])
  }))
  tab(recal, "11_recalibration_in_the_large.csv")

  # Decision curves
  th <- seq(0.05, 0.80, by = 0.01)
  dca <- do.call(rbind, lapply(names(risks), function(m)
    cbind(model = m, net_benefit(val$time, val$status, risks[[m]], HORIZON, th))))
  tab(dca, "12_decision_curves.csv")

  # Sensitivity: develop on node-positive Rotterdam patients only (matches GBSG eligibility)
  dev_np <- dev[dev$nodes > 0, ]
  cox_np <- fit_cox(dev_np)
  rp_np <- fit_rp(dev_np, sel$df, sel$tvc)
  rsf_np <- fit_rsf(dev_np, grid$nodesize[1], grid$mtry[1])
  sens <- do.call(rbind, list(
    data.frame(model = "Cox", t(performance(val$time, val$status, risk_cox(cox_np, val), HORIZON))),
    data.frame(model = "Royston-Parmar", t(performance(val$time, val$status, risk_rp(rp_np, val), HORIZON))),
    data.frame(model = "Random survival forest",
               t(performance(val$time, val$status, risk_rsf(rsf_np, val), HORIZON)))))
  sens$n_development <- nrow(dev_np)
  tab(sens, "13_sensitivity_node_positive_development.csv")

  # Time-varying hazard ratio for PgR from the Royston-Parmar model (figure data)
  q <- quantile(dev$lpgr, c(0.25, 0.75))
  ref <- dev[rep(1, 60), ]
  ref$time5 <- seq(0.25, 5, length.out = 60)
  hi <- lo <- ref
  hi$lpgr <- q[2]
  lo$lpgr <- q[1]
  hr <- predict(rp, newdata = lo, type = "hr", exposed = function(x) { x$lpgr <- q[2]; x }, se.fit = TRUE)
  tab(data.frame(years = ref$time5, hr = hr$Estimate, lower = hr$lower, upper = hr$upper),
      "14_rp_pgr_hazard_ratio_over_time.csv")

  saveRDS(list(dev = dev, val = val, risks = risks), file.path(root, "data/processed/predictions.rds"))

  summary <- list(
    n_dev = nrow(dev), events5_dev = sum(dev$status5), n_val = nrow(val), events5_val = sum(val$status5),
    rp_df = sel$df, rp_tvc = sel$tvc, rp_aic_gain = sel$path$aic[1] - tail(sel$path$aic, 1),
    rsf = list(nodesize = grid$nodesize[1], mtry = grid$mtry[1], ntree = 1000),
    internal = internal, external = ext, c_differences = cdiff, recalibration = recal,
    sensitivity_node_positive = sens,
    versions = c(R = R.version.string, survival = as.character(packageVersion("survival")),
                 rstpm2 = as.character(packageVersion("rstpm2")),
                 randomForestSRC = as.character(packageVersion("randomForestSRC")))
  )
  write_json(summary, file.path(root, "results/summary.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
  invisible(summary)
}
