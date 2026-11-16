# Performance measures for a survival prediction model at a fixed horizon.
# Follows McLernon et al. (2023, Ann Intern Med) and Austin, Harrell & van Klaveren (2020, Stat Med).
# Every function takes predicted RISK of the event by the horizon (not survival).

suppressPackageStartupMessages({
  library(survival)
  library(splines)
})

cloglog <- function(p) log(-log(1 - p))

# Inverse-probability-of-censoring weights from the reverse Kaplan-Meier
censoring_km <- function(time, status) {
  fit <- survfit(Surv(time, status == 0) ~ 1)
  function(t) {
    s <- summary(fit, times = t, extend = TRUE)$surv
    pmax(s, 1e-8)
  }
}

# Kaplan-Meier risk of the event by the horizon
km_risk <- function(time, status, horizon) {
  1 - summary(survfit(Surv(time, status) ~ 1), times = horizon, extend = TRUE)$surv
}

# Harrell's C and Uno's C, both truncated at the horizon
c_index <- function(time, status, risk, horizon) {
  y <- Surv(pmin(time, horizon), ifelse(time > horizon, 0, status))
  harrell <- concordance(y ~ risk, reverse = TRUE)
  uno <- concordance(y ~ risk, reverse = TRUE, timewt = "n/G2")
  c(harrell = unname(harrell$concordance), harrell_se = sqrt(harrell$var),
    uno = unname(uno$concordance), uno_se = sqrt(uno$var))
}

# IPCW Brier score and the index of prediction accuracy (IPA = scaled Brier)
brier <- function(time, status, risk, horizon) {
  G <- censoring_km(time, status)
  event_by_h <- time <= horizon & status == 1
  alive_at_h <- time > horizon
  w <- numeric(length(time))
  w[event_by_h] <- 1 / G(time[event_by_h] - 1e-8)
  w[alive_at_h] <- 1 / G(horizon)
  y <- as.numeric(event_by_h)
  b_model <- mean(w * (y - risk)^2)
  b_null <- mean(w * (y - km_risk(time, status, horizon))^2)
  c(brier = b_model, brier_null = b_null, ipa = 1 - b_model / b_null)
}

# Calibration: observed/expected, calibration slope, and a smoothed calibration curve with ICI
calibration <- function(time, status, risk, horizon, grid = NULL) {
  obs <- km_risk(time, status, horizon)
  oe <- obs / mean(risk)
  d <- data.frame(time = pmin(time, horizon), status = ifelse(time > horizon, 0, status),
                  lp = cloglog(risk))
  slope_fit <- coxph(Surv(time, status) ~ lp, data = d)
  # Austin et al. (2020): regress outcome on a spline of cloglog(predicted risk),
  # then read off the "observed" risk at the horizon for every patient.
  smooth <- coxph(Surv(time, status) ~ ns(lp, df = 3), data = d, x = TRUE)
  obs_i <- 1 - summary(survfit(smooth, newdata = d), times = horizon, extend = TRUE)$surv[1, ]
  err <- abs(obs_i - risk)
  curve <- NULL
  if (!is.null(grid)) {
    g <- data.frame(lp = cloglog(grid))
    curve <- data.frame(predicted = grid,
                        observed = 1 - summary(survfit(smooth, newdata = g), times = horizon,
                                               extend = TRUE)$surv[1, ])
  }
  list(stats = c(observed = obs, expected = mean(risk), oe = oe,
                 slope = unname(coef(slope_fit)), slope_se = unname(sqrt(vcov(slope_fit))),
                 ici = mean(err), e50 = unname(median(err)), e90 = unname(quantile(err, 0.9))),
       curve = curve)
}

# Decision curve analysis (net benefit) for survival outcomes at the horizon
net_benefit <- function(time, status, risk, horizon, thresholds) {
  n <- length(time)
  all_risk <- km_risk(time, status, horizon)
  out <- lapply(thresholds, function(pt) {
    pos <- risk >= pt
    nb_model <- if (any(pos)) {
      r_pos <- km_risk(time[pos], status[pos], horizon)
      share <- mean(pos)
      r_pos * share - (1 - r_pos) * share * pt / (1 - pt)
    } else 0
    nb_all <- all_risk - (1 - all_risk) * pt / (1 - pt)
    data.frame(threshold = pt, net_benefit = nb_model, treat_all = nb_all, treat_none = 0)
  })
  do.call(rbind, out)
}

# All measures at once
performance <- function(time, status, risk, horizon) {
  cal <- calibration(time, status, risk, horizon)
  c(c_index(time, status, risk, horizon), brier(time, status, risk, horizon), cal$stats)
}

# Intercept-only recalibration: S_new = S^exp(a), with a chosen so that mean predicted risk
# equals the Kaplan-Meier risk in the new population. Keeps the ranking, shifts the level.
recalibrate_in_the_large <- function(risk, target) {
  f <- function(a) mean(1 - (1 - risk)^exp(a)) - target
  a <- uniroot(f, c(-5, 5))$root
  list(a = a, risk = 1 - (1 - risk)^exp(a))
}
