# Fit the three models on Rotterdam and estimate their optimism by bootstrap.
#   1. Cox proportional hazards
#   2. Royston-Parmar flexible parametric model (rstpm2), with time-varying effects
#   3. Random survival forest (randomForestSRC)

suppressPackageStartupMessages({
  library(survival)
  library(splines)
  library(rstpm2)
  library(randomForestSRC)
})

RHS <- "ns(age, df = 3) + meno + size_cat + grade3 + lnodes + lpgr + ler + hormon"
RSF_VARS <- c("age", "meno", "size_cat", "grade3", "nodes", "pgr", "er", "hormon")
F_COX <- as.formula(paste("Surv(time5, status5) ~", RHS))

# ---------- Cox ----------
fit_cox <- function(d) coxph(F_COX, data = d, x = TRUE)

risk_cox <- function(fit, newdata, horizon = HORIZON) {
  # Baseline survival at the mean covariate values, raised to exp(centred linear predictor)
  s0 <- summary(survfit(fit), times = horizon, extend = TRUE)$surv
  1 - s0^exp(predict(fit, newdata = newdata, type = "lp"))
}

# ---------- Royston-Parmar ----------
fit_rp <- function(d, df = 2, tvc = NULL) {
  args <- list(formula = F_COX, data = d, df = df)
  if (length(tvc)) args$tvc <- setNames(as.list(rep(1, length(tvc))), tvc)
  do.call(stpm2, args)
}

risk_rp <- function(fit, newdata, horizon = HORIZON) {
  nd <- newdata
  nd$time5 <- horizon
  1 - as.numeric(predict(fit, newdata = nd, type = "surv"))
}

# Choose baseline df (1-5) by AIC, then add time-varying effects one at a time (forward,
# by AIC) among the predictors whose proportional-hazards test in the Cox model has p < 0.01.
select_rp <- function(d, cox) {
  aics <- sapply(1:5, function(k) AIC(fit_rp(d, df = k)))
  df <- which.min(aics)
  zph <- cox.zph(cox)$table
  cand <- setdiff(rownames(zph)[zph[, "p"] < 0.01], "GLOBAL")
  cand <- intersect(cand, c("lnodes", "lpgr", "ler", "meno", "grade3", "hormon"))  # numeric terms only
  chosen <- character(0)
  best <- AIC(fit_rp(d, df = df))
  path <- data.frame(step = "baseline", tvc = "", aic = best)
  repeat {
    left <- setdiff(cand, chosen)
    if (!length(left)) break
    trial <- sapply(left, function(v) AIC(fit_rp(d, df = df, tvc = c(chosen, v))))
    if (min(trial) >= best - 2) break  # require a real improvement
    v <- names(which.min(trial))
    chosen <- c(chosen, v)
    best <- min(trial)
    path <- rbind(path, data.frame(step = paste("add", v), tvc = paste(chosen, collapse = ", "), aic = best))
  }
  list(df = df, df_aic = aics, tvc = chosen, path = path, zph = zph)
}

# ---------- Random survival forest ----------
fit_rsf <- function(d, nodesize = 15, mtry = 3, ntree = 1000, seed = 2026) {
  rfsrc(Surv(time5, status5) ~ ., data = d[, c("time5", "status5", RSF_VARS)],
        ntree = ntree, nodesize = nodesize, mtry = mtry, seed = -abs(seed), importance = "permute")
}

risk_rsf <- function(fit, newdata, horizon = HORIZON) {
  p <- predict(fit, newdata = newdata[, RSF_VARS])
  j <- max(which(p$time.interest <= horizon))
  1 - p$survival[, j]
}

tune_rsf <- function(d, grid = expand.grid(nodesize = c(5, 15, 30, 60), mtry = c(2, 3, 4))) {
  grid$oob_c <- NA
  for (i in seq_len(nrow(grid))) {
    f <- rfsrc(Surv(time5, status5) ~ ., data = d[, c("time5", "status5", RSF_VARS)],
               ntree = 500, nodesize = grid$nodesize[i], mtry = grid$mtry[i], seed = -2026)
    grid$oob_c[i] <- 1 - tail(f$err.rate[!is.na(f$err.rate)], 1)
  }
  grid[order(-grid$oob_c), ]
}

# ---------- Internal validation (bootstrap optimism, Harrell 1996) ----------
apparent <- function(d, risk) {
  y <- Surv(pmin(d$time, HORIZON), ifelse(d$time > HORIZON, 0, d$status))
  c_h <- concordance(y ~ risk, reverse = TRUE)$concordance
  lp <- cloglog(risk)
  c(c_harrell = unname(c_h), slope = unname(coef(coxph(y ~ lp))))
}

bootstrap_optimism <- function(d, fitter, risker, B = 200, seed = 11) {
  set.seed(seed)
  full <- apparent(d, risker(fitter(d), d))
  opt <- t(replicate(B, {
    b <- d[sample(nrow(d), replace = TRUE), ]
    f <- fitter(b)
    apparent(b, risker(f, b)) - apparent(d, risker(f, d))
  }))
  data.frame(measure = names(full), apparent = full, optimism = colMeans(opt),
             corrected = full - colMeans(opt), row.names = NULL)
}
