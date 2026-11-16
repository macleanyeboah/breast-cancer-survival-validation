# Build the development (Rotterdam) and validation (GBSG) datasets with identical predictors.
#
# Outcome: recurrence-free survival (RFS) = time to recurrence or death, whichever comes first,
# following McLernon et al. (2023). Women who died without a recorded recurrence count as an
# event at their death time. Follow-up is cut at 5 years (the prediction horizon).

suppressPackageStartupMessages(library(survival))

HORIZON <- 5  # years

harmonise <- function(d) {
  d$size_cat <- factor(d$size_cat, levels = c("<=20", "20-50", ">50"))
  d$grade3 <- as.integer(d$grade == 3)          # Rotterdam has no grade 1, so 1-2 vs 3
  d$lnodes <- log1p(pmin(d$nodes, 20))          # cap extreme counts, then log
  d$lpgr <- log1p(pmin(d$pgr, 1000))            # receptor values are very skewed
  d$ler <- log1p(pmin(d$er, 1000))
  d$time5 <- pmin(d$time, HORIZON)
  d$status5 <- ifelse(d$time > HORIZON, 0L, d$status)
  d[, c("pid", "time", "status", "time5", "status5", "age", "meno", "size_cat", "grade3",
        "nodes", "lnodes", "pgr", "lpgr", "er", "ler", "hormon")]
}

build_rotterdam <- function(raw) {
  d <- raw
  d$status <- pmax(d$recur, d$death)
  d$time <- ifelse(d$recur == 1, d$rtime, d$dtime) / 365.25
  d$size_cat <- as.character(d$size)
  harmonise(d)
}

build_gbsg <- function(raw) {
  d <- raw
  d$time <- d$rfstime / 365.25
  d$size_cat <- cut(d$size, c(-Inf, 20, 50, Inf), labels = c("<=20", "20-50", ">50"))
  harmonise(d)
}

load_data <- function(root = ".") {
  rott <- build_rotterdam(read.csv(file.path(root, "data/raw/rotterdam.csv")))
  gbsg <- build_gbsg(read.csv(file.path(root, "data/raw/gbsg.csv")))
  list(dev = rott, val = gbsg)
}

case_mix <- function(dev, val) {
  row <- function(name, f) data.frame(characteristic = name, rotterdam = f(dev), gbsg = f(val))
  pct <- function(x) sprintf("%.1f%%", 100 * mean(x))
  med <- function(x) sprintf("%.0f (%.0f-%.0f)", median(x), quantile(x, .25), quantile(x, .75))
  km5 <- function(d) sprintf("%.1f%%", 100 * (1 - summary(survfit(Surv(time, status) ~ 1, d),
                                                         times = HORIZON)$surv))
  rbind(
    row("Patients", function(d) format(nrow(d), big.mark = ",")),
    row("RFS events within 5 years", function(d) format(sum(d$status5), big.mark = ",")),
    row("5-year RFS event risk (Kaplan-Meier)", km5),
    row("Age, median (IQR)", function(d) med(d$age)),
    row("Postmenopausal", function(d) pct(d$meno == 1)),
    row("Tumour > 20 mm", function(d) pct(d$size_cat != "<=20")),
    row("Grade 3", function(d) pct(d$grade3 == 1)),
    row("Node-negative", function(d) pct(d$nodes == 0)),
    row("Positive nodes, median (IQR)", function(d) med(d$nodes)),
    row("PgR, median (IQR) fmol/l", function(d) med(d$pgr)),
    row("ER, median (IQR) fmol/l", function(d) med(d$er)),
    row("Hormone therapy", function(d) pct(d$hormon == 1))
  )
}
