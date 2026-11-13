# SqnPAN_12_helpers.R - functions for SqnPAN_12_clinical.R (clinical links, analysis plan section 9 and control C5); functions only, no side effects.
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Sourced by SqnPAN_12_clinical.R after an md5 check.
# Needs R (4.6.0 used) and the package survival.
# Globals used at call time (defined in SqnPAN_12_clinical.R before this file is sourced): COV_MAIN2, EQ_HR2.

fmt2 <- function(x) format(x, big.mark = ",", scientific = FALSE)
wr2 <- function(x, f) { utils::write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); unname(tools::md5sum(f)) }
cmd5_2 <- function(x) { tf <- tempfile(); on.exit(unlink(tf)); utils::write.table(x, tf, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); unname(tools::md5sum(tf)) }
rd2 <- function(f) utils::read.delim(f, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
need2 <- function(df, cols, what) { m <- setdiff(cols, names(df)); if (length(m)) stop(what, ": missing columns: ", paste(m, collapse = ", ")); invisible(TRUE) }
pr2 <- function(df, digits = 3) print(format(df, digits = digits), row.names = FALSE)

# ------------------------------------------------------------------ Cox model, stratified by cancer type ------------------------------------------------------------------
# d: data frame with columns tcol, ecol, type, the exposure columns xs (0/1 numbers or a number such as age10) and the covariates covs.
# One row per exposure term: HR, 95% CI, Wald p, TOST p for HR inside [1/EQ_HR2, EQ_HR2]. Complete cases on the columns used.
cox2 <- function(d, xs, covs, tcol = "time", ecol = "event", zph = FALSE) {
  need2(d, c(tcol, ecol, "type", xs, covs), "cox2 data")
  dd <- d[stats::complete.cases(d[, c(tcol, ecol, "type", xs, covs), drop = FALSE]), , drop = FALSE]
  blank <- function(note) data.frame(term = xs, n = nrow(dd), events = sum(dd[[ecol]]), n_group = NA_integer_, events_group = NA_integer_,
                                     b = NA_real_, se = NA_real_, HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_, p = NA_real_, p_tost = NA_real_,
                                     zph_term = NA_real_, zph_global = NA_real_, note = note, stringsAsFactors = FALSE)
  if (nrow(dd) < 10L || sum(dd[[ecol]]) < 5L) return(blank("too few cases or events"))
  if ("sex" %in% covs && length(unique(dd$sex)) < 2L) covs <- setdiff(covs, "sex")
  if (any(vapply(xs, function(x) length(unique(dd[[x]])) < 2L, TRUE))) return(blank("exposure has one value"))
  f <- stats::as.formula(paste0("Surv(", tcol, ", ", ecol, ") ~ ", paste(c(xs, covs), collapse = " + "), " + strata(type)"))
  fit <- tryCatch(withCallingHandlers(survival::coxph(f, data = dd, ties = "efron"), warning = function(w) invokeRestart("muffleWarning")),
                  error = function(e) NULL)
  if (is.null(fit)) return(blank("coxph error"))
  cf <- summary(fit)$coefficients
  zt <- NULL
  if (zph) zt <- tryCatch(survival::cox.zph(fit)$table, error = function(e) NULL)
  m <- log(EQ_HR2)
  rows <- lapply(xs, function(x) {
    r <- blank("coefficient missing")[1L, , drop = FALSE]; r$term <- x
    if (!(x %in% rownames(cf))) return(r)
    b <- cf[x, "coef"]; se <- cf[x, "se(coef)"]
    if (!is.finite(b) || !is.finite(se) || se > 10) { r$note <- "not estimable"; return(r) }
    if (all(dd[[x]] %in% c(0, 1))) { g <- dd[[x]] == 1; r$n_group <- sum(g); r$events_group <- sum(dd[[ecol]][g]) }
    r$b <- b; r$se <- se; r$HR <- exp(b)
    r$HR_lo <- exp(b - stats::qnorm(0.975) * se); r$HR_hi <- exp(b + stats::qnorm(0.975) * se)
    r$p <- cf[x, "Pr(>|z|)"]
    r$p_tost <- max(stats::pnorm((b + m) / se, lower.tail = FALSE), stats::pnorm((b - m) / se))
    if (!is.null(zt) && x %in% rownames(zt)) { r$zph_term <- unname(zt[x, "p"]); r$zph_global <- unname(zt["GLOBAL", "p"]) }
    r$note <- ""
    r
  })
  do.call(rbind, rows)
}

# ------------------------------------------------------------------ power simulation (Weibull on the real covariates, follow-up and cancer types) ------------------------------------------------------------------
# survreg is an AFT model: log T = mu + sigma * W. A planted hazard ratio exp(beta) on a 0/1 exposure means mu shifted by -sigma * beta * x.
# Baseline: type-specific intercepts, common shape, age, sex, log missense count and the other exposure columns at their fitted effects;
# the planted exposure's own fitted effect is removed and replaced by the planted one. Censoring: Weibull with type-specific intercepts (pooled if that fit fails).
sim_prep2 <- function(cc, xs, plant) {
  cc$type <- droplevels(factor(cc$type)); cc$sex <- droplevels(factor(cc$sex))
  f <- stats::as.formula(paste0("Surv(time, event) ~ type + age + sex + logmis + ", paste(xs, collapse = " + ")))
  bf <- withCallingHandlers(survival::survreg(f, data = cc, dist = "weibull"), warning = function(w) invokeRestart("muffleWarning"))
  lp <- unname(stats::predict(bf, type = "lp"))
  b_pl <- unname(stats::coef(bf)[plant])
  if (!is.finite(b_pl)) stop("sim_prep2: the planted exposure ", plant, " has no fitted coefficient")
  base <- lp - b_pl * cc[[plant]]
  cf <- tryCatch(withCallingHandlers(survival::survreg(Surv(time, 1 - event) ~ type, data = cc, dist = "weibull"), warning = function(w) invokeRestart("muffleWarning")),
                 error = function(e) NULL)
  kind <- "type-specific"
  if (is.null(cf) || any(!is.finite(stats::coef(cf))) || any(abs(stats::coef(cf)) > 50)) {
    cf <- withCallingHandlers(survival::survreg(Surv(time, 1 - event) ~ 1, data = cc, dist = "weibull"), warning = function(w) invokeRestart("muffleWarning"))
    kind <- "pooled"
  }
  list(base = base, sigma = bf$scale, muc = unname(stats::predict(cf, type = "lp")), sc = cf$scale, plant = plant, xs = xs, cens_kind = kind,
       hr_fitted = exp(-b_pl / bf$scale), cc = cc)
}
sim_one2 <- function(sf, target, beta, seed) {
  cc <- sf$cc; set.seed(seed); n <- nrow(cc)
  u1 <- stats::runif(n); u2 <- stats::runif(n)
  Tt <- exp(sf$base - sf$sigma * beta * cc[[sf$plant]]) * (-log(u1))^sf$sigma
  Cc <- exp(sf$muc) * (-log(u2))^sf$sc
  d <- cc; d$time <- pmin(Tt, Cc); d$event <- as.integer(Tt <= Cc)
  r <- cox2(d, sf$xs, COV_MAIN2)
  r[r$term == target, , drop = FALSE]
}
# seed = base + 20000 * spec index + 2000 * HR index + repeat (repeat <= 500): all inside base .. base + 99,999
run_power2 <- function(cc, sp, sidx, nrep, hrs, base, sf = NULL) {
  if (is.null(sf)) sf <- sim_prep2(cc, sp$xs, sp$plant)
  res <- lapply(seq_along(hrs), function(h) {
    r <- do.call(rbind, lapply(seq_len(nrep), function(i) sim_one2(sf, sp$target, log(hrs[h]), base + 20000 * sidx + 2000 * h + i)))
    ok <- !is.na(r$p)
    data.frame(spec = sp$name, HR_planted = hrs[h], nrep = nrep, n_failed = sum(!ok),
               power_05 = mean(ok & r$p < 0.05 & r$b > 0), false_pos_two_sided = mean(ok & r$p < 0.05),
               equiv_share = mean(!is.na(r$p_tost) & r$p_tost < 0.05), median_events_group = stats::median(r$events_group, na.rm = TRUE),
               median_HR = stats::median(r$HR, na.rm = TRUE), cens_kind = sf$cens_kind, stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

# ------------------------------------------------------------------ permutation null (real outcome, exposure labels permuted among carriers inside cancer type x hypermutated) ------------------------------------------------------------------
run_null2 <- function(cc3, base, nrep) {
  ix <- split(seq_len(nrow(cc3)), interaction(cc3$type, cc3$hyper, drop = TRUE))
  out <- vector("list", nrep)
  for (r in seq_len(nrep)) {
    set.seed(base + r); lab <- cc3$grp
    for (i in ix) lab[i] <- lab[i[sample.int(length(i))]]
    d <- cc3; d$x_lo <- as.numeric(lab == "loss_only"); d$x_both <- as.numeric(lab == "both")
    z <- cox2(d, c("x_lo", "x_both"), COV_MAIN2); z$rep <- r
    out[[r]] <- z[z$term == "x_lo", , drop = FALSE]
  }
  do.call(rbind, out)
}

# ------------------------------------------------------------------ labels ------------------------------------------------------------------
label2 <- function(p, p_tost, HR) ifelse(is.na(p), "Not_fitted", ifelse(p < 0.05, ifelse(HR > 1, "Associated_worse", "Associated_better"),
                                                                       ifelse(!is.na(p_tost) & p_tost < 0.05, "Equivalent", "Not_clear")))
final2 <- function(lab, sensitive, e1ok, nullok, posok) {
  lab <- ifelse(lab == "Equivalent" & !e1ok, "Equivalent_not_reported", lab)
  lab <- ifelse(lab == "Not_clear" & !sensitive, "Inconclusive_test_not_sensitive", lab)
  lab <- ifelse(grepl("^Associated", lab) & !nullok, paste0(lab, "_null_check_FAILED"), lab)
  if (!posok) lab <- paste0(lab, "_control_FAILED")
  lab
}
