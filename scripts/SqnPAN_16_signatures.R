# SqnPAN_16_signatures.R - A3 mutational signatures (Amendment 4; added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-10-31; tool MutationalPatterns). Descriptive only: no label, no test;
#   no Q2, FoldX, clinical or cross-check label changes. Not blind: Fan et al. 2018 reported the position-3 excess of somatic sequon losses.
#   A prior expectation (clock-like CpG, SBS1) is NOT a result and is not written as one.
# Question: which mutational processes (COSMIC SBS signatures) produce the somatic sequon losses (by sequon position) and gains, compared with background missense?
# Inputs (SqnPAN_14 section 1; md5 read from P2_14_section1_ok.txt): data\processed\P2_14_case_spectra.tsv (coding SNVs of each case's chosen file, 96 channels) and
#   P2_14_a3_mutation_channels.tsv.gz (channel of every chosen-file missense mutation, keyed as SqnPAN_06); data\processed\P2_06_mutations.tsv.gz (md5 3be28f92...).
# Tool: MutationalPatterns (Bioconductor; version printed). fit_to_signatures_strict(mut_matrix, signatures, max_delta = 0.004, method = "backwards") [arguments printed
#   and checked in the code: the script stops if max_delta or method is not an argument]. Signatures: get_known_signatures(muttype = "snv", source = "COSMIC_v3.2",
#   sig_type = "reference", incl_poss_artifacts = FALSE, genome = "GRCh38") [checked in the code; stops if COSMIC_v3.2 is not offered].
# Fits: catalogue per cancer type (all cases) and per cancer type x hypermutated (yes / no; > 1,000 missense rows, SqnPAN_06). A group with fewer than 300 coding SNVs uses
#   its cancer type's all-case fit. Catalogues are exome (CDS) counts, not rescaled to genome trinucleotide composition (stated limit; comparisons inside a group are
#   not affected because every set in a group uses the same exposures).
# Per mutation: probability of each signature = exposure x signature probability of its channel, normalised over signatures (group of its case). Multi-nucleotide,
#   ambiguous, invalid and unjoined rows have no channel and are left out (counted). STOP if fewer than 95% of the T1 loss mutations have a channel.
# Sets (chosen files): T1 loss at N (sequon position 1), T1 loss X to P (position 2), T1 loss at S/T (position 3) (a mutation enters each position class of its lost T1
#   sequons); loss of any tier; gain plausible; gain unknown (gain_state_best, Amendment 2); references: missense at background N (res_role nonseqN), missense to a
#   non-S/T residue at background S/T (nonseqST), all missense.
# Reported: expected share of each signature category (SBS1; SBS5; APOBEC = SBS2 + SBS13; tobacco SBS4; UV SBS7a-d; MMR SBS6, 14, 15, 20, 21, 26, 44; POLE SBS10a-d,
#   28; other) and the model-free share of CpG C>T channels, pooled with 95% CI from 2,000 bootstrap resamples of cases (percentile, seed 32920000), and per cancer type
#   for sets with >= 20 mutations (no CI). Differences with CI: T1 S/T minus T1 N; T1 N minus background N; T1 S/T minus background S/T; loss any, gain plausible and gain
#   unknown minus all missense.
# Controls (section 1, judged): (a) biology: SKCM all-case fit UV share (SBS7a-d) >= 0.50; LUAD and LUSC all-case fits include SBS4 (> 0). (b) fit stability: every fit
#   used is refitted on 20 multinomial redraws of its own reconstructed catalogue at its real SNV count (set.seed(32900000 + repeat)); stable if the cosine similarity
#   of relative exposures is >= 0.90 in >= 16 of 20. If (a) fails, every A3 table is marked _signature_control_FAILED; unstable groups are flagged "fit not stable".
# Seeds: 32900000 + repeat (stability), 32920000 (bootstrap); block 32900000 to 32999999 reserved for A3.
# Usage (R 4.6.0, after SqnPAN_14 has finished; one-time install: install.packages("BiocManager"); BiocManager::install("MutationalPatterns")):
#   source("scripts/SqnPAN_16_signatures.R")
# Logs: results\logs\P2_16_section1_checks_fits_controls.txt, P2_16_section2_results.txt (each ends "# ---- END OF SECTION k ----").
options(stringsAsFactors = FALSE, warn = 1)
ROOT <- Sys.getenv("SQNPAN_ROOT", unset = local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f) == 1) normalizePath(file.path(dirname(normalizePath(f)), ".."), winslash = "/") else normalizePath(getwd(), winslash = "/") }))
LOGS <- file.path(ROOT, "results", "logs"); TABLES <- file.path(ROOT, "results", "tables"); PROC <- file.path(ROOT, "data", "processed")
SPEC_T <- file.path(PROC, "P2_14_case_spectra.tsv"); A3_T <- file.path(PROC, "P2_14_a3_mutation_channels.tsv.gz"); OK14 <- file.path(PROC, "P2_14_section1_ok.txt")
MUTS <- file.path(PROC, "P2_06_mutations.tsv.gz"); MUTS_MD5 <- "3be28f9265312957c13ab2f0faa2b2b2"
VERDICT <- file.path(PROC, "P2_16_controls_verdict.txt"); FITS_RDS <- file.path(PROC, "P2_16_fits.rds")
MIN_SNV <- 300; NREP <- 20; COS_MIN <- 0.90; STAB_N <- 16; NBOOT <- 2000; SEED_STAB <- 32900000; SEED_BOOT <- 32920000; MIN_TYPE <- 20; COVER_MIN <- 0.95
CATS <- list(SBS1 = "SBS1", SBS5 = "SBS5", APOBEC = c("SBS2", "SBS13"), Tobacco = "SBS4", UV = c("SBS7a", "SBS7b", "SBS7c", "SBS7d"),
             MMR = c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44"), POLE = c("SBS10a", "SBS10b", "SBS10c", "SBS10d", "SBS28"))

LG <- new.env()
lg_open <- function(k, name) {
  dir.create(LOGS, showWarnings = FALSE, recursive = TRUE)
  LG$con <- file(file.path(LOGS, sprintf("P2_16_section%d_%s.txt", k, name)), open = "a", encoding = "UTF-8"); LG$k <- k
  lg(strrep("=", 70)); lg(sprintf("SqnPAN_16 section %d (%s) started %s", k, name, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}
lg <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); writeLines(s, LG$con); flush(LG$con) }
lg_close <- function() { lg(sprintf("# ---- END OF SECTION %d ---- %s", LG$k, format(Sys.time(), "%Y-%m-%d %H:%M:%S"))); close(LG$con) }
md5 <- function(p) unname(tools::md5sum(p))
wtsv <- function(df, p) { write.table(df, p, sep = "\t", quote = FALSE, row.names = FALSE); lg("written ", p, " | md5 ", md5(p)) }
cosine <- function(a, b) sum(a * b) / sqrt(sum(a * a) * sum(b * b))

read_inputs <- function() {
  if (!file.exists(OK14)) stop("SqnPAN_14 section 1 has not finished (P2_14_section1_ok.txt missing)")
  ok <- read.delim(OK14, header = FALSE, skip = 1, col.names = c("file", "md5"))
  for (p in c(SPEC_T, A3_T)) {
    rec <- ok$md5[ok$file == basename(p)]; got <- md5(p)
    lg("md5 ", basename(p), ": ", got, if (length(rec) == 1 && got == rec) " (same as SqnPAN_14 recorded)" else " (DIFFERENT or not recorded)")
    if (length(rec) != 1 || got != rec) stop("input md5 differs from the SqnPAN_14 record: stop and decide")
  }
  got <- md5(MUTS); lg("md5 ", basename(MUTS), ": ", got, if (got == MUTS_MD5) " (same as recorded)" else " (DIFFERENT)")
  if (got != MUTS_MD5) stop("SqnPAN_06 mutation table md5 differs")
  sp <- read.delim(SPEC_T, check.names = FALSE, colClasses = c(case_id = "character", project = "character", hypermutated = "character"))
  sp
}

get_sigs <- function() {
  f <- formals(MutationalPatterns::get_known_signatures)
  lg("get_known_signatures arguments: ", paste(names(f), collapse = ", "))
  srcs <- eval(f$source); lg("  source choices: ", paste(srcs, collapse = ", "))
  if (!("COSMIC_v3.2" %in% srcs)) stop("COSMIC_v3.2 is not offered by this MutationalPatterns version: stop and decide")
  a <- list(muttype = "snv", source = "COSMIC_v3.2", sig_type = "reference", incl_poss_artifacts = FALSE, genome = "GRCh38")
  a <- a[names(a) %in% names(f)]; lg("  call: ", paste(names(a), sapply(a, as.character), sep = " = ", collapse = ", "))
  sig <- do.call(MutationalPatterns::get_known_signatures, a)
  lg("  signature matrix ", nrow(sig), " x ", ncol(sig), " | ", paste(colnames(sig), collapse = " "))
  sig
}

fit_one <- function(mm, sig) {
  r <- MutationalPatterns::fit_to_signatures_strict(mm, sig, max_delta = 0.004, method = "backwards")
  ct <- r$fit_res$contribution
  out <- matrix(0, ncol(sig), ncol(mm), dimnames = list(colnames(sig), colnames(mm)))
  out[rownames(ct), colnames(mm)] <- ct[, colnames(mm), drop = FALSE]
  out
}

section1 <- function() {
  sp <- read_inputs()
  suppressPackageStartupMessages(library(MutationalPatterns))
  lg("R ", R.version.string, " | MutationalPatterns ", as.character(packageVersion("MutationalPatterns")))
  f <- formals(fit_to_signatures_strict); lg("fit_to_signatures_strict arguments: ", paste(names(f), sapply(f, function(x) paste(deparse(x), collapse = "")), sep = " = ", collapse = "; "))
  if (!all(c("max_delta", "method") %in% names(f))) stop("fit_to_signatures_strict has no max_delta or method argument: stop and decide")
  lg("fit_to_signatures_strict body:"); for (s in deparse(body(fit_to_signatures_strict))) lg("  ", s)
  sig <- get_sigs()
  ch <- colnames(sp)[5:100]
  if (is.null(rownames(sig)) && nrow(sig) == 96) {
    rownames(sig) <- paste0(rep(rep(c("A", "C", "G", "T"), each = 4), 6), "[", rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each = 16), "]", rep(c("A", "C", "G", "T"), 24))
    lg("WARNING: the signature matrix has no row names; the MutationalPatterns 96-channel order (substitution, then 5' base, then 3' base) was assigned [check]")
  }
  if (is.null(rownames(sig)) || !setequal(rownames(sig), ch)) stop("signature row names do not equal the 96 channel names of the SqnPAN_14 table")
  sig <- sig[ch, , drop = FALSE]
  miss <- setdiff(unlist(CATS), colnames(sig)); lg("category signatures not in the matrix (left out of their category): ", if (length(miss)) paste(miss, collapse = " ") else "none")
  X <- as.matrix(sp[, ch]); proj <- sp$project; hyp <- ifelse(sp$hypermutated == "1", "hyper", "nonhyper")
  allm <- t(rowsum(X, proj)); grp <- paste(proj, hyp, sep = "|"); gm <- t(rowsum(X, grp))
  storage.mode(allm) <- "double"; storage.mode(gm) <- "double"
  ng <- colSums(gm); use_own <- ng >= MIN_SNV
  lg("catalogues: cancer types ", ncol(allm), " | groups (type x hypermutated) ", ncol(gm), " | groups with >= ", MIN_SNV, " SNVs (own fit) ", sum(use_own))
  for (g in colnames(gm)) lg("  group ", g, ": cases ", sum(grp == g), " | SNVs ", ng[g], if (use_own[g]) " | own fit" else " | uses the cancer type's all-case fit")
  t0 <- Sys.time(); fa <- fit_one(allm[, 1, drop = FALSE], sig); lg("timing: one fit ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2), " s")
  fa <- fit_one(allm, sig); fg <- fit_one(gm[, use_own, drop = FALSE], sig)
  rel <- function(m) sweep(m, 2, pmax(colSums(m), 1e-300), "/")
  catshare <- function(r) sapply(CATS, function(s) colSums(r[intersect(s, rownames(r)), , drop = FALSE]))
  ra <- rel(fa); cs <- catshare(ra); if (is.null(dim(cs))) cs <- matrix(cs, nrow = 1, dimnames = list(colnames(ra), names(CATS)))
  for (p in colnames(fa)) lg("  fit ", p, ": ", paste(sprintf("%s %.3f", colnames(cs), cs[p, ]), collapse = " "), " | signatures used ", sum(fa[, p] > 0))
  uv <- cs["TCGA-SKCM", "UV"]; s4 <- c(LUAD = fa["SBS4", "TCGA-LUAD"], LUSC = fa["SBS4", "TCGA-LUSC"])
  bio_ok <- uv >= 0.5 && all(s4 > 0)
  lg("CONTROL (a) biology: SKCM UV share ", sprintf("%.3f", uv), " (limit >= 0.50) | SBS4 contribution LUAD ", sprintf("%.1f", s4[1]), " LUSC ", sprintf("%.1f", s4[2]),
     " (must be > 0) | ", if (bio_ok) "PASSED" else "FAILED")
  fits <- cbind(fa, fg); colnames(fits) <- c(colnames(fa), colnames(fg)); cats <- cbind(allm, gm[, use_own, drop = FALSE]); colnames(cats) <- colnames(fits)
  recon <- sig %*% fits; nsnv <- colSums(cats); relf <- rel(fits); good <- matrix(0L, ncol(fits), 1, dimnames = list(colnames(fits), "n"))
  t0 <- Sys.time()
  for (rep in seq_len(NREP)) {
    set.seed(SEED_STAB + rep)
    red <- sapply(seq_len(ncol(fits)), function(j) rmultinom(1, nsnv[j], recon[, j] / sum(recon[, j]))[, 1])
    dimnames(red) <- list(rownames(sig), colnames(fits)); storage.mode(red) <- "double"
    rf <- rel(fit_one(red, sig))
    for (j in seq_len(ncol(fits))) good[j, 1] <- good[j, 1] + as.integer(cosine(rf[, j], relf[, j]) >= COS_MIN)
  }
  lg("CONTROL (b) stability: ", NREP, " redraws per fit | ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  stable <- good[, 1] >= STAB_N
  for (j in seq_len(ncol(fits))) lg("  ", colnames(fits)[j], ": cosine >= ", COS_MIN, " in ", good[j, 1], " of ", NREP, if (stable[j]) " | stable" else " | fit not stable")
  ftab <- do.call(rbind, lapply(colnames(fits), function(g) data.frame(fit = g, n_snv = nsnv[g], stable = stable[g], signature = rownames(fits),
                                                                       contribution = round(fits[, g], 4), relative = signif(relf[, g], 5))))
  wtsv(ftab[ftab$contribution > 0, ], file.path(TABLES, "P2_16_fits.tsv"))
  gmap <- setNames(ifelse(use_own, names(use_own), sub("\\|.*$", "", names(use_own))), names(use_own))
  saveRDS(list(sig = sig, fits = fits, gmap = gmap, stable = stable, bio_ok = bio_ok, ch = ch), FITS_RDS)
  writeLines(c(if (bio_ok) "PASSED" else "FAILED", sprintf("biology: SKCM UV %.3f; SBS4 LUAD %.1f LUSC %.1f", uv, s4[1], s4[2]),
               sprintf("fits not stable: %s", if (any(!stable)) paste(names(stable)[!stable], collapse = ", ") else "none"),
               paste("written", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), VERDICT)
  lg("verdict ", VERDICT, ": ", if (bio_ok) "PASSED" else "FAILED (A3 tables marked _signature_control_FAILED)")
}

section2 <- function() {
  if (!file.exists(FITS_RDS) || !file.exists(VERDICT)) stop("section 1 has not finished")
  F <- readRDS(FITS_RDS); sig <- F$sig; fits <- F$fits; gmap <- F$gmap; ch <- F$ch; mark <- if (F$bio_ok) "" else "_signature_control_FAILED"
  sp <- read_inputs()
  t0 <- Sys.time()
  m <- read.delim(gzfile(MUTS), colClasses = "character", quote = "")
  m <- m[m$in_chosen_file == "1", c("case_id", "project", "accession", "position", "ref", "alt", "hypermutated", "res_role", "lost_sequon_N", "lost_tiers", "gain_state_best", "class")]
  a3 <- read.delim(gzfile(A3_T), colClasses = "character", quote = "")
  lg("read: chosen-file mutations ", nrow(m), " | channel table ", nrow(a3), " | ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  key <- function(d) paste(d$case_id, d$accession, d$position, d$ref, d$alt, sep = "|")
  m$channel <- a3$channel[match(key(m), key(a3))]; m$ci <- match(m$channel, ch)
  lg("channel join: ", sum(!is.na(m$channel)), " joined | single channel ", sum(!is.na(m$ci)), " | not joined ", sum(is.na(m$channel)),
     " | ", paste(names(table(m$channel[is.na(m$ci) & !is.na(m$channel)])), table(m$channel[is.na(m$ci) & !is.na(m$channel)]), collapse = ", "))
  m$grp <- unname(gmap[paste(m$project, ifelse(m$hypermutated == "1", "hyper", "nonhyper"), sep = "|")])
  if (any(is.na(m$grp))) stop("a case has no fit group")
  # sets
  n <- nrow(m); S <- list()
  isloss <- m$class %in% c("loss", "loss_and_gain"); rel1 <- rel2 <- rel3 <- rep(FALSE, n)
  for (i in which(isloss & grepl("T1", m$lost_tiers, fixed = TRUE))) {
    ns <- as.integer(strsplit(m$lost_sequon_N[i], ";", fixed = TRUE)[[1]]); ts <- strsplit(m$lost_tiers[i], ";", fixed = TRUE)[[1]]
    r <- as.integer(m$position[i]) - ns[ts == "T1"]; rel1[i] <- any(r == 0); rel2[i] <- any(r == 1); rel3[i] <- any(r == 2)
  }
  S$T1_loss_N <- rel1; S$T1_loss_XP <- rel2; S$T1_loss_ST <- rel3; S$loss_any <- isloss
  S$gain_plausible <- m$gain_state_best == "plausible"; S$gain_unknown <- m$gain_state_best == "unknown"
  S$bg_N <- m$res_role == "nonseqN"; S$bg_ST <- m$res_role == "nonseqST" & !(m$alt %in% c("S", "T")); S$all_missense <- rep(TRUE, n)
  lg("T1 loss counts by sequon position (SqnPAN_09 printed N 1,009, X to P 28, S/T 1,721 as its own count): N ", sum(rel1), " | X to P ", sum(rel2), " | S/T ", sum(rel3))
  cov <- data.frame(set = names(S), n = sapply(S, sum), with_channel = sapply(S, function(s) sum(s & !is.na(m$ci))))
  cov$share <- round(cov$with_channel / pmax(cov$n, 1), 4); for (i in seq_len(nrow(cov))) lg("  set ", cov$set[i], ": ", cov$n[i], " mutations | with a channel ", cov$with_channel[i], " (", cov$share[i], ")")
  if (any(cov$share[cov$set %in% c("T1_loss_N", "T1_loss_ST")] < COVER_MIN)) stop("fewer than 95% of the T1 loss mutations have a channel: stop and decide")
  wtsv(cov, file.path(TABLES, "P2_16_coverage.tsv"))
  # posterior category shares per group x channel
  cn <- c(names(CATS), "other"); csig <- lapply(CATS, function(s) intersect(s, rownames(fits))); other <- setdiff(rownames(fits), unlist(csig))
  CM <- sapply(c(csig, list(other = other)), function(s) as.numeric(rownames(fits) %in% s)); rownames(CM) <- rownames(fits)
  grps <- colnames(fits); post <- array(0, c(length(grps), 96, length(cn)), dimnames = list(grps, ch, cn))
  for (g in grps) { pm <- sweep(sig, 2, fits[, g], "*"); rs <- rowSums(pm); pm <- pm / ifelse(rs > 0, rs, 1); post[g, , ] <- pm %*% CM }
  cpg <- as.numeric(substr(ch, 3, 5) == "C>T" & substr(ch, 7, 7) == "G")
  ok <- !is.na(m$ci); mm <- m[ok, ]; gi <- match(mm$grp, grps)
  P <- matrix(0, nrow(mm), length(cn)); for (k in seq_along(cn)) P[, k] <- post[cbind(gi, mm$ci, k)]
  P <- cbind(P, cpg[mm$ci], 1); colnames(P) <- c(cn, "CpG_CtoT", "n")
  cases <- sort(unique(m$case_id)); ci <- match(mm$case_id, cases)
  blocks <- lapply(names(S), function(s) { sel <- S[[s]][ok]; r <- matrix(0, length(cases), ncol(P))
    if (any(sel)) { rs <- rowsum(P[sel, , drop = FALSE], ci[sel]); r[as.integer(rownames(rs)), ] <- rs }; r })
  X <- do.call(cbind, blocks); nb <- ncol(P); setn <- names(S)
  share_of <- function(tot) { o <- matrix(NA, length(setn), nb - 1, dimnames = list(setn, colnames(P)[-nb]))
    for (s in seq_along(setn)) { v <- tot[((s - 1) * nb + 1):(s * nb)]; o[s, ] <- v[-nb] / v[nb] }; o }
  obs <- share_of(colSums(X))
  DIFFS <- list(c("T1_loss_ST", "T1_loss_N"), c("T1_loss_N", "bg_N"), c("T1_loss_ST", "bg_ST"), c("loss_any", "all_missense"), c("gain_plausible", "all_missense"), c("gain_unknown", "all_missense"))
  set.seed(SEED_BOOT); t0 <- Sys.time(); BS <- array(NA, c(NBOOT, length(setn), nb - 1)); nc <- length(cases)
  for (b in seq_len(NBOOT)) { w <- tabulate(sample.int(nc, nc, replace = TRUE), nc); BS[b, , ] <- share_of(as.numeric(crossprod(w, X))) }
  lg("bootstrap: ", NBOOT, " resamples of ", nc, " cases | ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  q <- function(v, p) { v <- sort(v[is.finite(v)]); if (!length(v)) return(NA); v[max(1, ceiling(p * length(v)))] }
  rows <- list(); drows <- list()
  for (s in seq_along(setn)) for (k in seq_len(nb - 1)) rows[[length(rows) + 1]] <- data.frame(set = setn[s], n = sum(S[[s]] & !is.na(m$ci)), category = colnames(P)[k],
      share = round(obs[s, k], 4), lo95 = round(q(BS[, s, k], 0.025), 4), hi95 = round(q(BS[, s, k], 0.975), 4))
  for (d in DIFFS) { a <- match(d[1], setn); b <- match(d[2], setn)
    for (k in seq_len(nb - 1)) drows[[length(drows) + 1]] <- data.frame(contrast = paste(d[1], "minus", d[2]), category = colnames(P)[k], difference = round(obs[a, k] - obs[b, k], 4),
      lo95 = round(q(BS[, a, k] - BS[, b, k], 0.025), 4), hi95 = round(q(BS[, a, k] - BS[, b, k], 0.975), 4)) }
  sh <- do.call(rbind, rows); di <- do.call(rbind, drows); sh$note <- mark; di$note <- mark
  wtsv(sh, file.path(TABLES, "P2_16_shares_pooled.tsv")); wtsv(di, file.path(TABLES, "P2_16_differences.tsv"))
  trows <- list(); proj <- mm$project
  for (s in setn) for (p in sort(unique(proj))) { sel <- S[[s]][ok] & proj == p
    if (sum(sel) >= MIN_TYPE) { v <- colSums(P[sel, , drop = FALSE]); trows[[length(trows) + 1]] <- data.frame(set = s, project = p, n = sum(sel), t(round(v[-nb] / v[nb], 4)), check.names = FALSE) } }
  bt <- do.call(rbind, trows); bt$note <- mark; wtsv(bt, file.path(TABLES, "P2_16_shares_by_type.tsv"))
  lg("=== POOLED SHARES (expected share of each signature category; descriptive, no label)", mark, " ===")
  for (s in setn) lg("  ", s, " (n ", sum(S[[s]] & !is.na(m$ci)), "): ", paste(sprintf("%s %.3f", colnames(P)[-nb], obs[s, ]), collapse = " "))
  lg("=== DIFFERENCES with 95% CI (bootstrap over cases) ===")
  for (i in seq_len(nrow(di))) if (di$category[i] %in% c("SBS1", "SBS5", "APOBEC", "Tobacco", "UV", "MMR", "POLE", "CpG_CtoT"))
    lg("  ", di$contrast[i], " | ", di$category[i], ": ", sprintf("%.3f (%.3f to %.3f)", di$difference[i], di$lo95[i], di$hi95[i]))
}

run_sec <- function(k, name, fn) {
  lg_open(k, name); ok <- TRUE
  tryCatch(fn(), error = function(e) { ok <<- FALSE; lg("ERROR: ", conditionMessage(e)) })
  lg_close(); ok
}
if (run_sec(1, "checks_fits_controls", section1)) run_sec(2, "results", section2) else cat("SECTION 1 DID NOT FINISH: section 2 not run (see the log)\n")
# ---- END OF SCRIPT ----
