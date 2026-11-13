# SqnPAN_21_activedriver.R - B4 ActiveDriver part (Amendment 4; added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-11-07). Post hoc: no earlier label changes.
# Question (independent of the context model of SqnPAN_20): in which proteins are somatic missense mutations enriched at and around T1 sequon N, compared with the rest of
#   the protein, with ActiveDriver's adjustment for disordered regions? ActiveDriver (Reimand and Bader 2013; CRAN 1.0.0) has no trinucleotide-context model.
# Inputs (written by SqnPAN_20 section 4; md5 compared with the SqnPAN_20 section 4 log): proteins with AlphaFold status "use", one fragment, at least one T1 sequon and at least
#   one chosen-file missense mutation (main data, hypermutated cases included); sequences UniProt 2026_03; disorder = 1 where AlphaFold pLDDT < 50, else 0; mutations =
#   chosen-file missense rows (gene = UniProt accession, sample = GDC case); sites = T1 sequon N (residue N, kinase empty).
# Call: ActiveDriver defaults (flank 7, mid_flank 2, regression_type "poisson", enriched_only TRUE, simplified FALSE, skip_mismatch TRUE), mc.cores 1 (Windows). The defaults
#   are checked in the code; the function code is saved to a file (rule 4 of Amendment 4). Proteins are run in chunks (the 200 longest first, timed, then chunks of 500);
#   set.seed before each call. Gene p-values come from the element all_gene_based_fdr [check: printed on the package example first; the script stops if it is missing].
# Label: BH over ALL proteins of the set (a protein without a returned p gets p = 1). "Enriched glycosite region" if FDR < 0.05, else "Not enriched".
# Controls (run and judged before the real run; verdicts written first):
#   Positive control: sites = 8 known driver hotspot residues (CTNNB1 S33, S37, T41, S45; U2AF1 S34; FGFR3 S249; ERBB2 S310; PIK3CA N345) in their 5 proteins; BH over the
#     5 proteins; PASSED if at least 4 of 5 have FDR < 0.05. If FAILED, every label gets the suffix _positive_control_FAILED.
#   Null: sites moved to random N outside any sequon in the same protein (same number as its T1 sites), 3 repeats on 1,000 random proteins (draws made in SqnPAN_20, seeds
#     33340000 + 0 to 3); BH over the 1,000 in each repeat; PASSED if the mean share with FDR < 0.05 is <= 0.05. If FAILED, every label gets the suffix _null_control_FAILED.
# Descriptive: Bailey driver overlap (alias step of SqnPAN_12); agreement with the SqnPAN_20 T1 hotspots (proteins with >= 1 main Hotspot against Enriched proteins).
# Seeds (block 33300000 to 33399999): null runs 33341000 + 10 x repeat + chunk; real runs 33360000 + chunk; positive control 33380000.
# Usage (RStudio or R console):  source("scripts/SqnPAN_21_activedriver.R")
#   Needs install.packages("ActiveDriver") once (and readxl, already used by SqnPAN_12). A restart continues: finished chunks in data\processed\P2_21_parts are reused.
# Log: results\logs\P2_21_activedriver.txt. Outputs: results\tables\P2_21_activedriver_results.tsv, P2_21_activedriver_regions.tsv, P2_21_controls.tsv,
#   P2_21_concordance.tsv; data\processed\P2_21_controls_verdict.txt, P2_21_ActiveDriver_function_code.txt.
local({
ROOT <- Sys.getenv("SQNPAN_ROOT", unset = local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f) == 1) normalizePath(file.path(dirname(normalizePath(f)), ".."), winslash = "/") else normalizePath(getwd(), winslash = "/") })); PROC <- file.path(ROOT, "data", "processed"); TBL <- file.path(ROOT, "results", "tables"); LOGD <- file.path(ROOT, "results", "logs")
PARTS <- file.path(PROC, "P2_21_parts"); dir.create(PARTS, showWarnings = FALSE)
LOG <- file.path(LOGD, "P2_21_activedriver.txt"); con <- file(LOG, open = "a"); sink(con, split = TRUE); on.exit({ sink(); close(con) }, add = TRUE)
t00 <- Sys.time(); cat(strrep("=", 70), "\nP2_21 ActiveDriver, started", format(t00), "\n")
SEED <- list(null = 33341000, real = 33360000, pos = 33380000)
for (s in c(SEED$null + 10 * 3 + 9, SEED$real + 99, SEED$pos)) if (s < 33300000 || s > 33399999) stop("seed outside the B4 block")
FDR_CUT <- 0.05

# ---------------- 1. package and function check ----------------
if (!requireNamespace("ActiveDriver", quietly = TRUE)) stop("ActiveDriver is not installed: run install.packages(\"ActiveDriver\") once, then source this script again")
cat(R.version.string, "| ActiveDriver", as.character(packageVersion("ActiveDriver")), "| readxl", as.character(packageVersion("readxl")), "\n")
AD <- ActiveDriver::ActiveDriver; fm <- formals(AD); print(args(AD))
want <- list(flank = 7, mid_flank = 2, mc.cores = 1, simplified = FALSE, skip_mismatch = TRUE, regression_type = "poisson", enriched_only = TRUE)
for (k in names(want)) if (!identical(eval(fm[[k]]), want[[k]])) stop(paste("default of", k, "differs from the pre-specification:", deparse(fm[[k]])))
cat("defaults checked: flank 7, mid_flank 2, mc.cores 1, simplified FALSE, skip_mismatch TRUE, regression_type poisson, enriched_only TRUE\n")
FCODE <- file.path(PROC, "P2_21_ActiveDriver_function_code.txt")
writeLines(c(deparse(AD), "", "# internal functions of the package:", unlist(lapply(ls(asNamespace("ActiveDriver")), function(n) c(paste0("## ", n), deparse(get(n, asNamespace("ActiveDriver"))))))), FCODE)
cat("function code saved:", FCODE, "md5", unname(tools::md5sum(FCODE)), "\n")
ex <- new.env(); data("ActiveDriver_data", package = "ActiveDriver", envir = ex)
cat("example objects:", paste(ls(ex), collapse = ", "), "\n")
set.seed(SEED$pos - 1)
tm <- system.time(er <- AD(ex$sequences, ex$sequence_disorder, ex$mutations, ex$phosphosites))
cat("example run", round(tm[["elapsed"]], 2), "s; output elements:", paste(names(er), collapse = ", "), "\n")
for (n in names(er)) { x <- er[[n]]; cat(" -", n, "| class", class(x)[1], "| dim", paste(dim(x), collapse = " x "), "| columns", paste(colnames(x), collapse = ", "), "\n"); if (is.data.frame(x)) print(utils::head(x, 3)) }
if (!("all_gene_based_fdr" %in% names(er)) || !all(c("gene", "p") %in% colnames(er$all_gene_based_fdr)))
  stop("all_gene_based_fdr with columns gene and p is not in the output: see the lines above")

# ---------------- 2. inputs ----------------
l4 <- readLines(file.path(LOGD, "P2_20_section4_activedriver_inputs.txt"), warn = FALSE)
rec_md5 <- function(key) { m <- regmatches(l4, regexec(paste0("^ActiveDriver input ", key, ": .* md5 ([0-9a-f]{32})$"), l4)); m <- Filter(length, m); if (!length(m)) stop(paste("no md5 line for", key)); m[[length(m)]][2] }
F <- c(sequences = "P2_20_ad_sequences.tsv", disorder = "P2_20_ad_disorder.tsv", mutations = "P2_20_ad_mutations.tsv.gz", sites_real = "P2_20_ad_sites_real.tsv",
       sites_null = "P2_20_ad_sites_null.tsv", sites_pos = "P2_20_ad_sites_pos.tsv", proteins = "P2_20_ad_proteins.tsv")
for (k in names(F)) { got <- unname(tools::md5sum(file.path(PROC, F[[k]]))); r <- rec_md5(k); cat("md5", F[[k]], got, if (got == r) "(same as SqnPAN_20 log)" else paste("DIFFERENT, SqnPAN_20 log", r), "\n"); if (got != r) stop("input md5 differs") }
rd <- function(k, ...) utils::read.delim(if (grepl("\\.gz$", F[[k]])) gzfile(file.path(PROC, F[[k]])) else file.path(PROC, F[[k]]), stringsAsFactors = FALSE, ...)
sq <- rd("sequences", colClasses = "character"); ds <- rd("disorder", colClasses = "character")
SEQ <- setNames(sq$sequence, sq$gene); DIS <- setNames(ds$disorder, ds$gene)
MUT <- rd("mutations", colClasses = c("character", "character", "integer", "character", "character"))
PRO <- rd("proteins", colClasses = c("character", "character", "integer", "integer", "integer", "character"))
rds <- function(k) { x <- rd(k, colClasses = "character"); x$position <- as.integer(x$position); x$kinase <- ""; x }
S_REAL <- rds("sites_real"); S_NULL <- rds("sites_null"); S_NULL$rep <- as.integer(S_NULL$rep); S_POS <- rds("sites_pos")
stopifnot(!anyDuplicated(names(SEQ)), identical(names(SEQ), names(DIS)), all(nchar(SEQ) == nchar(DIS)), all(grepl("^[01]+$", DIS)))
mm <- sum(substr(SEQ[MUT$gene], MUT$position, MUT$position) != MUT$wt_residue)
sm <- sum(substr(SEQ[S_REAL$gene], S_REAL$position, S_REAL$position) != "N")
cat("proteins", length(SEQ), "| mutations", nrow(MUT), "| reference residue differs", mm, "(must be 0) | real sites", nrow(S_REAL), "in", length(unique(S_REAL$gene)),
    "proteins | site residue not N", sm, "(must be 0) | disordered share", round(sum(nchar(gsub("0", "", DIS))) / sum(nchar(DIS)), 4), "\n")
if (mm || sm) stop("residue checks failed")
REAL <- PRO$gene[grepl("real", PRO$set)]; NULLP <- sort(unique(S_NULL$gene)); POSP <- sort(unique(S_POS$gene))
cat("real proteins", length(REAL), "| null proteins", length(NULLP), "| positive-control proteins", paste(POSP, collapse = ", "), "\n")

run_ad <- function(prots, sites, seed, tag) {
  prots <- intersect(prots, unique(sites$gene)); m <- MUT[MUT$gene %in% prots, ]; s <- sites[sites$gene %in% prots, c("gene", "position", "residue", "kinase")]
  set.seed(seed); t0 <- Sys.time()
  r <- AD(SEQ[prots], DIS[prots], m, s)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  g <- if (!is.null(r) && !is.null(r$all_gene_based_fdr) && nrow(r$all_gene_based_fdr)) r$all_gene_based_fdr[, c("gene", "p")] else data.frame(gene = character(0), p = numeric(0))
  reg <- if (!is.null(r) && is.data.frame(r$all_region_based_pval) && nrow(r$all_region_based_pval)) cbind(tag = tag, r$all_region_based_pval) else NULL
  out <- data.frame(gene = prots, p = 1, returned = FALSE, stringsAsFactors = FALSE); i <- match(g$gene, out$gene); out$p[i[!is.na(i)]] <- g$p[!is.na(i)]; out$returned[i[!is.na(i)]] <- TRUE
  list(g = out, reg = reg, sec = el, n_mut = nrow(m))
}
chunked <- function(prots, sites, seed0, tag, first = NULL, size = 500, max_chunks = Inf) {
  L <- PRO$length[match(prots, PRO$gene)]; prots <- prots[order(-L, prots)]
  idx <- if (is.null(first)) split(seq_along(prots), ceiling(seq_along(prots) / size)) else c(list(seq_len(min(first, length(prots)))), split(seq_along(prots)[-seq_len(first)], ceiling(seq_along(prots[-seq_len(first)]) / size)))
  G <- list(); R <- list()
  for (ci in seq_along(idx)) {
    if (ci > max_chunks) break
    f1 <- file.path(PARTS, sprintf("%s_chunk%02d_genes.tsv", tag, ci)); f2 <- file.path(PARTS, sprintf("%s_chunk%02d_regions.rds", tag, ci))
    if (file.exists(f1) && file.exists(f2)) { G[[ci]] <- utils::read.delim(f1, stringsAsFactors = FALSE); R[[ci]] <- readRDS(f2); next }
    x <- run_ad(prots[idx[[ci]]], sites, seed0 + ci, paste0(tag, "_", ci))
    utils::write.table(x$g, f1, sep = "\t", quote = FALSE, row.names = FALSE); saveRDS(x$reg, f2); G[[ci]] <- x$g; R[[ci]] <- x$reg
    cat(sprintf("  %s chunk %d: %d proteins, %d mutations, %.0f s, p returned for %d\n", tag, ci, length(idx[[ci]]), x$n_mut, x$sec, sum(x$g$returned)))
    if (ci == 1 && !is.null(first)) {
      tot <- sum(PRO$length[match(REAL, PRO$gene)]) + 3 * sum(PRO$length[match(NULLP, PRO$gene)]); done <- sum(PRO$length[match(prots[idx[[1]]], PRO$gene)])
      cat(sprintf("  TIMING: the %d longest proteins took %.0f s; estimate for all real and null runs (by total length) about %.0f minutes\n", length(idx[[1]]), x$sec, x$sec * tot / done / 60))
    }
  }
  list(g = do.call(rbind, G), reg = do.call(rbind, Filter(Negate(is.null), R)))
}

# ---------------- 3. controls ----------------
V <- list(); ctl <- list()
pc <- run_ad(POSP, S_POS, SEED$pos, "positive")$g; pc$fdr <- p.adjust(pc$p, "BH"); pc$symbol <- PRO$symbol[match(pc$gene, PRO$gene)]
print(pc); npc <- sum(pc$fdr < FDR_CUT); V$positive <- npc >= 4
cat("positive control:", npc, "of", nrow(pc), "proteins FDR < 0.05 (rule >= 4 of 5) ->", if (V$positive) "PASSED" else "FAILED", "\n")
ctl[[1]] <- data.frame(control = "positive", rep = 0, proteins = nrow(pc), fdr_below = npc, share = npc / nrow(pc))
# timing on the real data first (the 200 longest real proteins = the first real chunk; its p-values are not looked at before the verdicts)
invisible(chunked(REAL, S_REAL, SEED$real, "real", first = 200, max_chunks = 1))
sh <- c()
for (rp in 1:3) {
  x <- chunked(NULLP, S_NULL[S_NULL$rep == rp, ], SEED$null + 10 * rp, paste0("null", rp))$g
  x <- x[match(NULLP, x$gene), ]; x$p[is.na(x$p)] <- 1; fd <- p.adjust(x$p, "BH"); sh[rp] <- mean(fd < FDR_CUT)
  cat("null repeat", rp, ": proteins", length(NULLP), "| FDR < 0.05", sum(fd < FDR_CUT), "| share", round(sh[rp], 4), "\n")
  ctl[[1 + rp]] <- data.frame(control = "null", rep = rp, proteins = length(NULLP), fdr_below = sum(fd < FDR_CUT), share = sh[rp])
}
V$null <- mean(sh) <= 0.05; cat("null control: mean share", round(mean(sh), 4), "(rule <= 0.05) ->", if (V$null) "PASSED" else "FAILED", "\n")
VF <- file.path(PROC, "P2_21_controls_verdict.txt"); writeLines(c(paste0("positive\t", if (V$positive) "PASSED" else "FAILED"), paste0("null\t", if (V$null) "PASSED" else "FAILED")), VF)
cat("verdict file", VF, "md5", unname(tools::md5sum(VF)), "\n")
wr <- function(x, nm) { f <- file.path(TBL, nm); utils::write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = ""); cat("wrote", f, nrow(x), "rows, md5", unname(tools::md5sum(f)), "\n") }
wr(do.call(rbind, ctl), "P2_21_controls.tsv")

# ---------------- 4. real result ----------------
RR <- chunked(REAL, S_REAL, SEED$real, "real", first = 200)
g <- RR$g[match(REAL, RR$g$gene), ]; g$gene <- REAL; g$p[is.na(g$p)] <- 1; g$returned[is.na(g$returned)] <- FALSE
g$fdr <- p.adjust(g$p, "BH"); i <- match(g$gene, PRO$gene); g$symbol <- PRO$symbol[i]; g$length <- PRO$length[i]; g$n_T1_sites <- PRO$n_T1_sites[i]; g$n_mutations <- PRO$n_mutations[i]
suf <- paste0(if (!V$null) "_null_control_FAILED" else "", if (!V$positive) "_positive_control_FAILED" else "")
g$label <- paste0(ifelse(g$fdr < FDR_CUT, "Enriched glycosite region", "Not enriched"), suf)
b <- as.data.frame(readxl::read_excel(file.path(ROOT, "data", "raw", "Bailey2018_TableS1_driver_genes.xlsx"), sheet = "Table S1", skip = 3, col_types = "text"))
if (nrow(b) != 739L || length(unique(b$Gene)) != 299L) stop("Bailey Table S1 counts differ from the recorded counts")
ALIAS <- data.frame(old = c("FAM46D", "H3F3A", "H3F3C", "HIST1H1C", "HIST1H1E", "RQCD1", "TCEB1", "WHSC1"), new = c("TENT5D", "H3-3A", "H3-5", "H1-2", "H1-4", "CNOT9", "ELOC", "NSD2"), stringsAsFactors = FALSE)
drv <- unique(ifelse(b$Gene %in% ALIAS$old, ALIAS$new[match(b$Gene, ALIAS$old)], b$Gene)); g$bailey_driver <- g$symbol %in% drv
HOT <- file.path(TBL, "P2_20_T1_hotspots.tsv"); hs <- if (file.exists(HOT)) utils::read.delim(HOT, stringsAsFactors = FALSE) else NULL
if (!is.null(hs)) { hm <- hs[startsWith(hs$label_main, "Hotspot"), ]; g$n_T1_hotspots_P2_20 <- as.integer(table(factor(hm$accession, levels = g$gene))) } else { cat("SqnPAN_20 hotspot table not found: concordance skipped\n"); g$n_T1_hotspots_P2_20 <- NA }
g <- g[order(g$p, g$gene), c("gene", "symbol", "length", "n_T1_sites", "n_mutations", "returned", "p", "fdr", "label", "bailey_driver", "n_T1_hotspots_P2_20")]
wr(g, "P2_21_activedriver_results.tsv")
if (!is.null(RR$reg)) wr(RR$reg, "P2_21_activedriver_regions.tsv")
en <- g$fdr < FDR_CUT
cat("\nREAL: proteins", nrow(g), "| p returned", sum(g$returned), "| Enriched glycosite region (FDR < 0.05)", sum(en), "| of them Bailey drivers", sum(en & g$bailey_driver), "\n")
print(utils::head(g[en, c("symbol", "gene", "n_T1_sites", "n_mutations", "p", "fdr", "bailey_driver", "n_T1_hotspots_P2_20")], 40), row.names = FALSE)
if (!is.null(hs)) {
  h <- g$n_T1_hotspots_P2_20 > 0; tab <- table(ActiveDriver_enriched = en, P2_20_hotspot = h); print(tab)
  wr(as.data.frame(tab), "P2_21_concordance.tsv")
}
cat("SqnPAN_21 finished", format(Sys.time()), "| total", round(as.numeric(difftime(Sys.time(), t00, units = "mins")), 1), "minutes\n")
})
# ---- END OF SCRIPT ----
