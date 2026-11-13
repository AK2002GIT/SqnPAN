# SqnPAN_23_dndscv.R - dNdScv at Tier 1 N-glycosite residues (added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# dNdScv at T1 N-glycosite residues. Method: Martincorena I, Raine KM, Gerstung M, Dawson KJ, Haase K, Van Loo P, Davies H, Stratton MR,
# Campbell PJ. Universal patterns of selection in cancer and somatic tissues. Cell 2017;171(5):1029-1041.e21. doi 10.1016/j.cell.2017.09.042,
# PMID 29056346. Package dndscv from GitHub (im3sanger/dndscv), installed 2026-11-09, version 0.0.1.0 (commit printed in section 2).
# Run SqnPAN_23_dndscv_prep.py first. Then (from the repository folder):  Rscript scripts/SqnPAN_23_dndscv.R
# (Rscript scripts/SqnPAN_23_dndscv.R 3 starts at section 3, and so on; finished dNdScv runs are saved in data\processed\P2_23_runs and not run again.)
# Logs: results\logs\P2_23_section2_structure.txt, P2_23_section3_controls.txt, P2_23_section4_result.txt, P2_23_section5_leaveout.txt
# (each ends with "---- END OF SECTION k ----").
#
# ===== PRE-SPECIFICATION (dated 2026-11-09; D1 A, D2 A, D3 B) =====
# The rules for added analyses (analysis plan Amendment 4) apply: added after the primary analysis; no existing label changes; reported beside Q2
# (SqnPAN_09) and A1 (SqnPAN_14), never replacing them. Not blind: Martincorena et al. 2017 reported near-neutral missense dN/dS in tumours, and the
# Q2 and A1 results (Equivalent) were known when this test was designed.
# 1. Question. Do missense mutations at T1 N-glycosite residues (the sequon N and its S/T at +2) occur more or less often than the dNdScv
#    neutral model expects, compared with the other N and S/T missense mutations in the same genes?
# 2. Mutations. One chosen file per case (rule 4A; P2_06_cases.tsv), GRCh38; all SNVs of that file (synonymous included: dNdScv fits its
#    neutral model on them). Indels and DNP/TNP/ONP rows are counted and left out (the site test needs SNVs only).
# 3. Tool. dndscv(refdb = "hg38"): built-in refcds_GRCh38_hg38, covariates switched by the package to hg19_hg38_epigenome_pcawg, sm
#    "192r_3w", outmats = TRUE; every other argument at its default (printed in section 2). D1 A: the primary run uses the dNdScv defaults
#    (max_coding_muts_per_sample 3000, max_muts_per_gene_per_sample 3); an uncapped run (both 1e9) and a run without the hypermutated cases
#    (SqnPAN_06 rule, more than 1,000 missense rows) are reported beside it.
# 4. Gene to UniProt. Strict rule (primary): a RefCDS gene is used if its translation (final stop removed) is identical to the UniProt
#    2026_03 sequence of the accession with the same primary gene name (codon i = residue i). Looser rule (beside): for the other
#    accessions of a RefCDS gene, residue i maps to codon k if the UniProt window i-10..i+10 (clipped at the ends, at least 15 residues)
#    occurs exactly once in the RefCDS translation. Fixed now: if the strict rule covers fewer than 60% of the T1 sequon positions of the
#    SqnPAN_08 list, the looser rule (strict plus window) becomes primary.
# 5. Counting (Q2 rules). At N every missense change counts; at S/T only changes to a residue that is not S or T (S<->T removed from both
#    classes); nonsense and synonymous changes are not counted. Observed = SNVs annotated Missense by dNdScv (annotmuts) at the mapped codon.
#    Expected at a position = relmr(gene) x the sum, over the counted single-base changes of the codon, of the absolute trinucleotide rate,
#    as sitednds computes it (relmr = genemuts exp_syn_cv / exp_syn; rates = mle_submodel x t, with TTT>TGT = 1).
# 6. Primary test (D2 A). Mantel-Haenszel rate ratio, strata = gene x residue type (N; S/T), exposure = expected count, glycosite positions
#    (role = tier) against positions in no sequon (role 0); log2 scale. 95% and 90% CI from 2,000 bootstrap resamples over genes (seed
#    33620000); p from the stratified score test. Labels as in Q2: Enriched or Depleted = BH FDR < 0.05 and |log2| >= 0.32; Equivalent =
#    90% CI inside -0.14..0.14; Inconclusive otherwise. BH family = the 9 pooled units (T1, T2, T3 x N, S/T, both) of a set-up.
#    PRIMARY UNIT: T1 N + S/T, main run, primary mapping. Same seed twice (runs 1 and 2) and seed 33630000 once (run 3); labels from run 1.
# 7. Beside (no new label). Absolute dN/dS at T1 positions and at background positions of the same genes (observed / expected, Poisson 95%
#    CI, divided by the synonymous calibration of C-Cal); the sitednds globaldnds_knownsites output for the T1 list; the dNdScv global
#    dN/dS (globaldnds). Set-ups: main run with the other mapping rule; uncapped run; no-hypermutated run; main run without the 299 Bailey
#    driver genes. D3 B: leave out, one at a time, each of the 5 projects with the most SNVs in the dNdScv input (a new dNdScv run each;
#    primary unit and its N and S/T parts only; descriptive).
# 8. Controls (section 3), run before the real result; none reads observed counts at real glycosite positions.
#    C-A code check: this script's own code gives the same observed and expected sums as sitednds for background, synonymous and hotspot
#        lists, and its mutation keys match at least 99.9% of the dNdScv annotation. FAIL = code error: the script stops; no label.
#    C-P1 drivers: TP53, KRAS and PIK3CA each with qglobal_cv < 0.1 in the main run. Bailey genes detected: descriptive.
#    C-P2 site route: all missense changes at KRAS G12, G13, Q61; NRAS Q61; BRAF V600; PIK3CA E542, E545, H1047; IDH1 R132; TP53 R175, R248,
#        R273 must be labelled Enriched by the primary statistic.
#    C-Cal scale: the synonymous list gives dN/dS between 0.8 and 1.25; gates only the absolute readout of point 7.
#    C-N null: 20 random sets of background positions matched in number per gene and residue type (seeds 33600001 to 33600020). Changed
#        in at most 1 of 20 repeats = pass; Equivalent in at least 16 of 20 = Equivalent allowed for that unit.
#    C-Pl planted: ratios 1.5 and 0.67 labelled Changed in the right direction in at least 16 of 20 repeats (seeds 33610001 to 33610020);
#        1.25 and 0.80 reported, not judged.
#    If C-P1, C-P2 or C-N fails, the result it guards reads "Inconclusive, test not sensitive", with the raw label beside it.
# 9. Reading. Equivalent = a shift of at most about 10% in either direction, never "no selection"; the wording rule for Equivalent results still applies.
# 10. Seeds: block 33600000 to 33699999, checked in the code.
# 11. Files: SqnPAN_23_dndscv_prep.py (section 1) and this script (sections 2 to 5); one command each; tables results\tables\P2_23_*.tsv.
#     Time: a guess, not measured (20 to 60 minutes per dNdScv run; 8 runs); section 3 prints the measured time of the first run.
#
# Implementation details fixed while writing this script, before any result (2026-11-09; within the pre-specified text):
#  a. Strict rule: a gene name used by more than one RefCDS entry, or a RefCDS gene identical to more than one accession, is left out;
#     window rule only for genes with one candidate accession and no strict match. A codon reached by two positions is left out.
#  b. Stop codons: changes at a stop codon are in no class; changes to a stop are nonsense (not counted).
#  c. "Inside" the Equivalent margin includes a bound exactly on +-0.14. A unit is "Not tested (too few data)" if the expected glycosite
#     events under H0 are below 5. Bootstrap resamples the informative genes of the unit (genes with a stratum that has expected counts in
#     both classes and at least one event); CIs need at least 90% finite resamples. p: stratified score test (glycosite events in a stratum
#     ~ binomial(all events, Eg / (Eg + Eb))).
#  d. C-A (a): coding base, trinucleotide, codon number and reference residue agree for >= 99.9% of the Missense and Synonymous annotmuts rows
#     of the genes in the structures. C-A (b): lists = background changes in 200 random genes with T1 positions (seed 33640000), all
#     synonymous changes of the same genes, all missense changes at the hotspot codons; expected sums equal (relative difference < 1e-6)
#     and observed counts differ by at most 0.1%.
#  e. C-P2: strata gene x reference residue; single test (no BH); Enriched = p < 0.05 and log2 >= 0.32. A hotspot codon whose RefCDS residue
#     differs is left out and listed.
#  f. C-N: in each stratum the number of tier-t glycosite positions is drawn from its background positions (a stratum with too few background
#     positions is skipped and counted); BH over the 9 units of a repeat; bootstrap seed 33620000. Capability is judged per unit.
#  g. C-Pl: glycosite events per stratum ~ Poisson(Eg x Ob / Eb x ratio), ratios 1.5, 0.67, 1.25, 0.80 and 1.0; T1 units only, BH over the 3;
#     judged on the primary unit. A C-Pl failure is reported beside the result (the approved text names no consequence for it).
#  h. D3 B: the 5 projects come from the SNV counts of the dNdScv input (section 2 prints them).
#  i. Dated correction 2026-11-09 (found in the first section 3 log: the C-P1 line printed NULL for all three genes and passed with 0 genes
#     detected, because sel_cv has no column qglobal_cv when the input has no indels). C-P1 now uses qglobal_cv if present, else qallsubs_cv (the
#     combined missense + truncating q-value), and stops if any of the three genes has no q value. Rule unchanged: each q < 0.1. The first section 3
#     run was stopped before it finished; its log lines stay in the log file and are not used.
# ========================================================================================================================================

options(scipen = 999, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
START <- if (length(args) >= 1) suppressWarnings(as.integer(args[1])) else 2L
if (is.na(START) || !(START %in% 2:5)) stop("the argument must be 2, 3, 4 or 5")
SCRIPT_DIR <- local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1) dirname(normalizePath(f)) else getwd() })
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "data", "processed"))) ROOT <- Sys.getenv("SQNPAN_ROOT", unset = normalizePath(getwd(), winslash = "/"))
PROC <- file.path(ROOT, "data", "processed"); TABLES <- file.path(ROOT, "results", "tables"); LOGS <- file.path(ROOT, "results", "logs")
RUNS <- file.path(PROC, "P2_23_runs")
for (d in c(TABLES, LOGS, RUNS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
MUT <- file.path(PROC, "P2_23_mutations.tsv.gz"); CASES <- file.path(PROC, "P2_23_cases.tsv"); POSF <- file.path(PROC, "P2_23_positions.tsv.gz")
UNIF <- file.path(PROC, "P2_23_uniprot.tsv"); BAILF <- file.path(PROC, "P2_23_bailey.txt"); GATE1 <- file.path(PROC, "P2_23_section1_ok.txt")
MAPF <- file.path(PROC, "P2_23_map.tsv.gz"); GATE2 <- file.path(PROC, "P2_23_section2_ok.txt"); VERD <- file.path(PROC, "P2_23_controls_verdict.txt")
CAPF <- file.path(TABLES, "P2_23_capability.tsv")
SEED_NULL <- 33600000L + 1:20; SEED_PLANT <- 33610000L + 1:20; SEED_BOOT <- 33620000L; SEED_BOOT3 <- 33630000L; SEED_ANCHOR <- 33640000L
NBOOT <- 2000L; EQ <- 0.14; CH <- 0.32; FDRCUT <- 0.05; INS <- "Inconclusive, test not sensitive"; NT <- "Not tested (too few data)"
B <- c("T", "C", "A", "G"); AA <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
AA_LET <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
COMP <- c(T = "A", C = "G", A = "T", G = "C")
HOT <- data.frame(gene = c("KRAS", "KRAS", "KRAS", "NRAS", "BRAF", "PIK3CA", "PIK3CA", "PIK3CA", "IDH1", "TP53", "TP53", "TP53"),
                  aa = c("G", "G", "Q", "Q", "V", "E", "E", "H", "R", "R", "R", "R"),
                  codon = c(12L, 13L, 61L, 61L, 600L, 542L, 545L, 1047L, 132L, 175L, 248L, 273L))
cidx <- function(b1, b2, b3) 16L * (b1 - 1L) + 4L * (b2 - 1L) + (b3 - 1L) + 1L

# ---------------- logging ----------------
.LOGFILE <- NULL
L <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); if (!is.null(.LOGFILE)) cat(s, "\n", sep = "", file = .LOGFILE, append = TRUE); invisible(NULL) }
LP <- function(x, indent = "    ") for (ln in capture.output(print(x))) L(indent, ln)
run_section <- function(k, name, fn) {
  .LOGFILE <<- file.path(LOGS, sprintf("P2_23_section%d_%s.txt", k, name))
  L(strrep("=", 70)); L(sprintf("SqnPAN_23 section %d (%s) started %s", k, name, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  t0 <- Sys.time(); ok <- FALSE; nwarn <- 0L
  tryCatch(withCallingHandlers({ fn(); ok <- TRUE },
      message = function(m) { L("  [message] ", trimws(conditionMessage(m))); invokeRestart("muffleMessage") },
      warning = function(w) { nwarn <<- nwarn + 1L; L("  [warning] ", conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) { L("ERROR: ", conditionMessage(e)); cl <- conditionCall(e); if (!is.null(cl)) L("  in: ", paste(deparse(cl), collapse = " ")) },
    finally = {
      L(sprintf("warnings in this section: %d", nwarn))
      L(sprintf("minutes: %.1f", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      L(sprintf("---- END OF SECTION %d ---- %s", k, format(Sys.time(), "%Y-%m-%d %H:%M:%S"))) })
  ok
}
md5 <- function(f) unname(tools::md5sum(f))
mins_since <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
read_kv <- function(f) { x <- read.delim(f, header = FALSE, colClasses = "character", quote = ""); setNames(x[[2]], x[[1]]) }
check_gate1 <- function() {
  if (!file.exists(GATE1)) stop("run SqnPAN_23_dndscv_prep.py first (gate file of section 1 missing)")
  g <- read.delim(GATE1, header = FALSE, colClasses = "character")
  for (i in seq_len(nrow(g))) { m <- md5(file.path(PROC, g[i, 1]))
    L(sprintf("  md5 %s %s %s", g[i, 1], m, if (identical(m, g[i, 2])) "(as in section 1)" else "DIFFERS from section 1"))
    if (!identical(m, g[i, 2])) stop("an input changed after section 1") }
}
check_seeds <- function() {
  s <- c(SEED_NULL, SEED_PLANT, SEED_BOOT, SEED_BOOT3, SEED_ANCHOR)
  ok <- all(s >= 33600000L & s <= 33699999L) && !anyDuplicated(s)
  L("seeds: ", length(s), " values, all in the block 33600000 to 33699999 and none repeated: ", ok); if (!ok) stop("seed check failed")
}

# ---------------- RefCDS helpers ----------------
load_refcds <- function() { e <- new.env(); utils::data("refcds_GRCh38_hg38", package = "dndscv", envir = e); e }
seqchar <- function(x) { if (is.null(x)) return(NA_character_); s <- as.character(x); if (length(s) > 1) s <- paste(s, collapse = ""); toupper(s) }
translate_gene <- function(s) {
  if (is.na(s)) return(NA_character_); b <- match(strsplit(s, "")[[1]], B); n <- length(b)
  if (n < 3L || n %% 3L != 0L || anyNA(b)) return(NA_character_)
  m <- matrix(b, nrow = 3L); aa <- AA[cidx(m[1, ], m[2, ], m[3, ])]; if (aa[length(aa)] == "*") aa <- aa[-length(aa)]
  paste(aa, collapse = "")
}
win_map <- function(useq, prot, pos) {
  n <- nchar(useq); a <- pmax(1L, pos - 10L); b <- pmin(n, pos + 10L); out <- rep(NA_integer_, length(pos))
  for (q in seq_along(pos)) { if (b[q] - a[q] + 1L < 15L) next
    m <- gregexpr(substr(useq, a[q], b[q]), prot, fixed = TRUE)[[1]]
    if (length(m) == 1L && m[1] > 0L) out[q] <- as.integer(m[1]) + pos[q] - a[q] }
  out
}
prep_genes <- function(R, gidx) {
  t0 <- Sys.time(); gidx <- gidx[!is.na(gidx)]
  sq <- vapply(gidx, function(j) seqchar(R[[j]]$seq_cds), ""); up <- vapply(gidx, function(j) seqchar(R[[j]]$seq_cds1up), "")
  dn <- vapply(gidx, function(j) seqchar(R[[j]]$seq_cds1down), "")
  ok <- !is.na(sq) & nchar(sq) %% 3L == 0L & nchar(up) == nchar(sq) & nchar(dn) == nchar(sq) & !grepl("[^ACGT]", sq) & !grepl("[^ACGT]", up) & !grepl("[^ACGT]", dn)
  L("  gene structures: ", sum(ok), " genes (", sum(!ok), " left out: length not a multiple of 3, flank length differs or a base other than A, C, G, T)")
  gidx <- gidx[ok]; sq <- sq[ok]; up <- up[ok]; dn <- dn[ok]
  nm <- vapply(gidx, function(j) as.character(R[[j]]$gene_name), "")
  len <- nchar(sq); ncod <- as.integer(len %/% 3L); coff <- c(0L, cumsum(ncod))[seq_along(ncod)]
  ref <- match(strsplit(paste(sq, collapse = ""), "")[[1]], B); upv <- match(strsplit(paste(up, collapse = ""), "")[[1]], B)
  dnv <- match(strsplit(paste(dn, collapse = ""), "")[[1]], B)
  gpos <- vector("list", length(gidx)); minus <- logical(length(gidx)); chr <- character(length(gidx)); bad <- 0L
  for (i in seq_along(gidx)) { x <- R[[gidx[i]]]; iv <- as.matrix(x$intervals_cds); iv <- iv[order(iv[, 1]), , drop = FALSE]
    p <- unlist(mapply(function(a, b) seq.int(a, b), iv[, 1], iv[, 2], SIMPLIFY = FALSE), use.names = FALSE)
    minus[i] <- as.character(x$strand) %in% c("-1", "-")
    if (minus[i]) p <- rev(p)
    if (length(p) != len[i]) bad <- bad + 1L
    gpos[[i]] <- as.integer(p); chr[i] <- as.character(x$chr) }
  if (bad > 0L) stop("CDS interval length differs from the sequence length in ", bad, " genes")
  nc <- sum(ncod); i3 <- 3L * seq_len(nc)
  b1 <- ref[i3 - 2L]; b2 <- ref[i3 - 1L]; b3 <- ref[i3]
  G <- list(name = nm, chr = chr, minus = minus, gpos = gpos, gposv = unlist(gpos, use.names = FALSE), ncod = ncod, coff = coff,
            ref = ref, upv = upv, dnv = dnv, ctx = 16L * (upv - 1L) + 4L * (ref - 1L) + (dnv - 1L) + 1L,
            b1 = b1, b2 = b2, b3 = b3, refaa = AA[cidx(b1, b2, b3)], cgene = rep.int(seq_along(gidx), ncod), cnum = sequence(ncod), nc = nc)
  L(sprintf("  structures built: %d codons in %.1f min", nc, mins_since(t0)))
  G
}
smat_from <- function(ms) {
  if (!all(c("name", "mle") %in% names(ms))) stop("mle_submodel has no name and mle columns [check]")
  sm <- setNames(ms$mle, ms$name); sm["TTT>TGT"] <- 1; sm <- sm * sm["t"]; sm <- sm[setdiff(names(sm), c("wmis", "wnon", "wspl", "t"))]
  S <- matrix(NA_real_, 64L, 4L)
  for (u in 1:4) for (r in 1:4) for (d in 1:4) for (a in 1:4) if (a != r)
    S[16L * (u - 1L) + 4L * (r - 1L) + (d - 1L) + 1L, a] <- sm[paste0(B[u], B[r], B[d], ">", B[u], B[a], B[d])]
  if (sum(is.na(S)) != 64L) stop("trinucleotide rate names missing in mle_submodel: ", sum(is.na(S)) - 64L)
  S
}
codon_rates <- function(G, S) {
  nc <- G$nc; Esyn <- EmisAll <- EmisST <- numeric(nc); bb <- list(G$b1, G$b2, G$b3); stopc <- G$refaa == "*"
  for (p in 1:3) { bi <- 3L * seq_len(nc) - 3L + p; rb <- bb[[p]]; cx <- G$ctx[bi]
    for (a in 1:4) { nb <- bb; nb[[p]] <- rep.int(a, nc); maa <- AA[cidx(nb[[1]], nb[[2]], nb[[3]])]
      r <- S[cbind(cx, a)]; r[rb == a] <- 0
      syn <- maa == G$refaa & !stopc; mis <- maa != G$refaa & maa != "*" & !stopc
      Esyn <- Esyn + r * syn; EmisAll <- EmisAll + r * mis; EmisST <- EmisST + r * (mis & !(maa %in% c("S", "T"))) } }
  list(Esyn = Esyn, EmisAll = EmisAll, EmisST = EmisST)
}
obs_codons <- function(G, am) {
  need <- c("gene", "pos", "aachange", "ref_cod", "ref3_cod", "impact"); miss <- setdiff(need, names(am))
  if (length(miss)) stop("annotmuts columns missing: ", paste(miss, collapse = ", "))
  am <- am[am$impact %in% c("Missense", "Synonymous") & am$gene %in% G$name, ]
  gi <- match(am$gene, G$name); bidx <- rep(NA_integer_, nrow(am)); posi <- as.integer(am$pos)
  sp <- split(seq_len(nrow(am)), gi)
  for (k in names(sp)) { i <- as.integer(k); rr <- sp[[k]]; bidx[rr] <- match(posi[rr], G$gpos[[i]]) }
  cod <- (bidx - 1L) %/% 3L + 1L; cg <- G$coff[gi] + cod; gb <- 3L * G$coff[gi] + bidx
  our_ref <- B[G$ref[gb]]; our_ref3 <- paste0(B[G$upv[gb]], our_ref, B[G$dnv[gb]])
  nch <- nchar(am$aachange); aa_ref <- substr(am$aachange, 1, 1); aa_mut <- substr(am$aachange, nch, nch)
  aa_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", am$aachange)))
  agree <- !is.na(bidx) & our_ref == toupper(am$ref_cod) & our_ref3 == toupper(am$ref3_cod) & aa_num == cod & aa_ref == G$refaa[cg]
  agree[is.na(agree)] <- FALSE
  mis <- am$impact == "Missense"; okc <- agree & G$refaa[cg] != "*"; okc[is.na(okc)] <- FALSE
  ex <- data.frame(gene = am$gene, pos = am$pos, aachange = am$aachange, ref_cod = am$ref_cod, ref3_cod = am$ref3_cod, impact = am$impact,
                   script_ref = our_ref, script_ref3 = our_ref3, script_codon = cod, script_aa = G$refaa[cg])[!agree, ][1:min(5, sum(!agree)), ]
  list(Omis = tabulate(cg[okc & mis], G$nc), OmisST = tabulate(cg[okc & mis & !(aa_mut %in% c("S", "T"))], G$nc),
       Osyn = tabulate(cg[okc & !mis], G$nc), n = nrow(am), n_agree = sum(agree), frac = if (nrow(am)) mean(agree) else 0, examples = ex)
}
run_data <- function(run, G, label) {
  t0 <- Sys.time(); S <- smat_from(run$mle_submodel); CR <- codon_rates(G, S); OB <- obs_codons(G, run$annotmuts)
  gm <- run$genemuts; relmr <- setNames(gm$exp_syn_cv / gm$exp_syn, gm$gene_name)
  L(sprintf("  run %s: code check C-A (a): %d of %d Missense and Synonymous rows of the structure genes agree (%.5f); %.1f min",
            label, OB$n_agree, OB$n, OB$frac, mins_since(t0)))
  if (OB$frac < 0.999) { L("  first rows that disagree:"); LP(OB$examples); stop("code check C-A (a) FAILED (below 0.999)") }
  list(CR = CR, OB = OB, relmr = relmr)
}
build_pos <- function(G, MAP, mapping, RD, drop_genes = character(0), quiet = FALSE) {
  cod <- if (mapping == "strict") MAP$codon_strict else MAP$codon_window
  gi <- match(MAP$gene, G$name); ok <- !is.na(cod) & !is.na(gi); ok[ok] <- cod[ok] <= G$ncod[gi[ok]]
  X <- MAP[ok, ]; gi <- gi[ok]; cg <- G$coff[gi] + cod[ok]
  bad <- G$refaa[cg] != X$residue; if (any(bad)) stop("the residue at the mapped codon differs from UniProt in ", sum(bad), " positions")
  rr <- unname(RD$relmr[G$name[gi]]); okr <- is.finite(rr) & rr > 0; drp <- G$name[gi] %in% drop_genes
  dup <- duplicated(cg) | duplicated(cg, fromLast = TRUE); isN <- X$residue == "N"
  P <- data.frame(gi = gi, res = as.integer(!isN), role = X$role, cg = cg,
                  E = rr * ifelse(isN, RD$CR$EmisAll[cg], RD$CR$EmisST[cg]), O = ifelse(isN, RD$OB$Omis[cg], RD$OB$OmisST[cg]))
  keep <- okr & !drp & !dup
  if (!quiet) L(sprintf("  positions (%s mapping): %d mapped; left out: no gene rate %d, codon reached twice %d, driver genes %d; kept %d (role 1 %d, role 2 %d, role 3 %d, role 0 %d)",
                        mapping, nrow(P), sum(!okr), sum(dup), sum(drp), sum(keep), sum(keep & P$role == 1L), sum(keep & P$role == 2L), sum(keep & P$role == 3L), sum(keep & P$role == 0L)))
  P[keep, ]
}

# ---------------- statistics ----------------
agg <- function(P, sel, g) {
  k <- P$gi[sel] * 2L + P$res[sel]; gg <- g[sel]; O <- P$O[sel]; E <- P$E[sel]
  m <- rowsum(cbind(O * gg, E * gg, O * (!gg), E * (!gg)), k); kk <- as.integer(rownames(m))
  list(Og = m[, 1], Eg = m[, 2], Ob = m[, 3], Eb = m[, 4], sg = kk %/% 2L, sres = kk %% 2L)
}
subA <- function(A, i) lapply(A, function(v) v[i])
core_stats <- function(A, seed, nboot = NBOOT) {
  inf <- A$Eg > 0 & A$Eb > 0 & (A$Og + A$Ob) > 0
  out <- list(tested = FALSE, log2rr = NA_real_, lo95 = NA_real_, hi95 = NA_real_, lo90 = NA_real_, hi90 = NA_real_, p = NA_real_,
              ev_g = sum(A$Og[inf]), ev_b = sum(A$Ob[inf]), exp_g_h0 = 0, n_strata = sum(inf), n_genes = length(unique(A$sg[inf])), boot_dropped = NA_integer_)
  if (!any(inf)) return(out)
  Og <- A$Og[inf]; Ob <- A$Ob[inf]; Eg <- A$Eg[inf]; Eb <- A$Eb[inf]; sg <- A$sg[inf]; Tt <- Eg + Eb
  num <- Og * Eb / Tt; den <- Ob * Eg / Tt; n <- Og + Ob; pe <- Eg / Tt; e0 <- sum(n * pe); v0 <- sum(n * pe * (1 - pe))
  out$exp_g_h0 <- e0; out$log2rr <- log2(sum(num) / sum(den))
  out$p <- if (v0 > 0) 2 * pnorm(-abs((sum(Og) - e0) / sqrt(v0))) else NA_real_; out$tested <- (e0 >= 5) && (v0 > 0)
  ug <- sort(unique(sg)); gi <- match(sg, ug); gn <- as.numeric(rowsum(num, gi)); gd <- as.numeric(rowsum(den, gi)); Gn <- length(ug)
  set.seed(seed); bl <- numeric(0); per <- 100L
  for (ch in seq_len(ceiling(nboot / per))) { k <- min(per, nboot - (ch - 1L) * per)
    idx <- sample.int(Gn, Gn * k, replace = TRUE)
    bl <- c(bl, log2(colSums(matrix(gn[idx], nrow = Gn)) / colSums(matrix(gd[idx], nrow = Gn)))) }
  fin <- is.finite(bl); out$boot_dropped <- sum(!fin)
  if (sum(fin) >= 0.9 * nboot) { q <- quantile(bl[fin], c(0.025, 0.975, 0.05, 0.95), names = FALSE, type = 7)
    out$lo95 <- q[1]; out$hi95 <- q[2]; out$lo90 <- q[3]; out$hi90 <- q[4] }
  out
}
lab <- function(d) ifelse(!d$tested, NT,
  ifelse(!is.na(d$fdr) & d$fdr < FDRCUT & is.finite(d$log2rr) & abs(d$log2rr) >= CH, ifelse(d$log2rr > 0, "Enriched", "Depleted"),
  ifelse(is.finite(d$lo90) & is.finite(d$hi90) & d$lo90 >= -EQ & d$hi90 <= EQ, "Equivalent", "Inconclusive")))
SCOLS <- c("tested", "log2rr", "lo95", "hi95", "lo90", "hi90", "p", "ev_g", "ev_b", "exp_g_h0", "n_strata", "n_genes", "boot_dropped")
run_units <- function(P, GL, BG, seed, tiers = 1:3) {
  rows <- list()
  for (t in tiers) for (grp in c("N", "ST", "both")) {
    rs <- switch(grp, N = P$res == 0L, ST = P$res == 1L, both = rep(TRUE, nrow(P)))
    sel <- rs & (GL[[t]] | BG[[t]]); s <- core_stats(agg(P, sel, GL[[t]]), seed)
    rows[[length(rows) + 1L]] <- data.frame(tier = paste0("T", t), group = grp, s[SCOLS], n_pos_glyco = sum(sel & GL[[t]]), n_pos_bg = sum(sel & BG[[t]])) }
  d <- do.call(rbind, rows); d$fdr <- NA_real_; i <- which(d$tested & !is.na(d$p)); if (length(i)) d$fdr[i] <- p.adjust(d$p[i], "BH")
  d$raw_label <- lab(d); d
}
real_units <- function(P, seed, tiers = 1:3) run_units(P, lapply(1:3, function(t) P$role == t), lapply(1:3, function(t) P$role == 0L), seed, tiers)
fmt_row <- function(d) sprintf("%s %-4s log2RR %7.3f (95%% CI %7.3f to %7.3f; 90%% CI %7.3f to %7.3f) p %.3g FDR %s events %d vs %d, genes %d -> %s",
  d$tier, d$group, d$log2rr, d$lo95, d$hi95, d$lo90, d$hi90, d$p, ifelse(is.na(d$fdr), "NA", sprintf("%.3g", d$fdr)), as.integer(d$ev_g), as.integer(d$ev_b), as.integer(d$n_genes), d$raw_label)

# ---------------- dNdScv runs ----------------
reduce_run <- function(x, mins, nin) {
  am <- x$annotmuts; keepc <- intersect(c("sampleID", "chr", "pos", "ref", "mut", "gene", "aachange", "ref_cod", "mut_cod", "ref3_cod", "mut3_cod", "impact"), names(am))
  if (!all(c("gene_name", "exp_syn", "exp_syn_cv") %in% names(x$genemuts))) stop("genemuts has no exp_syn and exp_syn_cv columns [check]")
  list(annotmuts = am[am$impact %in% c("Missense", "Synonymous"), keepc], genemuts = x$genemuts[, c("gene_name", "exp_syn", "exp_syn_cv")],
       mle_submodel = x$mle_submodel, globaldnds = x$globaldnds, sel_cv = x$sel_cv, nexcl = length(x$exclsamples),
       nwrong = if (is.null(x$wrongmuts)) 0L else NROW(x$wrongmuts), n_input = nin, minutes = mins, impact_counts = table(x$annotmuts$impact))
}
get_run <- function(name, M, maxg, maxc, keep_full = FALSE) {
  fr <- file.path(RUNS, paste0(name, "_reduced.rds")); ff <- file.path(RUNS, paste0(name, "_full.rds"))
  if (file.exists(fr) && (!keep_full || file.exists(ff))) { L("dNdScv run ", name, ": saved run found, not run again"); return(readRDS(fr)) }
  L(sprintf("dNdScv run %s: %d SNV rows, %d samples, max_muts_per_gene_per_sample %g, max_coding_muts_per_sample %g; started %s",
            name, nrow(M), length(unique(M$sampleID)), maxg, maxc, format(Sys.time(), "%H:%M:%S")))
  t0 <- Sys.time()
  x <- dndscv::dndscv(M, refdb = "hg38", outmats = TRUE, max_muts_per_gene_per_sample = maxg, max_coding_muts_per_sample = maxc)
  mins <- mins_since(t0); L(sprintf("dNdScv run %s finished in %.1f minutes (measured)", name, mins))
  if (keep_full) { saveRDS(x, ff); L("  full output saved (for sitednds): ", ff) }
  red <- reduce_run(x, mins, nrow(M)); saveRDS(red, fr); rm(x); invisible(gc()); red
}
run_summary <- function(name, r) {
  L(sprintf("run %s: input SNVs %d; samples excluded by dNdScv %d; rows with a wrong reference base %d (%.4f); run minutes %.1f",
            name, r$n_input, r$nexcl, r$nwrong, r$nwrong / r$n_input, r$minutes))
  L("  annotated impacts: ", paste(names(r$impact_counts), as.integer(r$impact_counts), collapse = "; "))
  L("  global dN/dS (globaldnds):"); LP(r$globaldnds)
}
read_mut <- function(g) {
  M <- read.delim(gzfile(MUT), colClasses = c("character", "character", "integer", "character", "character"), quote = "")
  if (g[["chr_rule"]] == "strip") M$chr <- sub("^chr", "", M$chr)
  M
}
load_common <- function() {
  if (!file.exists(GATE2)) stop("section 2 has not finished (its gate file is missing)")
  g <- read_kv(GATE2); if (!identical(md5(MAPF), g[["map_md5"]])) stop("the mapping table changed after section 2")
  e <- load_refcds(); R <- get("RefCDS", envir = e); refnames <- vapply(R, function(x) as.character(x$gene_name), "")
  MAP <- read.delim(gzfile(MAPF), colClasses = c("character", "character", "integer", "character", "integer", "integer", "integer"), quote = "")
  G <- prep_genes(R, match(unique(c(MAP$gene, HOT$gene)), refnames)); rm(R, e); invisible(gc())
  list(g = g, MAP = MAP, G = G, primary = g[["primary_mapping"]], other = if (g[["primary_mapping"]] == "strict") "window" else "strict", bail = readLines(BAILF))
}

# ---------------- keys for sitednds ----------------
enum_keys <- function(G, cg, rule) {
  cg <- unique(cg); if (!length(cg)) return(character(0))
  out <- list(); gi <- G$cgene[cg]; cn <- G$cnum[cg]; raa <- G$refaa[cg]; bb <- list(G$b1[cg], G$b2[cg], G$b3[cg])
  for (p in 1:3) { gb <- 3L * (cg - 1L) + p; rb <- bb[[p]]
    for (a in 1:4) { nb <- bb; nb[[p]] <- rep.int(a, length(cg)); maa <- AA[cidx(nb[[1]], nb[[2]], nb[[3]])]
      syn <- maa == raa & raa != "*"; mis <- maa != raa & maa != "*" & raa != "*"
      sel <- switch(rule, misAll = mis, misST = mis & !(maa %in% c("S", "T")), syn = syn) & rb != a
      k <- which(sel); if (!length(k)) next
      g <- gi[k]; rr <- B[rb[k]]; mm <- rep(B[a], length(k)); mn <- G$minus[g]
      gref <- ifelse(mn, unname(COMP[rr]), rr); gmut <- ifelse(mn, unname(COMP[mm]), mm)
      up <- B[G$upv[gb[k]]]; dn <- B[G$dnv[gb[k]]]
      out[[length(out) + 1L]] <- paste(G$chr[g], G$gposv[gb[k]], gref, gmut, G$name[g], paste0(raa[k], cn[k], maa[k]), paste0(up, rr, dn), paste0(up, mm, dn), sep = ":") } }
  unique(unlist(out))
}
kn_extract <- function(x) {
  r <- list(obs = NA_real_, exp = NA_real_, est = NA_real_)
  if (inherits(x, "htest")) { r$obs <- as.numeric(x$statistic[1]); r$est <- as.numeric(x$estimate[1]); if (is.finite(r$est) && r$est > 0) r$exp <- r$obs / r$est; return(r) }
  v <- unlist(x); if (is.null(names(v))) return(r); nm <- tolower(names(v)); num <- suppressWarnings(as.numeric(v))
  pick <- function(pat) { i <- grep(pat, nm); i <- i[is.finite(num[i])]; if (length(i)) num[i[1]] else NA_real_ }
  r$obs <- pick("obs"); r$exp <- pick("exp"); r$est <- pick("dnds|estimate|mle|ratio")
  if (!is.finite(r$est) && is.finite(r$obs) && is.finite(r$exp) && r$exp > 0) r$est <- r$obs / r$exp
  if (!is.finite(r$exp) && is.finite(r$obs) && is.finite(r$est) && r$est > 0) r$exp <- r$obs / r$est
  r
}
anchor_one <- function(full, keys, oo, oe, label) {
  t0 <- Sys.time(); kn <- dndscv::sitednds(full, site_list = keys)$globaldnds_knownsites
  L(sprintf("  C-A (b) %s: %d listed changes; sitednds took %.1f min", label, length(keys), mins_since(t0)))
  L("    globaldnds_knownsites: ", paste(trimws(capture.output(str(kn))), collapse = " | "))
  x <- kn_extract(kn)
  ok_e <- is.finite(x$exp) && abs(oe / x$exp - 1) < 1e-6
  ok_o <- if (is.finite(x$obs)) abs(oo - x$obs) <= 0.001 * x$obs else NA
  pass <- isTRUE(ok_e) && !isFALSE(ok_o)
  L(sprintf("    this script: observed %d, expected %.6f; sitednds: observed %s, expected %s -> %s", as.integer(oo), oe,
            format(x$obs), format(x$exp, digits = 12), if (pass) "PASSED" else "FAILED"))
  data.frame(list = label, n_keys = length(keys), script_obs = oo, script_exp = oe, sitednds_obs = x$obs, sitednds_exp = x$exp, sitednds_est = x$est, passed = pass)
}
hot_codons <- function(G) {
  gi <- match(HOT$gene, G$name); cg <- ifelse(is.na(gi), NA_integer_, ifelse(HOT$codon > G$ncod[gi], NA_integer_, G$coff[gi] + HOT$codon))
  ok <- !is.na(cg) & G$refaa[cg] == HOT$aa; ok[is.na(ok)] <- FALSE
  for (i in which(!ok)) L("  hotspot left out: ", HOT$gene[i], " ", HOT$aa[i], HOT$codon[i], " (", if (is.na(gi[i])) "gene not in the structures" else paste0("RefCDS residue ", G$refaa[cg[i]]), ")")
  data.frame(gene = HOT$gene[ok], gi = gi[ok], aa = HOT$aa[ok], cg = cg[ok])
}

# ===================== section 2: structure check =====================
sec2 <- function() {
  L("R ", R.version.string, "; root ", ROOT)
  for (p in c("dndscv", "Biostrings", "GenomicRanges", "IRanges", "Rsamtools", "rtracklayer", "seqinr", "MASS", "poilog", "plyr"))
    L(sprintf("  package %-14s %s", p, tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")))
  d <- packageDescription("dndscv")
  L("  dndscv source: ", if (is.null(d$RemoteUsername)) "?" else paste0(d$RemoteUsername, "/", d$RemoteRepo), "; RemoteSha ", if (is.null(d$RemoteSha)) "(none)" else d$RemoteSha,
    "; Packaged ", if (is.null(d$Packaged)) "?" else d$Packaged)
  suppressPackageStartupMessages(library(dndscv))
  for (fn in c("dndscv", "sitednds")) { fa <- formals(get(fn, envir = asNamespace("dndscv"))); L("  arguments of ", fn, ":")
    for (n in names(fa)) L("    ", n, " = ", paste(deparse(fa[[n]]), collapse = " ")) }
  fcode <- file.path(LOGS, "P2_23_section2_function_code.txt")
  writeLines(c("dndscv::dndscv", capture.output(print(get("dndscv", envir = asNamespace("dndscv")))), "", "dndscv::sitednds",
               capture.output(print(get("sitednds", envir = asNamespace("dndscv"))))), fcode)
  L("  code of dndscv and sitednds written to ", fcode, " (md5 ", md5(fcode), ")")
  L("inputs from section 1:"); check_gate1()
  e <- load_refcds(); L("objects in refcds_GRCh38_hg38: ", paste(ls(e), collapse = ", "))
  if (!exists("RefCDS", envir = e, inherits = FALSE)) stop("RefCDS not found in refcds_GRCh38_hg38 [check]")
  R <- get("RefCDS", envir = e); L("RefCDS entries (genes): ", length(R)); L("fields of entry 1: ", paste(names(R[[1]]), collapse = ", "))
  miss <- setdiff(c("gene_name", "CDS_length", "chr", "strand", "intervals_cds", "seq_cds", "seq_cds1up", "seq_cds1down"), names(R[[1]]))
  if (length(miss)) stop("RefCDS fields missing: ", paste(miss, collapse = ", "))
  L("class of seq_cds ", class(R[[1]]$seq_cds)[1], "; intervals_cds ", class(R[[1]]$intervals_cds)[1], " ", paste(dim(as.matrix(R[[1]]$intervals_cds)), collapse = " x "),
    "; strand values ", paste(sort(unique(vapply(R, function(x) as.character(x$strand), ""))), collapse = " "))
  urchr <- sort(unique(vapply(R, function(x) as.character(x$chr), ""))); L("RefCDS chromosome names (", length(urchr), "): ", paste(urchr, collapse = " "))
  refnames <- vapply(R, function(x) as.character(x$gene_name), ""); dupn <- unique(refnames[duplicated(refnames)])
  L("gene names used by more than one RefCDS entry: ", length(dupn), " (left out of the mapping)")
  M <- read.delim(gzfile(MUT), colClasses = c("character", "character", "integer", "character", "character"), quote = "")
  L("mutation table: ", nrow(M), " SNV rows, ", length(unique(M$sampleID)), " samples; columns ", paste(names(M), collapse = " "))
  mchr <- sort(unique(M$chr)); L("MAF chromosome names: ", paste(mchr, collapse = " "))
  rule <- if (!any(grepl("^chr", urchr)) && all(grepl("^chr", mchr))) "strip" else "keep"
  mc <- if (rule == "strip") sub("^chr", "", M$chr) else M$chr
  L("chromosome rule: ", rule, "; SNV rows on chromosomes without a RefCDS gene: ", sum(!(mc %in% urchr))); rm(M, mc)
  CS <- read.delim(CASES, colClasses = c("character", "character", "integer", "integer"))
  bp <- sort(tapply(CS$n_snv, CS$project, sum), decreasing = TRUE); loo <- names(bp)[1:5]
  L("cases ", nrow(CS), "; hypermutated ", sum(CS$hypermutated == 1L), "; projects ", length(bp))
  L("D3 B: the 5 projects with the most SNVs in the dNdScv input: ", paste(sprintf("%s (%d)", loo, as.integer(bp[1:5])), collapse = ", "))
  U <- read.delim(UNIF, colClasses = "character", quote = "", na.strings = character(0))
  U$gi <- match(U$gene, refnames); U$gi[U$gene == "" | U$gene %in% dupn] <- NA_integer_
  L("UniProt accessions in the position list: ", nrow(U), "; with a RefCDS gene of the same name: ", sum(!is.na(U$gi)))
  tr <- rep(NA_character_, length(R)); t0 <- Sys.time(); ug <- unique(U$gi[!is.na(U$gi)])
  for (j in ug) tr[j] <- translate_gene(seqchar(R[[j]]$seq_cds))
  L(sprintf("  translated %d RefCDS genes (%d not translatable) in %.1f min", sum(!is.na(tr[ug])), sum(is.na(tr[ug])), mins_since(t0)))
  U$strict <- !is.na(U$gi) & !is.na(tr[U$gi]) & tr[U$gi] == U$seq; U$strict[is.na(U$strict)] <- FALSE
  dg <- U$gi[U$strict][duplicated(U$gi[U$strict])]
  if (length(dg)) { L("  RefCDS genes identical to more than one accession (left out): ", length(unique(dg))); U$strict[U$gi %in% dg] <- FALSE }
  L("strict rule: accessions mapped ", sum(U$strict))
  P <- read.delim(gzfile(POSF), colClasses = c("character", "character", "integer", "character", "integer"), quote = "")
  ui <- match(P$accession, U$accession); P$gene_ref <- refnames[U$gi[ui]]
  P$codon_strict <- ifelse(U$strict[ui], P$position, NA_integer_)
  sgenes <- unique(U$gi[U$strict]); cand <- which(!U$strict & !is.na(U$gi) & !is.na(tr[U$gi]) & !(U$gi %in% sgenes))
  dup2 <- U$gi[cand][duplicated(U$gi[cand])]; cand <- cand[!(U$gi[cand] %in% dup2)]
  L("window rule: candidate accessions ", length(cand), " (", length(unique(dup2)), " genes with more than one candidate accession left out)")
  P$codon_window <- P$codon_strict; t0 <- Sys.time(); pia <- split(seq_len(nrow(P)), P$accession)
  for (i in cand) { rr <- pia[[U$accession[i]]]; if (is.null(rr)) next; P$codon_window[rr] <- win_map(U$seq[i], tr[U$gi[i]], P$position[rr]) }
  L(sprintf("  window search took %.1f min", mins_since(t0)))
  for (ro in 0:3) for (rs in c("N", "S/T")) { s <- P$role == ro & (if (rs == "N") P$residue == "N" else P$residue != "N")
    L(sprintf("  positions role %d %-3s: %d; strict %d; strict + window %d", ro, rs, sum(s), sum(s & !is.na(P$codon_strict)), sum(s & !is.na(P$codon_window)))) }
  t1 <- P$role == 1L; cs <- mean(!is.na(P$codon_strict[t1])); cw <- mean(!is.na(P$codon_window[t1]))
  primary <- if (cs >= 0.60) "strict" else "window"
  L(sprintf("T1 sequon positions covered: strict %.3f, strict + window %.3f; fixed rule: strict if >= 0.60 -> PRIMARY MAPPING = %s", cs, cw, primary))
  MAPX <- P[!is.na(P$codon_window), c("accession", "gene_ref", "position", "residue", "role", "codon_strict", "codon_window")]; names(MAPX)[2] <- "gene"
  zz <- gzfile(MAPF, "w"); write.table(MAPX, zz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(zz)
  L("mapping table: ", nrow(MAPX), " rows (md5 ", md5(MAPF), ")")
  writeLines(c(paste0("primary_mapping\t", primary), paste0("chr_rule\t", rule), paste0("loo_projects\t", paste(loo, collapse = ",")),
               sprintf("cov_strict\t%.4f", cs), sprintf("cov_window\t%.4f", cw), paste0("map_md5\t", md5(MAPF))), GATE2)
  L("gate file written: ", GATE2)
}

# ===================== section 3: dNdScv runs and controls =====================
sec3 <- function() {
  check_seeds(); suppressPackageStartupMessages(library(dndscv))
  C0 <- load_common(); G <- C0$G; MAP <- C0$MAP; L("primary mapping ", C0$primary, "; chromosome rule ", C0$g[["chr_rule"]])
  M <- read_mut(C0$g); CS <- read.delim(CASES, colClasses = c("character", "character", "integer", "integer")); hyp <- CS$case_id[CS$hypermutated == 1L]
  RM <- get_run("main", M, 3, 3000, keep_full = TRUE)
  L(sprintf("measured time of the first run: %.1f minutes; the next two runs should take about as long", RM$minutes))
  RC <- get_run("nocaps", M, 1e9, 1e9); RH <- get_run("nohyper", M[!(M$sampleID %in% hyp), ], 3, 3000); rm(M); invisible(gc())
  run_summary("main", RM); run_summary("nocaps", RC); run_summary("nohyper", RH)
  # C-P1
  sc <- RM$sel_cv; L("C-P1 columns of sel_cv: ", paste(names(sc), collapse = ", "))
  qcol <- intersect(c("qglobal_cv", "qallsubs_cv"), names(sc))[1]
  if (is.na(qcol)) stop("C-P1: sel_cv has neither qglobal_cv nor qallsubs_cv")
  L("C-P1 uses column ", qcol, " (dated correction i: qglobal_cv exists only when indels are in the input; without indels the combined missense + truncating q is qallsubs_cv)")
  q <- sc[[qcol]][match(c("TP53", "KRAS", "PIK3CA"), sc$gene_name)]
  if (length(q) != 3L || anyNA(q)) stop("C-P1: a driver gene has no q value (guard against an empty pass): ", paste(format(q), collapse = " "))
  cp1 <- all(q < 0.1)
  L(sprintf("C-P1 drivers: %s TP53 %s, KRAS %s, PIK3CA %s (rule: each < 0.1) -> %s", qcol, format(q[1], digits = 3), format(q[2], digits = 3), format(q[3], digits = 3), if (cp1) "PASSED" else "FAILED"))
  det <- sc[!is.na(sc[[qcol]]) & sc[[qcol]] < 0.1, ]; det$bailey <- det$gene_name %in% C0$bail
  L("  genes with ", qcol, " < 0.1: ", nrow(det), "; of them Bailey driver genes ", sum(det$bailey), " (of ", length(C0$bail), "; descriptive)")
  write.table(det, file.path(TABLES, "P2_23_drivers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  # structures for the main run
  RD <- run_data(RM, G, "main"); P <- build_pos(G, MAP, C0$primary, RD)
  t1genes <- sort(unique(P$gi[P$role == 1L])); HC <- hot_codons(G)
  # C-A (b)
  full <- readRDS(file.path(RUNS, "main_full.rds")); set.seed(SEED_ANCHOR); ag <- sort(sample(t1genes, min(200L, length(t1genes))))
  bsel <- P$role == 0L & P$gi %in% ag
  k1 <- c(enum_keys(G, P$cg[bsel & P$res == 0L], "misAll"), enum_keys(G, P$cg[bsel & P$res == 1L], "misST"))
  cgA <- which(G$cgene %in% ag); rA <- RD$relmr[G$name[G$cgene[cgA]]]
  k3 <- enum_keys(G, HC$cg, "misAll"); rH <- RD$relmr[G$name[G$cgene[HC$cg]]]
  an <- rbind(anchor_one(full, k1, sum(P$O[bsel]), sum(P$E[bsel]), "background changes, 200 genes"),
              anchor_one(full, enum_keys(G, cgA, "syn"), sum(RD$OB$Osyn[cgA]), sum(rA * RD$CR$Esyn[cgA]), "synonymous changes, same genes"),
              anchor_one(full, k3, sum(RD$OB$Omis[HC$cg]), sum(rH * RD$CR$EmisAll[HC$cg]), "hotspot codons"))
  rm(full); invisible(gc()); write.table(an, file.path(TABLES, "P2_23_anchor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!all(an$passed)) stop("code check C-A (b) FAILED: no label; the script stops here")
  L("C-A code check PASSED")
  # C-Cal
  cgT <- which(G$cgene %in% t1genes); wsyn <- sum(RD$OB$Osyn[cgT]) / sum(RD$relmr[G$name[G$cgene[cgT]]] * RD$CR$Esyn[cgT]); ccal <- is.finite(wsyn) && wsyn >= 0.8 && wsyn <= 1.25
  L(sprintf("C-Cal synonymous dN/dS over the codons of the %d genes with T1 positions: %.4f (rule 0.8 to 1.25) -> %s", length(t1genes), wsyn, if (ccal) "PASSED" else "FAILED (absolute readout not calibrated)"))
  # C-P2
  cgH <- which(G$cgene %in% unique(HC$gi)); aaH <- G$refaa[cgH]; giH <- G$cgene[cgH]; kp <- paste(giH, aaH) %in% paste(HC$gi, HC$aa)
  cgH <- cgH[kp]; aaH <- aaH[kp]; giH <- giH[kp]; gH <- cgH %in% HC$cg; OH <- RD$OB$Omis[cgH]; EH <- RD$relmr[G$name[giH]] * RD$CR$EmisAll[cgH]
  kH <- giH * 32L + match(aaH, AA_LET); mH <- rowsum(cbind(OH * gH, EH * gH, OH * (!gH), EH * (!gH)), kH)
  sH <- core_stats(list(Og = mH[, 1], Eg = mH[, 2], Ob = mH[, 3], Eb = mH[, 4], sg = as.integer(rownames(mH)) %/% 32L), SEED_BOOT)
  cp2 <- isTRUE(is.finite(sH$p) && sH$p < 0.05 && is.finite(sH$log2rr) && sH$log2rr >= CH)
  L(sprintf("C-P2 hotspot codons (%d of %d usable): log2RR %.3f, p %.3g, events %d vs %d (rule Enriched: p < 0.05 and log2 >= 0.32) -> %s",
            nrow(HC), nrow(HOT), sH$log2rr, sH$p, as.integer(sH$ev_g), as.integer(sH$ev_b), if (cp2) "PASSED" else "FAILED"))
  # C-N null
  t0 <- Sys.time(); k <- P$gi * 2L + P$res; nb <- max(k) + 1L; bg <- which(P$role == 0L); nbg <- tabulate(k[bg] + 1L, nb); nl <- list()
  for (r in seq_along(SEED_NULL)) {
    set.seed(SEED_NULL[r]); u <- runif(length(bg)); o <- order(k[bg], u); bo <- bg[o]; ko <- k[bo]; rk <- sequence(rle(ko)$lengths); skipped <- 0L
    GL <- lapply(1:3, function(t) { need <- tabulate(k[P$role == t] + 1L, nb); want <- need[ko + 1L]; ok <- nbg[ko + 1L] > want
      skipped <<- skipped + sum(need > 0L & nbg <= need); v <- logical(nrow(P)); v[bo[rk <= want & ok]] <- TRUE; v })
    BG <- lapply(GL, function(v) P$role == 0L & !v); d <- run_units(P, GL, BG, SEED_BOOT); d$rep <- r; nl[[r]] <- d
    L(sprintf("  C-N repeat %d: Changed units %d of %d tested; strata skipped (too few background positions) %d", r, sum(d$raw_label %in% c("Enriched", "Depleted")), sum(d$tested), skipped)) }
  nl <- do.call(rbind, nl); write.table(nl, file.path(TABLES, "P2_23_null_runs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cap <- do.call(rbind, lapply(split(nl, paste(nl$tier, nl$group)), function(x) data.frame(tier = x$tier[1], group = x$group[1], n_tested = sum(x$tested),
          n_changed = sum(x$raw_label %in% c("Enriched", "Depleted")), n_equivalent = sum(x$raw_label == "Equivalent"))))
  cap$fp_passed <- cap$n_changed <= 1L; cap$equivalent_allowed <- cap$n_equivalent >= 16L
  cap <- cap[order(cap$tier, match(cap$group, c("N", "ST", "both"))), ]; write.table(cap, CAPF, sep = "\t", quote = FALSE, row.names = FALSE)
  L(sprintf("C-N null (%.1f min): per unit, Changed in at most 1 of 20 = passed; Equivalent in at least 16 of 20 = Equivalent allowed", mins_since(t0))); LP(cap)
  # C-Pl planted
  t0 <- Sys.time(); Pz <- P; Pz$O[P$role == 1L] <- 0L; A <- agg(Pz, P$role %in% c(0L, 1L), P$role == 1L); rm(Pz)
  lam <- ifelse(A$Eb > 0, A$Ob / A$Eb, 0); LEV <- c(1.5, 0.67, 1.25, 0.80, 1.0); pl <- list()
  for (r in seq_along(SEED_PLANT)) { set.seed(SEED_PLANT[r]); draws <- lapply(LEV, function(lv) rpois(length(A$Eg), A$Eg * lam * lv))
    for (j in seq_along(LEV)) { Ap <- A; Ap$Og <- draws[[j]]
      dd <- do.call(rbind, lapply(c("N", "ST", "both"), function(grp) { i <- switch(grp, N = Ap$sres == 0L, ST = Ap$sres == 1L, both = rep(TRUE, length(Ap$sg)))
              data.frame(group = grp, core_stats(subA(Ap, i), SEED_BOOT)[SCOLS]) }))
      dd$fdr <- NA_real_; i <- which(dd$tested & !is.na(dd$p)); if (length(i)) dd$fdr[i] <- p.adjust(dd$p[i], "BH")
      dd$raw_label <- lab(dd); dd$rep <- r; dd$ratio <- LEV[j]; pl[[length(pl) + 1L]] <- dd } }
  pl <- do.call(rbind, pl); write.table(pl, file.path(TABLES, "P2_23_planted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  L(sprintf("C-Pl planted (%.1f min), T1 units, labels per ratio over 20 repeats:", mins_since(t0)))
  for (grp in c("N", "ST", "both")) for (lv in LEV) { x <- pl[pl$group == grp & pl$ratio == lv, ]
    L(sprintf("  T1 %-4s ratio %.2f: Enriched %2d, Depleted %2d, Equivalent %2d, Inconclusive %2d, not tested %2d", grp, lv, sum(x$raw_label == "Enriched"),
              sum(x$raw_label == "Depleted"), sum(x$raw_label == "Equivalent"), sum(x$raw_label == "Inconclusive"), sum(x$raw_label == NT))) }
  pb <- pl[pl$group == "both", ]; cpl <- sum(pb$ratio == 1.5 & pb$raw_label == "Enriched") >= 16L && sum(pb$ratio == 0.67 & pb$raw_label == "Depleted") >= 16L
  L("C-Pl (primary unit: Enriched at 1.5 and Depleted at 0.67 in at least 16 of 20 each) -> ", if (cpl) "PASSED" else "FAILED")
  pv <- function(x) if (x) "PASSED" else "FAILED"
  writeLines(c("anchor\tPASSED", paste0("cp1\t", pv(cp1)), paste0("cp2\t", pv(cp2)), paste0("ccal\t", pv(ccal)), sprintf("wsyn\t%.6f", wsyn),
               paste0("cpl\t", pv(cpl)), paste0("capability_table\t", CAPF)), VERD)
  L("controls verdict file written: ", VERD)
}

# ===================== section 4: real result =====================
sec4 <- function() {
  if (!file.exists(VERD)) stop("section 3 has not finished (controls verdict missing)")
  V <- read_kv(VERD); if (V[["anchor"]] != "PASSED") stop("code check C-A did not pass")
  L("controls: C-P1 ", V[["cp1"]], "; C-P2 ", V[["cp2"]], "; C-Cal ", V[["ccal"]], " (", V[["wsyn"]], "); C-Pl ", V[["cpl"]])
  CAP <- read.delim(CAPF); C0 <- load_common(); G <- C0$G; MAP <- C0$MAP
  RM <- readRDS(file.path(RUNS, "main_reduced.rds")); RC <- readRDS(file.path(RUNS, "nocaps_reduced.rds")); RH <- readRDS(file.path(RUNS, "nohyper_reduced.rds"))
  DM <- run_data(RM, G, "main"); DC <- run_data(RC, G, "nocaps"); DH <- run_data(RH, G, "nohyper")
  P <- build_pos(G, MAP, C0$primary, DM)
  d1 <- real_units(P, SEED_BOOT); d2 <- real_units(P, SEED_BOOT); d3 <- real_units(P, SEED_BOOT3)
  same12 <- isTRUE(all.equal(d1, d2, tolerance = 0)); L("runs 1 and 2 (same seed) identical: ", same12)
  L(sprintf("run 3 (seed %d): labels that differ from run 1: %d; largest change of a 90%% CI bound %.4f", SEED_BOOT3, sum(d3$raw_label != d1$raw_label),
            max(abs(c(d3$lo90 - d1$lo90, d3$hi90 - d1$hi90)), na.rm = TRUE)))
  guard <- V[["cp1"]] == "PASSED" && V[["cp2"]] == "PASSED"
  ci <- match(paste(d1$tier, d1$group), paste(CAP$tier, CAP$group))
  d1$final_label <- ifelse(d1$raw_label == NT, NT, ifelse(!guard | !CAP$fp_passed[ci], INS, ifelse(d1$raw_label == "Equivalent" & !CAP$equivalent_allowed[ci], INS, d1$raw_label)))
  L("MAIN RUN, ", C0$primary, " mapping (labels from run 1; final label after the control rules):")
  for (i in seq_len(nrow(d1))) L("  ", fmt_row(d1[i, ]), " | FINAL ", d1$final_label[i])
  pr <- d1[d1$tier == "T1" & d1$group == "both", ]
  L(sprintf("PRIMARY UNIT (T1 N + S/T, main run, %s mapping): log2RR %.3f (95%% CI %.3f to %.3f; 90%% CI %.3f to %.3f), FDR %.3g, %d glycosite and %d background events in %d genes -> FINAL %s",
            C0$primary, pr$log2rr, pr$lo95, pr$hi95, pr$lo90, pr$hi90, pr$fdr, as.integer(pr$ev_g), as.integer(pr$ev_b), as.integer(pr$n_genes), pr$final_label))
  res <- rbind(cbind(setup = "main", mapping = C0$primary, run = 1L, d1), cbind(setup = "main", mapping = C0$primary, run = 2L, d2, final_label = NA),
               cbind(setup = "main", mapping = C0$primary, run = 3L, d3, final_label = NA))
  beside <- list(list("main", C0$other, DM, character(0)), list("nocaps", C0$primary, DC, character(0)), list("nohyper", C0$primary, DH, character(0)),
                 list("nodrivers", C0$primary, DM, C0$bail))
  for (b in beside) { L("SET-UP ", b[[1]], ", ", b[[2]], " mapping (beside; raw labels, BH within the set-up):")
    Pb <- build_pos(G, MAP, b[[2]], b[[3]], b[[4]]); db <- real_units(Pb, SEED_BOOT); for (i in seq_len(nrow(db))) L("  ", fmt_row(db[i, ]))
    res <- rbind(res, cbind(setup = b[[1]], mapping = b[[2]], run = 1L, db, final_label = NA)) }
  write.table(res, file.path(TABLES, "P2_23_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  # absolute readout
  wsyn <- as.numeric(V[["wsyn"]]); cal <- V[["ccal"]] == "PASSED"; tg <- unique(P$gi[P$role == 1L]); ab <- list()
  for (grp in c("N", "ST", "both")) { rs <- switch(grp, N = P$res == 0L, ST = P$res == 1L, both = rep(TRUE, nrow(P)))
    for (cl in c("T1 glycosite", "background, same genes")) { s <- rs & (if (cl == "T1 glycosite") P$role == 1L else P$role == 0L & P$gi %in% tg)
      o <- sum(P$O[s]); e <- sum(P$E[s]); ct <- poisson.test(as.integer(round(o)), e)$conf.int
      ab[[length(ab) + 1L]] <- data.frame(group = grp, class = cl, observed = o, expected = e, dnds_raw = o / e, lo95_raw = ct[1], hi95_raw = ct[2],
        dnds_cal = if (cal) o / e / wsyn else NA_real_, lo95_cal = if (cal) ct[1] / wsyn else NA_real_, hi95_cal = if (cal) ct[2] / wsyn else NA_real_) } }
  ab <- do.call(rbind, ab); write.table(ab, file.path(TABLES, "P2_23_absolute.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  L("ABSOLUTE dN/dS (beside, no label; divided by the synonymous calibration ", sprintf("%.4f", wsyn), if (cal) "" else ", NOT calibrated: raw values only", "):"); LP(ab)
  L("dNdScv global dN/dS of the main run (beside):"); LP(RM$globaldnds)
  tryCatch({ full <- readRDS(file.path(RUNS, "main_full.rds"))
    kt <- c(enum_keys(G, P$cg[P$role == 1L & P$res == 0L], "misAll"), enum_keys(G, P$cg[P$role == 1L & P$res == 1L], "misST"))
    t0 <- Sys.time(); kn <- dndscv::sitednds(full, site_list = kt)$globaldnds_knownsites
    L(sprintf("sitednds globaldnds_knownsites for the T1 list (%d changes; %.1f min; package output, no correction factor, beside):", length(kt), mins_since(t0))); LP(kn)
    x <- kn_extract(kn); s <- P$role == 1L
    L(sprintf("  this script for the same list: observed %d, expected %.4f, ratio %.4f; sitednds ratio %s", as.integer(sum(P$O[s])), sum(P$E[s]), sum(P$O[s]) / sum(P$E[s]), format(x$est, digits = 6)))
    rm(full); invisible(gc()) }, error = function(e) L("  sitednds for the T1 list failed: ", conditionMessage(e), " (beside only; result above unchanged)"))
}

# ===================== section 5: leave one project out (D3 B) =====================
sec5 <- function() {
  if (!file.exists(file.path(TABLES, "P2_23_results.tsv"))) stop("section 4 has not finished")
  suppressPackageStartupMessages(library(dndscv)); V <- read_kv(VERD); C0 <- load_common(); G <- C0$G; MAP <- C0$MAP
  loo <- strsplit(C0$g[["loo_projects"]], ",")[[1]]; L("projects left out one at a time: ", paste(loo, collapse = ", "))
  M <- read_mut(C0$g); CS <- read.delim(CASES, colClasses = c("character", "character", "integer", "integer")); out <- list()
  RM <- readRDS(file.path(RUNS, "main_reduced.rds")); d0 <- real_units(build_pos(G, MAP, C0$primary, run_data(RM, G, "main"), quiet = TRUE), SEED_BOOT, tiers = 1)
  out[[1]] <- cbind(left_out = "(none: main run)", d0)
  for (pj in loo) { ids <- CS$case_id[CS$project == pj]
    r <- get_run(paste0("loo_", sub("^TCGA-", "", pj)), M[!(M$sampleID %in% ids), ], 3, 3000)
    d <- real_units(build_pos(G, MAP, C0$primary, run_data(r, G, pj), quiet = TRUE), SEED_BOOT, tiers = 1)
    L("WITHOUT ", pj, " (", length(ids), " cases; run ", sprintf("%.1f", r$minutes), " min):"); for (i in seq_len(nrow(d))) L("  ", fmt_row(d[i, ]))
    out[[length(out) + 1L]] <- cbind(left_out = pj, d) }
  out <- do.call(rbind, out); write.table(out, file.path(TABLES, "P2_23_leave_one_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  L("leave-one-out table written (descriptive; BH over the 3 T1 units of each run)")
}

if (START <= 2 && !run_section(2, "structure", sec2)) quit(save = "no", status = 1)
if (START <= 3 && !run_section(3, "controls", sec3)) quit(save = "no", status = 1)
if (START <= 4 && !run_section(4, "result", sec4)) quit(save = "no", status = 1)
if (START <= 5 && !run_section(5, "leaveout", sec5)) quit(save = "no", status = 1)
cat("SqnPAN_23 finished: all sections ran. Last log line should start with ---- END OF SECTION 5 ----\n")
