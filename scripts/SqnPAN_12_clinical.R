# SqnPAN_12_clinical.R - clinical links (analysis plan section 9, control C5)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Functions are in SqnPAN_12_helpers.R (same folder).
# One command runs everything:  source("scripts/SqnPAN_12_clinical.R")
# Logs:   results\logs\P2_12_section1_checks_controls.txt, P2_12_section2_results.txt, P2_12_section3_drivers.txt
# Tables: results\tables\P2_12_*.tsv ; verdict: data\processed\P2_12_controls_verdict.txt
#
# PRESPEC-START
# PRE-SPECIFICATION for SqnPAN_12_clinical.R (2026-10-31). Decisions: D1 A (primary contrast = loss-only against control-only, "both" its own level, as in analysis plan
# section 9), D2 A (Equivalent label allowed, TOST margin HR 0.80 to 1.25, guarded by rule E1), D3 A (control = all background N and S/T changes in the loss proteins, S<->T excluded; the
# stricter residue-class-matched version is sensitivity S1). Written before any survival value was read: only file structures and exposure group sizes per cancer type (no outcome crossed with exposure).
# 0. Inputs (md5 checked, the script stops on a mismatch): P2_06_mutations.tsv.gz, P2_06_cases.tsv, P2_06_proteins.tsv, TCGA-CDR-SupplementalTableS1.xlsx (sheet TCGA-CDR; Liu et al. 2018, PMID 29625055),
#    Bailey2018_TableS1_driver_genes.xlsx (sheet Table S1; Bailey et al. 2018, PMID 29625053). Packages survival and readxl (versions and the arguments of coxph and survreg are printed).
# 1. Cohort (analysis plan Amendment 3): PFI for the 29 types ACC BLCA BRCA CESC CHOL COAD ESCA GBM HNSC KIRC KIRP LGG LIHC LUAD LUSC MESO OV PAAD PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM.
#    LAML and PCPG are left out and counted; DLBC and KICH enter only sensitivity S7. A case is used if it is in the CDR (joined by TCGA patient barcode = SqnPAN_06 "submitter"), has >= 1 missense row
#    in its chosen MAF file (rule 4A), PFI (or OS) and time > 0, age and sex. Cases dropped at each step are counted. Covariates: age (age_at_initial_pathologic_diagnosis), sex (gender),
#    log(missense rows in the chosen file). Every Cox model is stratified by cancer type, Efron ties.
# 2. Exposure per case, from mutation rows in the chosen file (in_chosen_file = 1). LOSS = a T1 loss mutation (class loss or loss_and_gain with T1 in lost_tiers; any of N loss, S/T loss, X to P loss).
#    Loss protein set = proteins with >= 1 loss mutation in the cohort types (fixed from the 29 types). CONTROL = a missense row with res_role nonseqN or nonseqST (a residue in no sequon of any tier,
#    the Q2 background) in a protein of the loss protein set; S<->T substitutions at nonseqST are not counted. Groups: loss_only, control_only, both, neither. Decision D3 A.
#    Guards (the script stops if the cases differ): 1,460,792 chosen rows; T1 loss rows in 1,468 cases (SqnPAN_09 run; a different row count is only printed).
# 3. Primary test (D1 A): Cox model on control_only, loss_only and both cases (reference control_only), terms x_lo (loss_only) and x_both. HR of loss_only against control_only = PRIMARY, one test.
#    Labels: Associated_worse / Associated_better if two-sided Wald p < 0.05; else Equivalent if the TOST p for HR inside 0.80 to 1.25 is < 0.05; else Not_clear. The final label takes the control
#    results into account (section 4): Not_clear reads Inconclusive_test_not_sensitive if power at planted HR 1.5 is < 0.80; Equivalent reads Equivalent_not_reported if the Equivalent share at planted
#    HR 1.5 is > 0.05; Associated reads ..._null_check_FAILED if the null check fails; any label reads ..._control_FAILED if the age pipeline control fails.
# 4. Secondary family (Holm over 3, same label rule, each with its own power control): (a) both against control_only (term x_both of the primary model), (b) any loss (loss_only + both) against
#    control_only (term x_any), (c) any loss adjusted for carrying a control mutation, with the neither cases as reference (terms x_L and x_C; x_L is judged).
# 5. Reported, not judged: S1 stricter controls (a control row counts only if its residue class, N or S/T, occurs among the loss mutations of the same protein); S2 without hypermutated cases
#    (> 1,000 missense rows in the chosen file, analysis plan section 5); S3 without rescued rows (only_rescued = 1) in the exposure; S4 N losses only (cases with S/T or X to P losses removed);
#    S5 S/T losses only (cases with N or X to P losses removed); S6 OS instead of PFI; S7 plus DLBC and KICH; S8 unadjusted (stratified by type only); leave-one-cancer-type-out (29 fits of
#    the primary model and of x_any); per cancer type: any loss against control_only, tested only if >= 15 any-loss cases and >= 10 events among them, BH within the tested set, descriptive, no label;
#    proportional-hazards check (cox.zph) of the main model.
# 6. Controls (section 1, judged before any real HR is read):
#    PC pipeline control: age per 10 years in the main model must give HR > 1 with one-sided p < 0.025 (checks the join, coding and strata; it is not a sensitivity control). If it fails, every label
#       carries _control_FAILED.
#    PW power (analysis plan C5): Weibull model of PFI on type + age + sex + log missense count + the exposure columns, fitted on the real analysis cohort; event times simulated with a planted HR 1.0, 1.25, 1.5
#       on the judged exposure, censoring from a Weibull fitted to the reversed status (type-specific; pooled if that fit fails); the real Cox model refitted. Repeats: 500 for the primary contrast,
#       300 for each of the three secondary contrasts. Set A decides, set B shows seed stability. Sensitive if power (p < 0.05 and HR > 1) at planted HR 1.5 is >= 0.80; the weak level 1.25 is reported,
#       not judged. E1: share labelled Equivalent (TOST p < 0.05) at planted HR 1.5 must be <= 0.05. E2 (Equivalent at HR 1.0) and the false-positive share at HR 1.0 are reported, not judged.
#    NL null (primary contrast): the group labels of the control_only, loss_only and both cases are permuted among those cases inside cancer type x hypermutated status, real outcome and covariates
#       kept, primary model refitted, 500 permutations per set. The share with p < 0.05 for x_lo must be <= 0.07 (set A decides). The same null result is used for the secondary contrasts.
#    Determinism: 10 simulated repeats are recomputed and must be identical; section 2 computes all tables twice and requires identical content md5.
# 7. Driver genes (section 3, descriptive only): the genes of all proteins with a chosen-file T1 loss mutation (all 33 types) against the 299 driver genes of Bailey et al. Table S1. Alias step for
#    old gene symbols: FAM46D=TENT5D (Q8NEK8), H3F3A=H3-3A (P84243), H3F3C=H3-5 (Q6NXT2), HIST1H1C=H1-2 (P16403), HIST1H1E=H1-4 (P10412), RQCD1=CNOT9 (Q92600), TCEB1=ELOC (Q15369), WHSC1=NSD2 (O96028)
#    [from memory; each is guarded by the UniProt accession in P2_06_proteins.tsv, a failed guard leaves the symbol unmatched and prints it]; CBWD3 has no mapping. COADREAD = COAD and READ; PANCAN = any type.
#    Output: one row per overlapping gene (Bailey cancers, role, loss rows and cases, cancer types, whether a loss occurs in a type where Bailey lists the gene, protein changes) and the share of drivers
#    among the loss genes, among all reviewed genes and among genes with a T1 sequon (context only, no test).
# 8. Seeds: power A 32300000, null A 32400000, power B 32500000, null B 32600000 (seed = base + 20000 x spec + 2000 x HR index + repeat; null: base + permutation). Overlap check in the code.
# 9. Placement (analysis plan rule 8): a test that answered goes into Methods and Results; a test that could not answer gets one line in the main text, methods in the Supplement.
# PRESPEC-END

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(survival))
ROOT2 <- Sys.getenv("SQNPAN_ROOT", unset = local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f) == 1) normalizePath(file.path(dirname(normalizePath(f)), ".."), winslash = "/") else normalizePath(getwd(), winslash = "/") }))
RAW2 <- file.path(ROOT2, "data/raw"); PROC2 <- file.path(ROOT2, "data/processed"); TBL2 <- file.path(ROOT2, "results/tables"); LOG2 <- file.path(ROOT2, "results/logs")
F2 <- c(mut = file.path(PROC2, "P2_06_mutations.tsv.gz"), cases = file.path(PROC2, "P2_06_cases.tsv"), prot = file.path(PROC2, "P2_06_proteins.tsv"),
        cdr = file.path(RAW2, "TCGA-CDR-SupplementalTableS1.xlsx"), bailey = file.path(RAW2, "Bailey2018_TableS1_driver_genes.xlsx"))
MD5_2 <- c(mut = "3be28f9265312957c13ab2f0faa2b2b2", cases = "1a131193724b29368975e7cc8d02c396", prot = "f77f5ffabf9bdce910e4aa71f8099eb4",
           cdr = "a4591b2dcee39591f59e5e25a6ce75fa", bailey = "144a9ccbbf038236decc1226df24620b")
MD5_HELPERS2 <- "2514a4ea15145093214649c871daf46a"
FVER2 <- file.path(PROC2, "P2_12_controls_verdict.txt")

TYPES29 <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PRAD", "READ",
             "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
COV_MAIN2 <- c("age", "sex", "logmis"); EQ_HR2 <- 1.25
HRS2 <- c(1, 1.25, 1.5); NNULL2 <- 500L; POWER_MIN2 <- 0.80; E1_MAX2 <- 0.05; FWER_MAX2 <- 0.07
SEED2 <- c(powA = 32300000, nullA = 32400000, powB = 32500000, nullB = 32600000)
SPECS <- list(
  list(name = "P_loss_only_vs_control_only", data = "cc3", xs = c("x_lo", "x_both"), plant = "x_lo", target = "x_lo", nrep = 500L),
  list(name = "R_both_vs_control_only", data = "cc3", xs = c("x_lo", "x_both"), plant = "x_both", target = "x_both", nrep = 300L),
  list(name = "Q_any_loss_vs_control_only", data = "cc3", xs = "x_any", plant = "x_any", target = "x_any", nrep = 300L),
  list(name = "M_any_loss_adj_control_carrier", data = "cc4", xs = c("x_L", "x_C"), plant = "x_L", target = "x_L", nrep = 300L))
SPECNAMES <- vapply(SPECS, function(sp) sp$name, "")

HP2 <- file.path(ROOT2, "scripts/SqnPAN_12_helpers.R")
if (!identical(unname(tools::md5sum(HP2)), MD5_HELPERS2)) stop("SqnPAN_12_helpers.R md5 differs from the one written in this script. Stopping.")
source(HP2)
E2 <- new.env()

# ------------------------------------------------------------------ cohort and exposure builders ------------------------------------------------------------------
build_cohort <- function(types, ep) {
  d <- E2$cases; st <- list()
  add <- function(s, dd) st[[length(st) + 1L]] <<- data.frame(step = s, n = nrow(dd), stringsAsFactors = FALSE)
  d <- d[d$type %in% types, , drop = FALSE]; add("cases in the selected types", d)
  d <- d[d$in_cdr, , drop = FALSE]; add("... found in TCGA-CDR", d)
  d <- d[d$missense >= 1L, , drop = FALSE]; add("... with >= 1 missense row in the chosen file", d)
  d$event <- d[[ep]]; d$time <- d[[paste0(ep, ".time")]]
  d <- d[!is.na(d$event) & d$event %in% c(0L, 1L) & !is.na(d$time) & d$time > 0, , drop = FALSE]; add(paste0("... with ", ep, " and time > 0"), d)
  d <- d[!is.na(d$age) & !is.na(d$sex), , drop = FALSE]; add("... with age and sex (analysis cohort)", d)
  d$logmis <- log(d$missense); rownames(d) <- NULL
  list(d = d, steps = do.call(rbind, st))
}
assign_exp <- function(d, strict = FALSE, norescue = FALSE) {
  lo <- E2$loss; co <- E2$ctrl[E2$ctrl$accession %in% E2$protset, , drop = FALSE]
  if (norescue) { lo <- lo[lo$only_rescued == 0L, , drop = FALSE]; co <- co[co$only_rescued == 0L, , drop = FALSE] }
  if (strict) co <- co[(co$res_role == "nonseqN" & co$accession %in% E2$protN) | (co$res_role == "nonseqST" & co$accession %in% E2$protST), , drop = FALSE]
  cnt <- function(ids) { m <- match(ids, d$case_id); m <- m[!is.na(m)]; tabulate(m, nbins = nrow(d)) }
  d$nL <- cnt(lo$case_id); d$nC <- cnt(co$case_id)
  d$L <- d$nL > 0L; d$C <- d$nC > 0L
  d$LN <- d$case_id %in% lo$case_id[lo$ltype == "N"]; d$LST <- d$case_id %in% lo$case_id[lo$ltype == "ST"]; d$LO <- d$case_id %in% lo$case_id[lo$ltype == "odd"]
  d$grp <- ifelse(d$L & d$C, "both", ifelse(d$L, "loss_only", ifelse(d$C, "control_only", "neither")))
  d
}
add_dummies <- function(d) {
  d$type <- factor(as.character(d$type)); d$sex <- factor(as.character(d$sex), levels = c("FEMALE", "MALE"))
  d$x_lo <- as.numeric(d$grp == "loss_only"); d$x_both <- as.numeric(d$grp == "both"); d$x_any <- as.numeric(d$grp %in% c("loss_only", "both"))
  d$x_L <- as.numeric(d$L); d$x_C <- as.numeric(d$C); d$age10 <- d$age / 10
  d
}

run_section <- function(k, name, fun) {
  f <- file.path(LOG2, sprintf("P2_12_section%d_%s.txt", k, name)); dir.create(LOG2, showWarnings = FALSE, recursive = TRUE)
  sink(f, split = TRUE)
  ok <- tryCatch({ fun(); TRUE }, error = function(e) { cat("\nERROR in section", k, ":", conditionMessage(e), "\n"); FALSE },
                 finally = { cat(sprintf("# ---- END OF SECTION %d ----\n", k)); sink() })
  ok
}

# ------------------------------------------------------------------ section 1: checks, cohort, controls ------------------------------------------------------------------
sec1 <- function() {
  t00 <- Sys.time(); cat("SqnPAN_12 section 1 (checks, cohort, controls), started", format(t00), "\n")
  cat(R.version.string, "| survival", as.character(utils::packageVersion("survival")), "| readxl", as.character(utils::packageVersion("readxl")), "\n")
  cat("coxph arguments:", paste(names(formals(survival::coxph)), collapse = ", "), "\n")
  cat("survreg arguments:", paste(names(formals(survival::survreg)), collapse = ", "), "\n")
  cat("Fixed: ties = efron; strata(type); survreg dist = weibull; TOST margin HR", EQ_HR2, "; helpers md5", MD5_HELPERS2, "\n\n")
  if (min(diff(sort(unname(SEED2)))) < 100000) stop("seed bases are less than 100,000 apart")
  mm <- unname(tools::md5sum(F2))
  print(data.frame(file = names(F2), md5 = mm, expected = MD5_2, ok = mm == MD5_2), row.names = FALSE)
  if (any(mm != MD5_2)) stop("md5 mismatch in an input file")

  # mutations
  hd <- strsplit(readLines(gzfile(F2["mut"]), n = 1L), "\t", fixed = TRUE)[[1]]
  exp_hd <- c("case_id", "submitter", "project", "accession", "gene", "position", "ref", "alt", "n_rows_all_files", "n_rows_chosen_file", "in_chosen_file", "only_rescued",
              "canonical_any", "hypermutated", "res_role", "lost_sequon_N", "lost_tiers", "gained_sequon_N", "gain_states", "gain_state_best", "class")
  if (!identical(hd, exp_hd)) stop("mutation table header differs from the expected 21 columns")
  mcls <- c("character", "NULL", "NULL", "character", "NULL", "character", "character", "character", "NULL", "NULL", "integer", "integer", "NULL", "NULL", "character",
            "NULL", "character", "NULL", "NULL", "NULL", "character")
  t0 <- Sys.time()
  mut <- utils::read.delim(gzfile(F2["mut"]), sep = "\t", header = TRUE, quote = "", comment.char = "", colClasses = mcls, na.strings = "NA", check.names = FALSE)
  cat("mutation rows read:", fmt2(nrow(mut)), "in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
  if (nrow(mut) != 1467936L) stop("mutation table has ", nrow(mut), " rows, expected 1,467,936")
  if (anyNA(mut$in_chosen_file) || anyNA(mut$only_rescued)) stop("missing in_chosen_file or only_rescued")
  mut <- mut[mut$in_chosen_file == 1L, , drop = FALSE]
  if (nrow(mut) != 1460792L) stop("chosen-file rows ", nrow(mut), ", expected 1,460,792")
  isloss <- mut$class %in% c("loss", "loss_and_gain") & grepl("(^|;)T1(;|$)", mut$lost_tiers)
  loss <- mut[isloss, , drop = FALSE]
  loss$ltype <- ifelse(loss$ref == "N", "N", ifelse(loss$ref %in% c("S", "T"), "ST", "odd"))
  cat("T1 loss rows (chosen files):", nrow(loss), "(SqnPAN_09 run: 2,757) in", length(unique(loss$case_id)), "cases (SqnPAN_09 run: 1,468)\n")
  if (length(unique(loss$case_id)) != 1468L) stop("T1 loss cases differ from the 1,468 of the SqnPAN_09 run")
  print(table(loss$ltype))
  isctrl <- mut$res_role %in% c("nonseqN", "nonseqST") & !(mut$res_role == "nonseqST" & mut$alt %in% c("S", "T"))
  ctrl <- mut[isctrl, , drop = FALSE]
  cat("background (control-type) rows, S<->T removed:", fmt2(nrow(ctrl)), "\n")
  rm(mut)

  # cases and CDR
  cases <- rd2(F2["cases"]); need2(cases, c("case_id", "submitter", "project", "missense_rows_chosen", "hypermutated"), "cases table")
  if (nrow(cases) != 10190L || anyDuplicated(cases$submitter) || anyDuplicated(cases$case_id)) stop("cases table: wrong row count or duplicated ids")
  cases$type <- sub("^TCGA-", "", cases$project); cases$missense <- as.integer(cases$missense_rows_chosen); cases$hyper <- as.integer(cases$hypermutated)
  cdr <- as.data.frame(readxl::read_excel(F2["cdr"], sheet = "TCGA-CDR", col_types = "text", na = c("#N/A", "[Not Available]", "[Not Applicable]", "[Discrepancy]", "")))
  need2(cdr, c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis", "gender", "OS", "OS.time", "PFI", "PFI.time"), "CDR sheet")
  if (anyDuplicated(cdr$bcr_patient_barcode)) stop("duplicated patient barcode in the CDR")
  j <- match(cases$submitter, cdr$bcr_patient_barcode); cases$in_cdr <- !is.na(j)
  cases$age <- suppressWarnings(as.numeric(cdr$age_at_initial_pathologic_diagnosis[j])); cases$sex <- factor(cdr$gender[j], levels = c("FEMALE", "MALE"))
  for (v in c("PFI", "OS")) cases[[v]] <- suppressWarnings(as.integer(cdr[[v]][j]))
  for (v in c("PFI.time", "OS.time")) cases[[v]] <- suppressWarnings(as.numeric(cdr[[v]][j]))
  cat("cases:", nrow(cases), "; in the CDR:", sum(cases$in_cdr), "; project type differs from CDR type in", sum(cdr$type[j] != cases$type, na.rm = TRUE), "cases\n")
  E2$cases <- cases; E2$loss <- loss; E2$ctrl <- ctrl

  base <- cases[cases$type %in% TYPES29 & cases$missense >= 1L, , drop = FALSE]
  E2$protset <- unique(loss$accession[loss$case_id %in% base$case_id])
  E2$protN <- unique(loss$accession[loss$case_id %in% base$case_id & loss$ltype == "N"]); E2$protST <- unique(loss$accession[loss$case_id %in% base$case_id & loss$ltype == "ST"])
  cat("loss protein set (29 types):", length(E2$protset), "proteins;", length(E2$protN), "with an N loss,", length(E2$protST), "with an S/T loss\n")
  cat("left out: LAML cases", sum(cases$type == "LAML"), ", PCPG cases", sum(cases$type == "PCPG"), "\n\n")

  bc <- build_cohort(TYPES29, "PFI"); cat("Cohort steps (PFI, 29 types):\n"); print(bc$steps, row.names = FALSE)
  A0 <- add_dummies(assign_exp(bc$d)); E2$A0 <- A0
  cc3 <- A0[A0$grp != "neither", , drop = FALSE]; cc3$grp <- as.character(cc3$grp); E2$cc3 <- cc3; E2$cc4 <- A0
  cat("\nExposure groups in the analysis cohort (sizes only):\n"); print(table(A0$grp))
  gc <- as.data.frame.matrix(table(A0$type, factor(A0$grp, levels = c("control_only", "loss_only", "both", "neither"))))
  gc$type <- rownames(gc); gc$n_cases <- rowSums(gc[, 1:4]); gc$events_in_type <- as.integer(tapply(A0$event, A0$type, sum)[gc$type]); gc <- gc[, c("type", "control_only", "loss_only", "both", "neither", "n_cases", "events_in_type")]
  cat("md5 group counts:", wr2(gc, file.path(TBL2, "P2_12_group_counts.tsv")), "\n"); pr2(gc)
  cat("hypermutated cases in cohort:", sum(A0$hyper == 1L), "\n\n")

  # ---- controls ----
  cat("=== CONTROLS (judged before any real exposure HR is read) ===\n")
  pc <- cox2(add_dummies(cc3), "age10", c("sex", "logmis", "x_lo", "x_both"))
  posok <- isTRUE(pc$HR > 1 && pc$p / 2 < 0.025)
  cat(sprintf("PC pipeline control: age HR per 10 years %.3f (95%% CI %.3f to %.3f), p %.3g -> %s\n\n", pc$HR, pc$HR_lo, pc$HR_hi, pc$p, if (posok) "PASSED" else "FAILED"))
  SF <- lapply(SPECS, function(sp) sim_prep2(if (sp$data == "cc3") cc3 else A0, sp$xs, sp$plant))
  for (si in seq_along(SPECS)) cat(sprintf("spec %s: data %s, fitted HR of the planted exposure %.3f, censoring model %s, Weibull shape sigma %.3f\n", SPECNAMES[si], SPECS[[si]]$data,
                                          SF[[si]]$hr_fitted, SF[[si]]$cens_kind, SF[[si]]$sigma))
  tP <- unname(system.time(for (i in 1:5) sim_one2(SF[[1]], "x_lo", log(1.5), 1L))["elapsed"]) / 5
  tM <- unname(system.time(for (i in 1:5) sim_one2(SF[[4]], "x_L", log(1.5), 1L))["elapsed"]) / 5
  est <- 2 * ((500 * 3 + 300 * 3 * 2) * tP + 300 * 3 * tM + NNULL2 * tP) / 60
  cat(sprintf("Timing from 5 measured fits each: %.3f s (cc3) and %.3f s (cc4) per simulated fit; estimated total of the power and null runs %.1f minutes (an estimate, not a measurement)\n\n", tP, tM, est))

  PW <- list(); NR <- list()
  for (set in c("A", "B")) {
    b <- unname(SEED2[paste0("pow", set)])
    for (si in seq_along(SPECS)) {
      sp <- SPECS[[si]]; t1 <- Sys.time()
      pw <- run_power2(NULL, sp, si, sp$nrep, HRS2, b, sf = SF[[si]]); pw$set <- set; PW[[length(PW) + 1L]] <- pw
      cat(sprintf("power set %s, %s done in %.1f min\n", set, sp$name, as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    }
    t1 <- Sys.time(); nr <- run_null2(cc3, unname(SEED2[paste0("null", set)]), NNULL2); nr$set <- set; NR[[set]] <- nr
    cat(sprintf("null set %s done in %.1f min\n", set, as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  }
  PW <- do.call(rbind, PW); NRall <- do.call(rbind, NR)
  cat("\nPower table (power_05 = p < 0.05 and HR > 1):\n"); pr2(PW[, c("set", "spec", "HR_planted", "nrep", "n_failed", "power_05", "false_pos_two_sided", "equiv_share", "median_events_group", "median_HR")])
  nullrate <- tapply(NRall$p < 0.05, NRall$set, function(v) mean(v, na.rm = TRUE))
  cat(sprintf("\nNL null: share with p < 0.05 for x_lo: set A %.3f, set B %.3f (rule: set A <= %.2f)\n", nullrate["A"], nullrate["B"], FWER_MAX2))
  nullok <- isTRUE(unname(nullrate["A"]) <= FWER_MAX2)
  PWA <- PW[PW$set == "A", ]
  sens <- vapply(SPECNAMES, function(n) isTRUE(PWA$power_05[PWA$spec == n & PWA$HR_planted == 1.5] >= POWER_MIN2), TRUE)
  e1ok <- vapply(SPECNAMES, function(n) isTRUE(PWA$equiv_share[PWA$spec == n & PWA$HR_planted == 1.5] <= E1_MAX2), TRUE)
  cat("\nSensitive (power >= 0.80 at HR 1.5, set A):\n"); print(sens); cat("E1 ok (Equivalent share at HR 1.5 <= 0.05, set A):\n"); print(e1ok)
  # determinism
  r1 <- do.call(rbind, lapply(1:10, function(i) sim_one2(SF[[1]], "x_lo", log(1.5), unname(SEED2["powA"]) + 20000 + 2000 * 3 + i)))
  r2 <- do.call(rbind, lapply(1:10, function(i) sim_one2(SF[[1]], "x_lo", log(1.5), unname(SEED2["powA"]) + 20000 + 2000 * 3 + i)))
  detok <- identical(cmd5_2(r1), cmd5_2(r2)); cat("Determinism (10 repeats recomputed):", if (detok) "identical" else "DIFFERENT", "\n")
  E2$V <- list(sens = sens, e1ok = e1ok, nullok = nullok, posok = posok, detok = detok, nullA = unname(nullrate["A"]), nullB = unname(nullrate["B"]))
  cat("md5 power table:", wr2(PW, file.path(TBL2, "P2_12_power_run1.tsv")), "\n")
  cat("md5 null table:", wr2(NRall, file.path(TBL2, "P2_12_null_run1.tsv")), "\n")
  vl <- c(paste("P2_12 controls verdict", format(Sys.time())), paste("pipeline_age_control_passed", posok), paste("null_check_passed", nullok), paste("determinism", detok),
          paste0("sensitive_", SPECNAMES, "=", sens), paste0("e1ok_", SPECNAMES, "=", e1ok))
  writeLines(vl, FVER2); cat("verdict written to", FVER2, "\n")
  if (!detok) stop("determinism check failed")
  cat("section 1 finished in", round(as.numeric(difftime(Sys.time(), t00, units = "mins")), 1), "minutes\n")
}

# ------------------------------------------------------------------ section 2: real results ------------------------------------------------------------------
compute_results <- function() {
  mk <- function(label, d, covs = COV_MAIN2, zph1 = FALSE) {
    d <- add_dummies(d); c3 <- d[d$grp != "neither", , drop = FALSE]
    out <- rbind(cbind(model = "M1_3level", cox2(c3, c("x_lo", "x_both"), covs, zph = zph1)), cbind(model = "M2_any_loss", cox2(c3, "x_any", covs)),
                 cbind(model = "M3_adj_control", cox2(d, c("x_L", "x_C"), covs)))
    cbind(analysis = label, out)
  }
  A0 <- assign_exp(build_cohort(TYPES29, "PFI")$d)
  sets <- list(S1_strict_controls = assign_exp(build_cohort(TYPES29, "PFI")$d, strict = TRUE), S2_no_hypermutated = A0[A0$hyper == 0L, , drop = FALSE],
               S3_no_rescued_rows = assign_exp(build_cohort(TYPES29, "PFI")$d, norescue = TRUE), S4_N_loss_only = A0[!(A0$LST | A0$LO), , drop = FALSE],
               S5_ST_loss_only = A0[!(A0$LN | A0$LO), , drop = FALSE], S6_OS_endpoint = assign_exp(build_cohort(TYPES29, "OS")$d),
               S7_plus_DLBC_KICH = assign_exp(build_cohort(c(TYPES29, "DLBC", "KICH"), "PFI")$d))
  long <- mk("main", A0, zph1 = TRUE)
  for (nm in names(sets)) long <- rbind(long, mk(nm, sets[[nm]]))
  long <- rbind(long, mk("S8_unadjusted", A0, covs = character(0)))
  types <- sort(unique(as.character(A0$type)))
  loo <- do.call(rbind, lapply(types, function(t) {
    d <- add_dummies(A0[A0$type != t, , drop = FALSE]); c3 <- d[d$grp != "neither", , drop = FALSE]
    cbind(left_out = t, rbind(cbind(model = "M1_3level", cox2(c3, c("x_lo", "x_both"), COV_MAIN2)), cbind(model = "M2_any_loss", cox2(c3, "x_any", COV_MAIN2))))
  }))
  pt <- do.call(rbind, lapply(types, function(t) {
    d <- add_dummies(A0[A0$type == t, , drop = FALSE]); c3 <- d[d$grp != "neither", , drop = FALSE]
    n_any <- sum(c3$x_any == 1); ev_any <- sum(c3$event[c3$x_any == 1])
    r <- data.frame(type = t, n_control_only = sum(c3$grp == "control_only"), n_any_loss = n_any, events_any_loss = ev_any, tested = (n_any >= 15 & ev_any >= 10),
                    HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_, p = NA_real_, stringsAsFactors = FALSE)
    if (r$tested) { z <- cox2(c3, "x_any", COV_MAIN2); r$HR <- z$HR; r$HR_lo <- z$HR_lo; r$HR_hi <- z$HR_hi; r$p <- z$p }
    r
  }))
  pt$p_bh <- NA_real_; pt$p_bh[pt$tested] <- stats::p.adjust(pt$p[pt$tested], "BH")
  list(long = long, loo = loo, pt = pt, A0 = A0)
}
sec2 <- function() {
  t00 <- Sys.time(); cat("SqnPAN_12 section 2 (real results), started", format(t00), "\n"); V <- E2$V
  cat("Controls in force: pipeline", V$posok, "| null", V$nullok, "(A", round(V$nullA, 3), ", B", round(V$nullB, 3), ") | sensitive:", paste(names(V$sens), V$sens, collapse = "; "), "\n\n")
  R1 <- compute_results(); R2 <- compute_results()
  same <- c(long = identical(cmd5_2(R1$long), cmd5_2(R2$long)), loo = identical(cmd5_2(R1$loo), cmd5_2(R2$loo)), pt = identical(cmd5_2(R1$pt), cmd5_2(R2$pt)))
  cat("Computed twice, identical content:", paste(names(same), same, collapse = "; "), "\n"); if (!all(same)) stop("two computations differ")
  L <- R1$long; A0 <- R1$A0
  gt <- do.call(rbind, lapply(c("control_only", "loss_only", "both", "neither"), function(g) { s <- A0[A0$grp == g, ]
    data.frame(group = g, n = nrow(s), events = sum(s$event), hypermutated = sum(s$hyper == 1L), median_log_missense = round(stats::median(s$logmis), 2), stringsAsFactors = FALSE) }))
  cat("\nGroups (analysis cohort):\n"); pr2(gt)
  prim <- L[L$analysis == "main" & L$model == "M1_3level" & L$term == "x_lo", , drop = FALSE]
  prim$spec <- SPECNAMES[1]; prim$family <- "primary"; prim$p_adj <- prim$p; prim$tost_adj <- prim$p_tost
  sec <- rbind(L[L$analysis == "main" & L$model == "M1_3level" & L$term == "x_both", ], L[L$analysis == "main" & L$model == "M2_any_loss" & L$term == "x_any", ],
               L[L$analysis == "main" & L$model == "M3_adj_control" & L$term == "x_L", ])
  sec$spec <- SPECNAMES[2:4]; sec$family <- "secondary_holm3"; sec$p_adj <- stats::p.adjust(sec$p, "holm"); sec$tost_adj <- stats::p.adjust(sec$p_tost, "holm")
  ps <- rbind(prim, sec); ps$label_raw <- label2(ps$p_adj, ps$tost_adj, ps$HR)
  ps$label <- final2(ps$label_raw, unname(V$sens[ps$spec]), unname(V$e1ok[ps$spec]), V$nullok, V$posok)
  cat("\n=== PRIMARY and SECONDARY (main analysis, PFI, 29 types) ===\n")
  pr2(ps[, c("family", "spec", "term", "n", "events", "n_group", "events_group", "HR", "HR_lo", "HR_hi", "p", "p_adj", "p_tost", "tost_adj", "label_raw", "label", "zph_term", "zph_global")])
  cat("\n=== All models (reported, not judged) ===\n")
  pr2(L[L$term %in% c("x_lo", "x_both", "x_any", "x_L", "x_C"), c("analysis", "model", "term", "n", "events", "n_group", "events_group", "HR", "HR_lo", "HR_hi", "p", "note")])
  cat("\n=== Leave-one-cancer-type-out: range of HR ===\n")
  for (m in c("M1_3level", "M2_any_loss")) for (tm in unique(R1$loo$term[R1$loo$model == m])) {
    s <- R1$loo[R1$loo$model == m & R1$loo$term == tm, ]
    cat(sprintf("%s %s: HR %.2f to %.2f over %d fits (left out: min at %s, max at %s)\n", m, tm, min(s$HR, na.rm = TRUE), max(s$HR, na.rm = TRUE), nrow(s), s$left_out[which.min(s$HR)], s$left_out[which.max(s$HR)]))
  }
  cat("\n=== Per cancer type: any loss against control-only (descriptive; tested if >= 15 cases and >= 10 events) ===\n"); pr2(R1$pt)
  cat("\nmd5 results_long:", wr2(L, file.path(TBL2, "P2_12_results_long.tsv")), "\n")
  cat("md5 primary_secondary:", wr2(ps, file.path(TBL2, "P2_12_primary_secondary.tsv")), "\n")
  cat("md5 leave_one_type_out:", wr2(R1$loo, file.path(TBL2, "P2_12_leave_one_type_out.tsv")), "\n")
  cat("md5 per_type:", wr2(R1$pt, file.path(TBL2, "P2_12_per_type.tsv")), "\n")
  cat("md5 groups:", wr2(gt, file.path(TBL2, "P2_12_groups.tsv")), "\n")
  co <- A0[, c("submitter", "type", "grp", "nL", "nC", "LN", "LST", "LO", "hyper", "age", "sex", "logmis", "PFI", "PFI.time", "OS", "OS.time")]
  cat("md5 cohort:", wr2(co, file.path(TBL2, "P2_12_cohort.tsv")), "\n")
  cat("section 2 finished in", round(as.numeric(difftime(Sys.time(), t00, units = "mins")), 1), "minutes\n")
}

# ------------------------------------------------------------------ section 3: driver genes (descriptive) ------------------------------------------------------------------
sec3 <- function() {
  t00 <- Sys.time(); cat("SqnPAN_12 section 3 (driver gene overlap, descriptive), started", format(t00), "\n")
  prot <- rd2(F2["prot"]); need2(prot, c("accession", "gene", "n_T1_sequons"), "proteins table")
  if (anyDuplicated(prot$accession)) stop("duplicated accession in the proteins table")
  prot$n_T1_sequons <- as.integer(prot$n_T1_sequons)
  loss <- E2$loss; loss$gene <- prot$gene[match(loss$accession, prot$accession)]; loss$type <- E2$cases$type[match(loss$case_id, E2$cases$case_id)]
  cat("T1 loss rows (all 33 types, chosen files):", nrow(loss), "; without a gene symbol:", sum(is.na(loss$gene)), "\n")
  b <- as.data.frame(readxl::read_excel(F2["bailey"], sheet = "Table S1", skip = 3, col_types = "text"))
  need2(b, c("Gene", "Cancer", "Decision"), "Bailey Table S1"); tsn <- grep("^Tumor suppressor", names(b), value = TRUE)[1]
  b <- b[!is.na(b$Gene) & nzchar(b$Gene), , drop = FALSE]
  cat("Bailey Table S1 rows with a gene:", nrow(b), "(expected 739); unique genes:", length(unique(b$Gene)), "(expected 299)\n")
  if (nrow(b) != 739L || length(unique(b$Gene)) != 299L) stop("Bailey Table S1 counts differ from the recorded counts")
  ALIAS <- data.frame(old = c("FAM46D", "H3F3A", "H3F3C", "HIST1H1C", "HIST1H1E", "RQCD1", "TCEB1", "WHSC1"), new = c("TENT5D", "H3-3A", "H3-5", "H1-2", "H1-4", "CNOT9", "ELOC", "NSD2"),
                      acc = c("Q8NEK8", "P84243", "Q6NXT2", "P16403", "P10412", "Q92600", "Q15369", "O96028"), stringsAsFactors = FALSE)   # [from memory]
  g <- prot$gene[match(ALIAS$acc, prot$accession)]; ALIAS$guard_ok <- !is.na(g) & g == ALIAS$new; ALIAS$gene_in_table <- g
  cat("Alias step (guard = the accession in P2_06_proteins.tsv carries the new symbol):\n"); pr2(ALIAS)
  b$gene_cur <- b$Gene; for (i in which(ALIAS$guard_ok)) b$gene_cur[b$Gene == ALIAS$old[i]] <- ALIAS$new[i]
  drv <- unique(b$gene_cur); unm <- drv[!(drv %in% prot$gene)]
  cat("Driver genes after the alias step:", length(drv), "; found among reviewed proteins:", length(drv) - length(unm), "; not found:", paste(unm, collapse = ", "), "\n\n")
  lg <- unique(loss$gene[!is.na(loss$gene)]); ov <- sort(intersect(lg, drv))
  tab <- do.call(rbind, lapply(ov, function(gn) {
    ls <- loss[loss$gene == gn & !is.na(loss$gene), , drop = FALSE]; bb <- b[b$gene_cur == gn, , drop = FALSE]
    canc <- sort(unique(bb$Cancer)); bt <- unique(unlist(lapply(canc, function(x) if (x == "COADREAD") c("COAD", "READ") else x)))
    lt <- table(ls$type); vv <- sort(table(paste0(ls$ref, ls$position, ls$alt)), decreasing = TRUE); nv <- min(10L, length(vv)); role <- bb[[tsn]]; role <- sort(unique(role[!is.na(role)]))
    data.frame(gene = gn, accession = paste(unique(ls$accession), collapse = ";"), n_T1_sequons = prot$n_T1_sequons[match(ls$accession[1], prot$accession)],
               bailey_cancers = paste(canc, collapse = ","), bailey_role = paste(role, collapse = ","), n_loss_rows = nrow(ls), n_loss_cases = length(unique(ls$case_id)), n_types = length(lt),
               types = paste0(names(lt), ":", as.integer(lt), collapse = ";"), loss_in_listed_type = ("PANCAN" %in% bt) | any(names(lt) %in% bt),
               variants = paste0(names(vv)[seq_len(nv)], "x", as.integer(vv)[seq_len(nv)], collapse = ";"), stringsAsFactors = FALSE)
  }))
  bg <- unique(prot$gene[!is.na(prot$gene)]); t1g <- unique(prot$gene[!is.na(prot$gene) & prot$n_T1_sequons > 0L]); drv_in <- intersect(drv, bg)
  sm <- data.frame(set = c("genes with a T1 loss mutation", "genes with a T1 sequon", "all reviewed genes"), n_genes = c(length(lg), length(t1g), length(bg)),
                   n_drivers = c(length(intersect(lg, drv_in)), length(intersect(t1g, drv_in)), length(drv_in)), stringsAsFactors = FALSE)
  sm$share_drivers <- round(sm$n_drivers / sm$n_genes, 4)
  cat("Summary (context only, no test):\n"); pr2(sm)
  cat("\nDriver genes with a T1 loss mutation:", nrow(tab), "; with a loss in a cancer type where Bailey lists the gene:", sum(tab$loss_in_listed_type), "; loss rows in these genes:", sum(tab$n_loss_rows), "\n")
  pr2(tab[order(-tab$n_loss_cases), c("gene", "bailey_cancers", "bailey_role", "n_loss_rows", "n_loss_cases", "n_types", "loss_in_listed_type", "variants")])
  cat("\nmd5 driver_overlap:", wr2(tab, file.path(TBL2, "P2_12_driver_overlap.tsv")), "\n")
  cat("md5 driver_summary:", wr2(sm, file.path(TBL2, "P2_12_driver_summary.tsv")), "\n")
  cat("md5 alias_step:", wr2(ALIAS, file.path(TBL2, "P2_12_driver_alias_step.tsv")), "\n")
  cat("section 3 finished in", round(as.numeric(difftime(Sys.time(), t00, units = "mins")), 1), "minutes\n")
}

ok1 <- run_section(1, "checks_controls", sec1); if (!ok1) stop("Section 1 failed; see results/logs/P2_12_section1_checks_controls.txt")
ok2 <- run_section(2, "results", sec2); if (!ok2) stop("Section 2 failed; see results/logs/P2_12_section2_results.txt")
ok3 <- run_section(3, "drivers", sec3); if (!ok3) stop("Section 3 failed; see results/logs/P2_12_section3_drivers.txt")
cat("SqnPAN_12 finished: sections 1 to 3 completed.\n")
