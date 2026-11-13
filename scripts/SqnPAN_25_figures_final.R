## ===================================================================================================================
## SqnPAN_25_figures_final.R - figures as submitted (8 main and 17 supplementary figures)
## Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
## Merged main figures, supplementary figures from every remaining panel, and every
## working term replaced by plain wording. Run from the repository folder:
##   Rscript scripts/SqnPAN_25_figures_final.R final
##   draft -> results\figures\P2_25_draft (PNG only); final -> results\figures\P2_25_final\pdf, \tiff (600 dpi LZW), \png (300 dpi)
##   Only some figures: add file names, e.g.  ... final Figure_2 Figure_S3
## Reads SqnPAN_24_figures.R (not changed) for theme, palette and panel code, and the same tables. Changes no table or label.
## Text control: the working terms are replaced in the panel code before the panels are built; then every text element of
##   every panel (titles, subtitles, captions, axis titles and tick labels, strips, legends, text layers) is collected,
##   written to P2_25_text_scan.tsv and checked against a list of working-term patterns; hits are logged.
## Layout fixes (2026-11-12): free() releases Figure 1A, 1C, 4A, 4C, 5A, 6D, 6E, 7E, 8A and 8C from the column
##   alignment that pushed them right; Figure 4B wider with overlapping ticks dropped; Figure 7A title starts at the plot edge.
##   Second pass: Figure 3 built from its panels (3A, 3B released, 3C labels wrapped); Figure 6A released and 6E kept aligned
##   (two free() panels in one row gave an error); Figure 8B labels wrapped so its letter sits by its panel.
##   Third pass: every panel letter on the same line as its panel title (TAGTH), all 25 figures.
##   Panel content unchanged.
## Seeds: no new random numbers (the SqnPAN_24 seeds 33700001 and 33700002 act inside the panel code).
## ===================================================================================================================

t_start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
MODE <- "draft"; ONLY <- character(0)
if (length(args) >= 1) { if (tolower(args[1]) %in% c("draft", "final")) MODE <- tolower(args[1]) else stop("The first argument must be draft or final") }
if (length(args) >= 2) ONLY <- args[-1]
ROOT <- Sys.getenv("SQNPAN_ROOT", unset = local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f) == 1) normalizePath(file.path(dirname(normalizePath(f)), ".."), winslash = "/") else normalizePath(getwd(), winslash = "/") }))
TAB <- file.path(ROOT, "results", "tables"); PRO <- file.path(ROOT, "data", "processed")
SRC <- file.path(ROOT, "scripts", "SqnPAN_24_figures.R")
OUT <- file.path(ROOT, "results", "figures", if (MODE == "draft") "P2_25_draft" else "P2_25_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (MODE == "final") for (s in c("pdf", "tiff", "png")) dir.create(file.path(OUT, s), showWarnings = FALSE)
LOGF <- file.path(OUT, "P2_25_log.txt"); cat("", file = LOGF)
lg <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(...)); message(m); cat(m, "\n", sep = "", file = LOGF, append = TRUE) }
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
lg("SqnPAN_25 start, mode ", MODE, if (length(ONLY)) paste0(", only ", paste(ONLY, collapse = " ")) else "", "; R ", getRversion(),
   ", ggplot2 ", as.character(packageVersion("ggplot2")), ", patchwork ", as.character(packageVersion("patchwork")))
lg("SqnPAN_24_figures.R md5 ", unname(tools::md5sum(SRC)))
FONT <- "Times New Roman"
fp <- tryCatch(systemfonts::match_fonts(FONT)$path, error = function(e) tryCatch(systemfonts::match_font(FONT)$path, error = function(e2) NA_character_))
lg("Font ", FONT, " -> ", fp)

## ---------------- panel code from SqnPAN_24, with the working terms replaced ----------------
src <- readLines(SRC, warn = FALSE)
i0 <- grep("^## ---------------- PALETTE", src)[1]; i1 <- which(src == "## Output") - 2
if (is.na(i0) || length(i1) != 1 || i1 <= i0) stop("could not find the panel code block in SqnPAN_24_figures.R")
code <- paste(src[i0:i1], collapse = "\n")
SUBS <- list(
  ## Corrections (2026-11-12): Figure 3B without the empty disulfide His rows, Figure 8C null share to one decimal, Figure S6C rule text as in SM7
  c('c1 <- cs[cs$config == "main" & cs$scope == "pooled" & !is.na(cs$est), ]', 'c1 <- cs[cs$config == "main" & cs$scope == "pooled" & !is.na(cs$est) & !(cs$tier == "T3" & cs$group != "N"), ]'),
  c('pct(n0["main"], 0), pct(n0["nohyper"], 0)', 'pct(n0["main"], 1), pct(n0["nohyper"], 1)'),
  c('" 16 of 20;\\nno Changed call in any null set"', '" 16 of 20 Equivalent;\\nat most 1 of 20 Changed"'),
  ## tier names
  c('TIER <- c(T1 = "T1 confirmed", T2 = "T2 predicted", T3 = "T3 motif only")',
    'TIER <- c(T1 = "Tier 1", T2 = "Tier 2", T3 = "Tier 3"); TIER_LONG <- c(T1 = "Tier 1 (experimental)", T2 = "Tier 2 (annotated)", T3 = "Tier 3 (motif only)")'),
  c('tc$lab <- factor(TIER[tc$tier], levels = rev(TIER))', 'tc$lab <- factor(TIER_LONG[tc$tier], levels = rev(TIER_LONG))'),
  c('loss_T1 = "T1 confirmed", loss_T2 = "T2 predicted", loss_T3 = "T3 motif only"',
    'loss_T1 = "Tier 1 (experimental)", loss_T2 = "Tier 2 (annotated)", loss_T3 = "Tier 3 (motif only)"'),
  c('"T1: experimental UniProt', '"Tier 1: experimental UniProt'), c('T2: other UniProt', 'Tier 2: other UniProt'), c('; T3: motif only"', '; Tier 3: motif only"'),
  c('short <- c(T1 = "T1", T2 = "T2", T3 = "T3")', 'short <- c(T1 = "Tier 1", T2 = "Tier 2", T3 = "Tier 3")'),
  c('"Planted effects, T1"', '"Planted effects, Tier 1"'), c('"dNdScv (T1 only)"', '"dNdScv (Tier 1 only)"'),
  ## F01 boxes and notes
  c('completes a new N-X-S/T (B3)"', 'completes a new N-X-S/T"'),
  c('"Q1 Catalogue of', '"Catalogue of'), c('"Q2 Selection test:', '"Selection test:'), c('"Q3 Structure: FoldX "', '"Structure: FoldX "'),
  c('"Q4 Clinical links:', '"Clinical links:'), c('"A1 Context-aware', '"Context-aware'), c('"A2 Germline control;', '"Germline control;'),
  c('"A3 Mutational', '"Mutational'), c('"B1 AlphaMissense', '"AlphaMissense'), c('"B2 CPTAC protein', '"CPTAC protein'),
  c('"B3 Gain side', '"Gain side'), c('"B4 Hotspots and', '"Hotspots and'), c('"C1 Aromatic', '"Aromatic'),
  c('"Added analyses are post hoc (Amendment 4 and dNdScv), each with its own controls; no primary label changed"',
    '"Added analyses are post hoc, each with its own controls; no primary label changed"'),
  ## selection test, per type, germline
  c('M1 <- "Position-count model (Q2, primary)"; M2 <- "Context-aware model (A1)"', 'M1 <- "Position-count model (primary)"; M2 <- "Context-aware model"'),
  c('"Q2 position-count"', '"Position-count"'), c('"A1 context-aware"', '"Context-aware"'),
  c('title = "Calibration of the expected rate (A1)"', 'title = "Calibration of the expected rate"'),
  c('"context-aware model (A1)"', '"context-aware model"'), c('"position-count model (Q2)"', '"position-count model"'),
  c('S1 <- "Tumour somatic, glycosites (A1)"; S2 <- "Germline, glycosites (A2)"', 'S1 <- "Tumour somatic, glycosites"; S2 <- "Germline, glycosites"'),
  c('so the tumour Equivalent reads Equivalent_biology_control_FAILED"', 'so the tumour Equivalent carries the note germline biology control failed"'),
  c('"Tumour: chosen-file SNVs per channel (counts pooled over types; the model uses per-type rates); germline: option A variant set"',
    '"Tumour: SNVs per channel from one selected file per tumour (counts pooled over types; the model uses per-type rates); germline: primary variant set"'),
  ## clinical
  c('the test is the stratified Cox model in A"', 'the test is the stratified Cox model"'),
  c('"S1 strict controls"', '"Strict controls"'), c('"S2 without hypermutated"', '"Without hypermutated"'), c('"S3 without rescued rows"', '"Without rescued rows"'),
  c('"S4 N losses only"', '"N losses only"'), c('"S5 S/T losses only"', '"S/T losses only"'), c('"S6 overall survival"', '"Overall survival"'),
  c('"S7 plus DLBC and KICH"', '"Plus DLBC and KICH"'), c('"S8 unadjusted"', '"Unadjusted"'),
  ## hotspots, ActiveDriver, CPTAC
  c('(label Hotspot_calibration_FAILED)"', '(label Hotspot, calibration failed)"'), c('"Bailey driver gene"', '"Known driver gene"'),
  c('"Agreement with B4 hotspots"', '"Agreement with single-position hotspots"'),
  c('x = "LOSS carrier-protein pairs"', 'x = "Loss carrier-protein pairs"'), c('"%s LOSS pairs in CPTAC tumours"', '"%s loss carrier-protein pairs in CPTAC tumours"'),
  c('Primary label Equivalent_positive_control_FAILED:', 'Primary label Equivalent, positive control failed;'),
  ## synthesis
  c('"Q2 position-count model (primary), T1 losses"', '"Position-count model (primary), Tier 1 losses"'),
  c('"A1 context-aware model, T1 losses"', '"Context-aware model, Tier 1 losses"'), c('"B3 sequon gains, plausible"', '"Sequon gains, plausible regions"'),
  c('"A2 germline, T1 glycosites (depletion expected)"', '"Germline, Tier 1 glycosites (depletion expected)"'),
  c('paste0("Q3 FoldX ", U_DD', 'paste0("FoldX ", U_DD'), c('paste0("B1 ", c(AM = "AlphaMissense", ESM = "ESM1b")[k]', 'paste0(c(AM = "AlphaMissense", ESM = "ESM1b")[k]'),
  c('paste0("Q4 PFI, "', 'paste0("Progression-free interval, "'), c('(labels as written by the analysis scripts)"', '(labels as pre-specified)"'),
  c('s <- do.call(rbind, rows); s$label <- vclass(s$lab);',
    's <- do.call(rbind, rows); if (exists("F22_BLOCKS", envir = globalenv())) s <- s[s$block %in% get("F22_BLOCKS", envir = globalenv()), , drop = FALSE]; s$label <- vclass(s$lab);'),
  ## label and capability maps
  c('"Q2 position-count model (tumour)"', '"Position-count model (tumour)"'), c('"A1 context-aware model (tumour)"', '"Context-aware model (tumour)"'),
  c('"A2 germline glycosites"', '"Germline glycosites"'), c('"Germline constrained residues (N = Cys, S/T = His)"', '"Germline Cys (N) and His (S/T); tiers: metal, active site, disulfide"'),
  c('"B3 sequon gains (T1 to T3 = plausibility)"', '"Sequon gains (tiers: plausible, unknown, implausible)"'), c('"B3 sequon gains"', '"Sequon gains"'),
  c('uniform_restricted = "Position counts, A1 protein set"', 'uniform_restricted = "Position counts, context-model protein set"'),
  ## layout: wider label wrap for the Fan et al. counts
  c('scale_y_discrete(labels = scales::label_wrap(40))', 'scale_y_discrete(labels = scales::label_wrap(75))')
)
for (s in SUBS) {
  n <- lengths(regmatches(code, gregexpr(s[1], code, fixed = TRUE)))
  if (n == 0) stop("text to replace not found in SqnPAN_24_figures.R: ", s[1])
  code <- gsub(s[1], s[2], code, fixed = TRUE)
}
lg(length(SUBS), " text replacements applied")
## remaining tier codes inside display strings ("T1 losses" -> "Tier 1 losses"); code such as "T2 = " is left alone
code <- gsub('(["\\s])T([123])(\\s)(?!=)', '\\1Tier \\2\\3', code, perl = TRUE)
eval(parse(text = code), envir = globalenv())
lg(length(FIGS), " panel builders loaded from SqnPAN_24")

panels <- function(id) { f <- FIGS[[id]]$fun; b <- body(f)
  b[[length(b)]] <- quote(mget(grep("^p[A-F]$|^p[12]$", ls(), value = TRUE))); body(f) <- b; f() }
whole <- function(id) FIGS[[id]]$fun()
COMPOSITE <- c("F01", "F02", "F03", "F06", "F07", "F10", "F11", "F12", "F14", "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F22")

## ---------------- text scan of every panel ----------------
BAD <- paste0("\\b(Q[1-4]|A[1-3]|B[1-4]|C[12]|T[1-3]|S[1-8]|F[0-9]{2})\\b|P2_|FAILED|Amendment|nohyper|nocaps|nodrivers|norescue|",
              "uniform_restricted|sensB|option A|\\bLOSS\\b|\\bCONTROL\\b|\\bTRUNC\\b|\\bchat\\b|\\bitem\\b|Not_clear|_worse|_better|calib_|\\bM[123]_|\\brun[0-9]")
txt_of <- function(p) {
  out <- character(0)
  lb <- tryCatch(ggplot2::get_labs(p), error = function(e) tryCatch(p$labels, error = function(e2) list()))
  out <- c(out, unlist(lapply(lb, function(z) if (is.character(z)) z else NULL)))
  b <- tryCatch(suppressWarnings(ggplot_build(p)), error = function(e) NULL)
  if (!is.null(b)) {
    for (d in b$data) if ("label" %in% names(d)) out <- c(out, as.character(d$label))
    for (pp in b$layout$panel_params) for (ax in c("x", "y")) { s <- pp[[ax]]
      if (!is.null(s)) out <- c(out, tryCatch(as.character(unlist(s$get_labels())), error = function(e) character(0))) }
    ly <- as.data.frame(b$layout$layout)
    for (v in setdiff(names(ly), c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y", "AXIS_X", "AXIS_Y"))) out <- c(out, as.character(ly[[v]]))
    for (s in b$plot$scales$scales) out <- c(out, tryCatch(as.character(s$get_labels()), error = function(e) character(0)),
                                             if (is.character(s$name)) s$name else character(0))
  }
  out <- unique(out[!is.na(out)]); out[nzchar(out)]
}
scan <- list(); nhit <- 0L
for (id in names(FIGS)) {
  lv <- tryCatch(suppressWarnings(if (id %in% COMPOSITE) panels(id) else list(main = whole(id))), error = function(e) { lg("  scan could not build ", id, ": ", conditionMessage(e)); list() })
  for (pn in names(lv)) { tx <- txt_of(lv[[pn]]); if (!length(tx)) next
    h <- grepl(BAD, tx, perl = TRUE); nhit <- nhit + sum(h)
    scan[[length(scan) + 1]] <- data.frame(panel_source = paste(id, pn), text = gsub("[\r\n\t]+", " | ", tx), working_term = ifelse(h, "HIT", ""), stringsAsFactors = FALSE) }
}
scan <- do.call(rbind, scan)
write.table(scan, file.path(OUT, "P2_25_text_scan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
lg("Text scan: ", nrow(scan), " text elements, ", nhit, " with a working-term pattern", if (nhit) " (see P2_25_text_scan.tsv, column working_term)" else "")
if (nhit) for (i in which(scan$working_term == "HIT")) lg("  HIT ", scan$panel_source[i], ": ", scan$text[i])

## ---------------- figure layouts ----------------
## Panel letter and panel title on one line: the letter sits in the top-left corner of each panel and the
## title and subtitle start just after it; used for every lettered figure, also those built inside SqnPAN_24 (tagit).
TITLE_LINE <- theme(plot.title.position = "plot", plot.title = element_text(margin = margin(0, 0, 2, 20)),
                    plot.subtitle = element_text(margin = margin(0, 0, 3, 20)))
TAGTH <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT, hjust = 0, vjust = 1),
               plot.tag.position = "topleft", plot.tag.location = "plot") + TITLE_LINE
tagp <- function(p) p + plot_annotation(tag_levels = "A") & TAGTH
tagit <- tagp
collectp <- function(p) p + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical")
nolegend <- function(p) p + theme(legend.position = "none")
OF <- list()
addo <- function(file, title, w, h, src, fun) OF[[file]] <<- list(file = file, title = title, w = w, h = h, src = src, fun = fun)

addo("Figure_1", "Study design, sequon evidence tiers and the landscape of sequon-changing mutations", 180, 230, "F01 A-C; F02 A-B", function() {
  a <- panels("F01"); b <- panels("F02")
  tagp((free(a$pA) | a$pB) / free(a$pC) / (b$pA | b$pB) + plot_layout(heights = c(1, 2.3, 1.15))) })
addo("Figure_2", "Selection at sequon positions estimated by three frameworks", 180, 230, "F03 A; F21 A-B; F22 selection block", function() {
  a <- panels("F03"); d <- panels("F21")
  assign("F22_BLOCKS", "Selection (log2 rate ratio)", envir = globalenv()); s <- panels("F22"); rm("F22_BLOCKS", envir = globalenv())
  p1 <- nolegend(s$p1 + labs(title = "Selection estimates on one scale",
                             subtitle = "Each estimate divided by its equivalence margin; grey band: inside the margin; germline and constrained Cys/His for comparison") +
    theme(strip.text = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(), plot.margin = margin(3, 5, 14, 3)) + TITLE_LINE)
  syn <- wrap_elements(full = (p1 | s$p2) + plot_layout(widths = c(2, 1.15)))
  tagp(collectp(a$pA / (d$pA | d$pB) / syn + plot_layout(heights = c(1.1, 1, 0.85)))) })
addo("Figure_3", "Germline comparison and positive control at constrained Cys/His", 180, 195, "F06 A-D", function() {
  x <- panels("F06"); xD <- x$pD + scale_y_discrete(labels = scales::label_wrap(22))
  tagp(collectp(free(x$pA) / ((free(x$pB) | xD) + plot_layout(widths = c(1, 1))) / x$pC + plot_layout(heights = c(1.3, 0.9, 0.6)))) })
addo("Figure_4", "Mutational processes behind sequon losses and the cross-check with Fan et al. 2018", 180, 230, "F07 A-C; F10 A-B", function() {
  a <- panels("F07"); f <- panels("F10")
  aB <- a$pB + guides(x = guide_axis(check.overlap = TRUE), y = guide_axis(check.overlap = TRUE))
  tagp(((free(a$pA) | aB) + plot_layout(widths = c(2.5, 1))) / free(a$pC) / ((f$pA | f$pB) + plot_layout(widths = c(1.6, 1))) + plot_layout(heights = c(0.95, 1, 1.2))) })
addo("Figure_5", "Predicted structural and functional impact of sequon losses", 180, 230, "F11 A; F12 A-C", function() {
  x <- panels("F11"); y <- panels("F12")
  tagp(collectp(free(x$pA) / y$pA / ((y$pC | y$pB) + plot_layout(widths = c(2.4, 1))) + plot_layout(heights = c(0.85, 1, 1.1)))) })
addo("Figure_6", "Protein amount in CPTAC carriers of sequon losses and clinical links", 180, 230, "F20 A-D; F14 A, C", function() {
  c2 <- panels("F20"); k <- panels("F14")
  tagp((free(c2$pA) | c2$pB) / c2$pC / ((free(c2$pD) | k$pC) + plot_layout(widths = c(1, 1.1))) / nolegend(k$pA) + plot_layout(heights = c(1, 1.15, 0.9, 0.62))) })
addo("Figure_7", "Sequon-creating substitutions, shortlist and glycan modelling", 180, 230, "F16 A-B; F17 A-D", function() {
  g <- panels("F16"); s <- panels("F17")
  gB <- g$pB + guides(colour = guide_legend(nrow = 3), shape = guide_legend(nrow = 3))
  sD <- s$pD + coord_cartesian(clip = "off") + theme(plot.margin = margin(3, 16, 3, 3))
  gA <- nolegend(g$pA)
  tagp((gA | gB) / (s$pA | s$pB) / ((free(s$pC) | sD) + plot_layout(widths = c(1, 1.6))) + plot_layout(heights = c(1.1, 0.8, 1.25))) })
addo("Figure_8", "Single-position hotspots at glycosites and ActiveDriver test of glycosite regions", 180, 230, "F18 A-D; F19 A-C", function() {
  h <- panels("F18"); a <- panels("F19")
  hB <- h$pB + labs(title = "Tier 1 glycosites against background") + scale_y_discrete(labels = scales::label_wrap(24))
  tagp(((free(h$pA) | hB) + plot_layout(widths = c(1.15, 1))) / (free(h$pC) | h$pD) / (((a$pA / a$pC) | a$pB) + plot_layout(widths = c(1, 1.35))) + plot_layout(heights = c(1, 0.8, 1.25))) })

addo("Figure_S1", "Tier 1 sequon losses by cancer type and the most recurrent losses", 180, 125, "F02 C-D", function() { b <- panels("F02"); tagp(b$pC | b$pD) })
addo("Figure_S2", "Substitution spectra behind the expected rates", 180, 110, "F09", function() whole("F09"))
addo("Figure_S3", "Selection test across set-ups, calibration of the expected rate and planted effects", 180, 150, "F03 B-D", function() {
  a <- panels("F03"); tagp(collectp(((a$pB | a$pC) + plot_layout(widths = c(1.2, 1))) / a$pD + plot_layout(heights = c(1, 0.75)))) })
addo("Figure_S4", "Tier 1 sequon losses by cancer type, position-count model", 180, 160, "F05", function() whole("F05"))
addo("Figure_S5", "Tier 1 sequon losses by cancer type, context-aware model", 180, 160, "F04", function() whole("F04"))
addo("Figure_S6", "dNdScv absolute dN/dS, planted effects and null control", 180, 90, "F21 C-E", function() { d <- panels("F21"); tagp(d$pC | d$pD | d$pE) })
addo("Figure_S7", "What each test can detect: planted effects in all pooled units", 180, 175, "F24", function() whole("F24"))
addo("Figure_S8", "Labels of every pooled unit in every set-up", 180, 200, "F23", function() whole("F23"))
addo("Figure_S9", "Signature shares of Tier 1 sequon losses by cancer type", 180, 175, "F08", function() whole("F08"))
addo("Figure_S10", "Structural location of loss and matched control residues", 180, 70, "F11 B", function() panels("F11")$pB)
addo("Figure_S11", "Aromatic residue before the sequon", 180, 95, "F13", function() whole("F13"))
addo("Figure_S12", "Protein amount by CPTAC cancer type and glycopeptides at lost sites", 180, 80, "F20 E-F", function() { c2 <- panels("F20"); tagp(c2$pE | c2$pF) })
addo("Figure_S13", "Progression-free interval by cancer type", 180, 100, "F14 B", function() panels("F14")$pB)
addo("Figure_S14", "Clinical analysis: sensitivity analyses, leave-one-type-out and power", 180, 215, "F15 A-C", function() whole("F15"))
addo("Figure_S15", "Calibration of the gain model and FoldX of sequon-creating substitutions", 180, 80, "F16 C-D", function() { g <- panels("F16"); tagp(collectp(g$pC | g$pD)) })
addo("Figure_S16", "Top background hotspots", 180, 90, "F18 E", function() panels("F18")$pE)
addo("Figure_S17", "All headline estimates scaled to their equivalence margins", 180, 125, "F22", function() whole("F22"))

## ---------------- output ----------------
dev_print <- function(open_fun, p) { open_fun(); tryCatch(print(p), finally = grDevices::dev.off()); invisible(TRUE) }
save_fig <- function(p, base, w, h) {
  if (MODE == "draft") { f <- file.path(OUT, paste0(base, ".png"))
    dev_print(function() ragg::agg_png(f, width = w, height = h, units = "mm", res = 300, background = "white"), p); return(f) }
  f1 <- file.path(OUT, "pdf", paste0(base, ".pdf")); f2p <- file.path(OUT, "tiff", paste0(base, ".tiff")); f3 <- file.path(OUT, "png", paste0(base, ".png"))
  dev_print(function() grDevices::cairo_pdf(f1, width = w / 25.4, height = h / 25.4, family = FONT), p)
  dev_print(function() ragg::agg_tiff(f2p, width = w, height = h, units = "mm", res = 600, compression = "lzw", background = "white"), p)
  dev_print(function() ragg::agg_png(f3, width = w, height = h, units = "mm", res = 300, background = "white"), p)
  paste(f1, f2p, f3, sep = " | ")
}
idx <- list(); n_ok <- 0L
for (f in OF) {
  if (length(ONLY) && !(f$file %in% ONLY)) next
  t0 <- Sys.time(); nw <- 0L
  st <- tryCatch(withCallingHandlers({ p <- f$fun(); save_fig(p, f$file, f$w, f$h); "made" },
                                     warning = function(w) { nw <<- nw + 1L; lg("  warning in ", f$file, ": ", conditionMessage(w)); invokeRestart("muffleWarning") }),
                 error = function(e) { while (grDevices::dev.cur() > 1) grDevices::dev.off(); lg("  ERROR in ", f$file, ": ", conditionMessage(e)); paste0("ERROR: ", conditionMessage(e)) })
  if (identical(st, "made")) n_ok <- n_ok + 1L
  lg(f$file, ": ", st, ", ", nw, " warnings, ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
  idx[[length(idx) + 1]] <- data.frame(file = f$file, title = f$title, width_mm = f$w, height_mm = f$h, panels_from_P2_24 = f$src, status = st, warnings = nw, stringsAsFactors = FALSE)
}
if (length(idx)) write.table(do.call(rbind, idx), file.path(OUT, "P2_25_figure_index.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
lg("SqnPAN_25 finished: ", n_ok, " of ", length(idx), " figures made in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min; folder ", OUT)
