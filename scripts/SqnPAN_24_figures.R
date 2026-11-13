## ===================================================================================================================
## SqnPAN_24_figures.R - figures: theme, palette and the code of every panel (draft renditions; the submitted set is made by SqnPAN_25_figures_final.R)
## Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
##
## Run from the repository folder:
##   Rscript scripts/SqnPAN_24_figures.R draft
##       -> results\figures\P2_24_draft\  one PNG (300 dpi) per figure, P2_24_figure_index.tsv, P2_24_log.txt
##   Rscript scripts/SqnPAN_24_figures.R final
##       -> results\figures\P2_24_final\pdf (cairo PDF), \tiff (600 dpi, LZW), \png (300 dpi)   [after the corrections]
##   Only some figures: add their ids, e.g.  ... SqnPAN_24_figures.R draft F03 F14
##   From RStudio: source("scripts/SqnPAN_24_figures.R") makes the draft set.
##
## Reads only: results\tables\P2_*.tsv (outputs of SqnPAN_09 to SqnPAN_23), data\processed\P2_06_sites.tsv,
##   data\processed\P2_22_pairs.tsv.gz. Writes only inside results\figures\P2_24_draft or P2_24_final.
##   No table, label or result is changed; every label shown is the label written by the analysis scripts,
##   except F15 A/B, where the sensitivity rows get the pre-specified SqnPAN_12 rule (Wald p, then TOST p) as a descriptive label.
## Hard-coded numbers (no table holds them): the B3 shortlist funnel steps 1 to 5 in F17 A (from the SqnPAN_19_gain_side.py log).
## Style: Times New Roman, base 7 pt, one theme, one palette (block PALETTE), width 180 mm.
##   Categorical colours checked with the dataviz palette validator (CVD and normal-vision separation on white).
## Random numbers: point jitter in F20 B (seed 33700001) and label placement by ggrepel (seed 33700002);
##   both in seed block 33700000 to 33799999, used by no other script.
## Packages: ggplot2, patchwork, ggrepel, scales, ragg, systemfonts, survival (installed from CRAN once if missing).
## Each figure is built inside its own error trap: one failing figure is logged and the others are still made.
## ===================================================================================================================

t_start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
MODE <- "draft"; ONLY <- character(0)
if (length(args) >= 1) {
  if (tolower(args[1]) %in% c("draft", "final")) MODE <- tolower(args[1]) else stop("The first argument must be draft or final")
}
if (length(args) >= 2) ONLY <- toupper(args[-1])

ROOT <- Sys.getenv("SQNPAN_ROOT", unset = local({ a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f) == 1) normalizePath(file.path(dirname(normalizePath(f)), ".."), winslash = "/") else normalizePath(getwd(), winslash = "/") }))
if (!dir.exists(ROOT)) stop("repository folder not found: ", ROOT)
TAB <- file.path(ROOT, "results", "tables"); PRO <- file.path(ROOT, "data", "processed")
OUT <- file.path(ROOT, "results", "figures", if (MODE == "draft") "P2_24_draft" else "P2_24_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (MODE == "final") for (s in c("pdf", "tiff", "png")) dir.create(file.path(OUT, s), showWarnings = FALSE)
LOGF <- file.path(OUT, "P2_24_log.txt"); cat("", file = LOGF)
lg <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(...)); message(m); cat(m, "\n", sep = "", file = LOGF, append = TRUE) }

need <- c("ggplot2", "patchwork", "ggrepel", "scales", "ragg", "systemfonts", "survival")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) { message("Installing: ", paste(miss, collapse = ", ")); install.packages(miss, repos = "https://cloud.r-project.org") }
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
lg("SqnPAN_24 start, mode ", MODE, if (length(ONLY)) paste0(", only ", paste(ONLY, collapse = " ")) else "", "; R ", getRversion(),
   ", ggplot2 ", as.character(packageVersion("ggplot2")), ", patchwork ", as.character(packageVersion("patchwork")))

FONT <- "Times New Roman"
fp <- tryCatch(systemfonts::match_fonts(FONT)$path, error = function(e) tryCatch(systemfonts::match_font(FONT)$path, error = function(e2) NA_character_))
lg("Font ", FONT, " -> ", fp)
if (is.na(fp) || !grepl("times", tolower(basename(fp)))) lg("WARNING: Times New Roman was not found by systemfonts; check the figures")

## ---------------- PALETTE (one palette for every figure) ----------------
COL <- list(up = "#B2182B", down = "#2166AC", equiv = "#1A1A1A", inconc = "#9E9E9E", nt = "#C8C8C8", desc = "#5A5A5A",
            N = "#6250d6", ST = "#1baf7a", XP = "#eda100",               # sequon positions 1, 3, 2
            glyco = "#eb6834", control = "#8C8C8C", trunc = "#e87ba4",     # glycosite or loss side, control side, truncating
            T1 = "#A63603", T2 = "#E6550D", T3 = "#FDAE6B",                # evidence tiers (one hue, dark to light)
            band = "#EBEBEB", ref = "#4D4D4D", grid = "#EEEEEE")
CAT8 <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")   # categorical order (fixed)
SEQ_LO <- "#F3F1FC"; SEQ_HI <- "#3F2FA8"     # sequential, violet (counts and shares of planted sets)
SEQ2_LO <- "#FFF3E8"; SEQ2_HI <- "#A63603"   # sequential, orange (signature shares)
BOX <- c(Data = "#EDEDED", Primary = "#FBE1D3", Added = "#E4E0F7")
VLEV <- c("Changed, higher", "Changed, lower", "Equivalent", "Inconclusive or not clear", "Not tested", "Descriptive (no label)")
VCOL <- setNames(c(COL$up, COL$down, COL$equiv, COL$inconc, COL$nt, COL$desc), VLEV)
VSHP <- setNames(c(16, 16, 16, 21, 4, 18), VLEV)

## ---------------- text constants (ASCII source, Unicode escapes) ----------------
U_DD <- "\u0394\u0394G"; U_X <- "\u00d7"; U_PM <- "\u00b1"; U_GE <- "\u2265"; U_LE <- "\u2264"; U_DOT <- "\u00b7"; U_AR <- "\u2192"
EQ <- 0.14; CHG <- 0.32
SUB_EQ <- paste0("Grey band: equivalence region ", U_PM, "0.14; dashed: change threshold ", U_PM, "0.32; thick bar 90% CI, thin bar 95% CI")
XL_RR <- "log2 rate ratio, glycosite against other positions"
SUB_EQ2 <- paste0("Grey band: equivalence ", U_PM, "0.14; dashed: change ", U_PM, "0.32;\nthick bar 90% CI, thin bar 95% CI")

## ---------------- helpers ----------------
fmtn <- function(x) ifelse(is.na(x), "", formatC(round(as.numeric(x)), format = "d", big.mark = ","))
f2 <- function(x, d = 2) ifelse(is.na(x), "", formatC(as.numeric(x), format = "f", digits = d))
pct <- function(x, d = 1) paste0(formatC(100 * as.numeric(x), format = "f", digits = d), "%")
.cache <- new.env()
rt <- function(f, dir = TAB) {
  key <- file.path(dir, f)
  if (!exists(key, envir = .cache, inherits = FALSE)) {
    if (!file.exists(key)) stop("missing input ", key)
    con <- if (grepl("\\.gz$", f)) gzfile(key) else key
    assign(key, read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", "nan", "NaN", "")), envir = .cache)
  }
  get(key, envir = .cache, inherits = FALSE)
}
istrue <- function(x) as.character(x) %in% c("TRUE", "True", "true", "1")
vclass <- function(x) {
  x <- as.character(x); out <- rep(VLEV[4], length(x))
  out[grepl("^Equivalent", x)] <- VLEV[3]
  out[grepl("^(Enriched|Higher|Associated_worse|Hotspot)", x)] <- VLEV[1]
  out[grepl("^(Depleted|Lower|Associated_better)", x)] <- VLEV[2]
  out[grepl("^Descriptive", x)] <- VLEV[6]
  out[is.na(x) | x == "" | grepl("^Not tested|too few", x)] <- VLEV[5]
  factor(out, levels = VLEV)
}
desc_lab <- function(n) factor(rep(VLEV[6], n), levels = VLEV)

GRP <- c(N = "N (position 1)", ST = "S/T (position 3)", both = "N + S/T")
TIER <- c(T1 = "T1 confirmed", T2 = "T2 predicted", T3 = "T3 motif only")
CS_TIER <- c(T1 = "Metal-binding", T2 = "Active site", T3 = "Disulfide"); CS_GRP <- c(N = "Cys", ST = "His", both = "Cys + His")
G_TIER <- c(T1 = "Plausible", T2 = "Plausible or unknown", T3 = "Implausible (comparison)"); G_GRP <- c(N = "new N", ST = "new S/T", both = "new N + new S/T")
unit_lab <- function(tt, gg, TT = TIER, GG = GRP) paste0(TT[tt], " ", U_DOT, " ", GG[gg])
unit_levels <- function(TT = TIER, GG = GRP) as.vector(t(outer(TT, GG, function(a, b) paste0(a, " ", U_DOT, " ", b))))
ufac <- function(x, lev = unit_levels()) factor(x, levels = rev(lev[lev %in% x]))
CFG <- c(main = "Main", strict = "Structure-matched strata", nohyper = "Without hypermutated tumours", union = "All primary-tumour MAF files",
         norescue = "Without rescued rows", uniform_restricted = "Position counts, A1 protein set", sensB = "Germline: two or more resources")
cfac <- function(x) factor(CFG[x], levels = rev(CFG[names(CFG) %in% x]))
eng <- function(d) data.frame(config = d$config, tier = d$tier, group = d$group, scope = d$scope, tested = d$tested,
                              est = d$log2rr, lo95 = d$lo95, hi95 = d$hi95, lo90 = d$lo90, hi90 = d$hi90, fdr = d$fdr,
                              final = d$final_label, label = vclass(d$final_label), nglyco = d$n_events_glyco, stringsAsFactors = FALSE)

## ---------------- theme ----------------
TXT <- 2.2
gt <- function(...) geom_text(..., family = FONT)
theme_p2 <- function(base = 7) theme_classic(base_size = base, base_family = FONT) + theme(
  text = element_text(family = FONT, colour = "black"),
  axis.text = element_text(size = base - 0.6, colour = "black"), axis.title = element_text(size = base - 0.2),
  axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
  plot.title = element_text(size = base + 0.8, face = "bold", hjust = 0, margin = margin(0, 0, 2, 0)),
  plot.subtitle = element_text(size = base - 1, colour = "grey25", hjust = 0, lineheight = 0.95, margin = margin(0, 0, 3, 0)),
  plot.tag = element_text(size = 10, face = "bold"),
  strip.background = element_rect(fill = "#F2F2F2", colour = NA), strip.text = element_text(size = base - 0.6, face = "bold", margin = margin(2, 2, 2, 2)),
  legend.title = element_text(size = base - 0.8, face = "bold"), legend.text = element_text(size = base - 0.8),
  legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1.5, "mm"),
  panel.grid = element_blank(), plot.margin = margin(3, 5, 3, 3))
theme_set(theme_p2())
TT_VOID <- theme(plot.title = element_text(size = 7.8, face = "bold", hjust = 0, family = FONT, margin = margin(0, 0, 2, 0)),
                 plot.subtitle = element_text(size = 6, colour = "grey25", hjust = 0, family = FONT, lineheight = 0.95, margin = margin(0, 0, 3, 0)))
tagit <- function(p) p + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
sc_v <- function() list(scale_colour_manual(values = VCOL, limits = VLEV, drop = FALSE, name = NULL, guide = guide_legend(nrow = 1)),
                        scale_shape_manual(values = VSHP, limits = VLEV, drop = FALSE, name = NULL, guide = guide_legend(nrow = 1)))

## Forest plot used by most figures. d needs y (factor), est, label (vclass factor); optional lo90/hi90, lo95/hi95.
forest <- function(d, xlab, band = c(-EQ, EQ), lines = c(-CHG, CHG), ref = 0, bands = NULL, facet = NULL, ncol = NULL, scales = "fixed",
                   xlim = NULL, log_x = FALSE, title = NULL, subtitle = NULL, legend = "bottom", ntext = NULL, pre = NULL) {
  p <- ggplot(d, aes(x = est, y = y))
  bd <- bands
  if (is.null(bd) && !is.null(band)) bd <- data.frame(xmin = band[1], xmax = band[2])
  if (!is.null(bd)) p <- p + geom_rect(data = bd, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf, inherit.aes = FALSE, fill = COL$band, colour = NA)
  if (length(lines)) p <- p + geom_vline(xintercept = lines, linetype = "22", linewidth = 0.25, colour = "grey60")
  p <- p + geom_vline(xintercept = ref, linewidth = 0.3, colour = COL$ref)
  if (!is.null(pre)) p <- p + pre
  if (all(c("lo95", "hi95") %in% names(d))) p <- p + geom_linerange(aes(xmin = lo95, xmax = hi95, colour = label), linewidth = 0.3, na.rm = TRUE, show.legend = TRUE)
  if (all(c("lo90", "hi90") %in% names(d))) p <- p + geom_linerange(aes(xmin = lo90, xmax = hi90, colour = label), linewidth = 0.95, na.rm = TRUE, show.legend = TRUE)
  p <- p + geom_point(aes(colour = label, shape = label), size = 1.5, fill = "white", stroke = 0.5, na.rm = TRUE, show.legend = TRUE) + sc_v()
  if (!is.null(ntext)) p <- p + gt(aes(x = Inf, label = .data[[ntext]]), hjust = 1.03, size = 1.85, colour = "grey30", na.rm = TRUE)
  if (log_x) p <- p + scale_x_log10()
  if (!is.null(facet)) p <- p + facet_wrap(facet, ncol = ncol, scales = scales)
  p + coord_cartesian(xlim = xlim) + labs(x = xlab, y = NULL, title = title, subtitle = subtitle) +
    theme(legend.position = legend, panel.grid.major.y = element_line(colour = COL$grid, linewidth = 0.2))
}

## Planted-effect capability (share of 20 planted sets called as expected)
SCEN <- setNames(c(paste0("0.67", U_X), paste0("0.8", U_X), paste0("1.25", U_X), paste0("1.5", U_X), "None"),
                 c("run1_changed_right_0.667", "run1_changed_right_0.8", "run1_changed_right_1.25", "run1_changed_right_1.5", "run1_equivalent_1.0"))
ALL_UNITS <- unique(c(unit_levels(), unit_levels(CS_TIER, CS_GRP), unit_levels(G_TIER, G_GRP)))
cap_long <- function(d, analysis, TT = TIER, GG = GRP, tiers = c("T1", "T2", "T3")) {
  d <- d[d$config == "main" & d$scope == "pooled" & d$tier %in% tiers, ]
  do.call(rbind, lapply(names(SCEN), function(k) data.frame(analysis = analysis, y = unit_lab(d$tier, d$group, TT, GG), scen = SCEN[[k]],
                                                          n = as.numeric(d[[k]]), stringsAsFactors = FALSE)))
}
captile <- function(d, title, subtitle = NULL, ncol = NULL) {
  d$analysis <- factor(d$analysis, levels = unique(d$analysis)); d$scen <- factor(d$scen, levels = SCEN); d$y <- factor(d$y, levels = rev(ALL_UNITS[ALL_UNITS %in% d$y]))
  ggplot(d, aes(scen, y, fill = n / 20)) + geom_tile(colour = "white", linewidth = 0.5) +
    gt(aes(label = fmtn(n), colour = n >= 12), size = 1.9) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = SEQ_LO, high = SEQ_HI, limits = c(0, 1), breaks = c(0, 0.5, 1), labels = scales::label_percent(), na.value = "#F7F7F7",
                        name = "Called as expected", guide = guide_colourbar(title.position = "top", barwidth = unit(28, "mm"), barheight = unit(2.5, "mm"))) +
    facet_wrap(~analysis, ncol = ncol, scales = "free_y") +
    labs(x = "Planted rate ratio (0.67 and 0.8: expected Depleted; 1.25 and 1.5: expected Enriched; None: expected Equivalent)", y = NULL, title = title, subtitle = subtitle) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 5.6))
}

## ---------------- figure registry ----------------
FIGS <- list()
add_fig <- function(id, name, title, w, h, src, fun) FIGS[[id]] <<- list(id = id, name = name, title = title, w = w, h = h, src = src, fun = fun)

## =================================================================================================================
## F01 Study design
## =================================================================================================================
add_fig("F01", "study_design", "Study design: data, sequon tiers and analyses", 180, 170,
        "P2_06_sites.tsv; P2_09_Q1_counts.tsv; P2_09_Q1_loss_position_in_sequon.tsv; P2_22_counts.tsv", function() {
  s <- rt("P2_06_sites.tsv", PRO)
  tc <- as.data.frame(table(tier = s$tier), stringsAsFactors = FALSE); tc <- tc[tc$tier %in% names(TIER), ]
  tc$lab <- factor(TIER[tc$tier], levels = rev(TIER))
  ps <- rt("P2_09_Q1_loss_position_in_sequon.tsv"); ps <- ps[ps$set == "main" & ps$scope == "pooled" & ps$tier == "T1", ]
  nl <- setNames(ps$n_sequon_losses, ps$sequon_position)
  q <- rt("P2_09_Q1_counts.tsv"); qa <- q[q$set == "main" & q$scope == "pooled", ]
  ncase <- qa$scope_cases_total[1]; nmis <- qa$n_mutations[qa$category == "all_missense"]
  cc <- rt("P2_22_counts.tsv"); tm <- cc[cc$item == "tumours_with_used_MAF", ]
  bx <- data.frame(x = 1:3, lab = c("N", "X", "S/T"), sub = c("position 1", "position 2\n(X is not P)", "position 3"),
                   loss = c(paste0("Loss: N to another residue\n", fmtn(nl["1 (N)"]), " T1 events"),
                            paste0("Loss: X to P\n", fmtn(nl["2 (X to P)"]), " T1 events"),
                            paste0("Loss: S/T to non-S/T\n", fmtn(nl["3 (S/T)"]), " T1 events")),
                   col = c(COL$N, COL$XP, COL$ST), stringsAsFactors = FALSE)
  pA <- ggplot(bx) + geom_tile(aes(x, 2.25, fill = col, colour = col), width = 0.86, height = 0.62, alpha = 0.22, linewidth = 0.8) +
    scale_fill_identity() + scale_colour_identity() +
    gt(aes(x, 2.25, label = lab), size = 4.2, fontface = "bold") +
    gt(aes(x, 1.8, label = sub), size = 2.0, lineheight = 0.9, vjust = 1) +
    gt(aes(x, 1.05, label = loss), size = 2.0, lineheight = 0.9, vjust = 1, colour = "grey15") +
    annotate("text", x = 2, y = 2.9, label = "N-glycosylation sequon N-X-S/T", size = 2.4, fontface = "bold", family = FONT) +
    annotate("text", x = 2, y = 0.12, label = "Gain: a substitution that creates a new N or a new S/T\ncompletes a new N-X-S/T (B3)",
             size = 2.0, colour = "grey15", lineheight = 0.9, family = FONT) +
    scale_x_continuous(limits = c(0.45, 3.55)) + scale_y_continuous(limits = c(-0.2, 3.05)) +
    theme_void(base_family = FONT) + TT_VOID + labs(title = "Sequon and mutation classes")
  pB <- ggplot(tc, aes(Freq, lab, fill = tier)) + geom_col(width = 0.62) + gt(aes(label = fmtn(Freq)), hjust = -0.15, size = TXT) +
    scale_fill_manual(values = c(T1 = COL$T1, T2 = COL$T2, T3 = COL$T3), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2)), labels = scales::label_comma()) +
    labs(x = "Sequons (N-X-S/T, X not P), UniProt 2026_03", y = NULL, title = "Sequon evidence tiers",
         subtitle = "T1: experimental UniProt evidence or N-GlycositeAtlas site\nT2: other UniProt N-glycosylation annotation; T3: motif only")
  kinds <- c(rep("Data", 6), rep("Primary", 5), rep("Added", 9))
  bxs <- data.frame(col = c(rep(1, 6), rep(2, 5), rep(3, 5), rep(4, 4)), row = c(1:6, 1:5, 1:5, 1:4), kind = kinds, stringsAsFactors = FALSE,
    txt = c(sprintf("TCGA somatic mutations\nGDC masked MAFs (GRCh38)\n%s tumours, 33 cancer types", fmtn(ncase)),
            sprintf("%s missense mutations\nmapped to UniProt positions", fmtn(nmis)),
            "Sequons and evidence tiers\nUniProt 2026_03,\nN-GlycositeAtlas",
            "AlphaFold models\n(pLDDT, exposure); FoldX,\nAlphaMissense, ESM1b",
            "Germline missense variants\n(gnomAD, TOPMed, ExAC,\n1000 Genomes, ESP)",
            sprintf("TCGA-CDR survival (PFI)\nCPTAC proteomes: %s tumours\nwith a MAF, %d cancer types", fmtn(sum(as.numeric(tm$value))), nrow(tm)),
            "Q1 Catalogue of\nsequon losses and gains", "Q2 Selection test:\nglycosite against\nother N, S/T positions",
            paste0("Q3 Structure: FoldX ", U_DD, "\nagainst matched\ncontrol substitutions"), "Q4 Clinical links:\nstratified Cox model\nof PFI",
            "Cross-check with\nFan et al. 2018",
            "A1 Context-aware\nselection test\n(trinucleotide rates)", "A2 Germline control;\nconstrained Cys/His\ncheck",
            "A3 Mutational\nsignatures", "B1 AlphaMissense\nand ESM1b", "B2 CPTAC protein\namount in carriers",
            "B3 Gain side, shortlist,\nNetNGlyc, Re-Glyco", "B4 Hotspots and\nActiveDriver", "C1 Aromatic\nsequon split", "dNdScv\n(Martincorena 2017)"))
  bxs$x <- c(0.62, 1.78, 2.94, 4.06)[bxs$col]; bxs$y <- 6.4 - bxs$row
  hd <- data.frame(x = c(0.62, 1.78, 3.5), y = 6.18, lab = c("Data", "Primary analysis (pre-specified)", "Added after the primary analysis"))
  pC <- ggplot(bxs) + geom_tile(aes(x, y, fill = kind), width = 1.08, height = 0.84, colour = "grey55", linewidth = 0.25) +
    gt(aes(x, y, label = txt), size = 1.95, lineheight = 0.88) + gt(data = hd, aes(x, y, label = lab), size = 2.3, fontface = "bold") +
    scale_fill_manual(values = BOX, guide = "none") + scale_x_continuous(limits = c(0.05, 4.63)) + scale_y_continuous(limits = c(-0.1, 6.35)) +
    theme_void(base_family = FONT) + TT_VOID + labs(title = "Data and analyses", subtitle = "Added analyses are post hoc (Amendment 4 and dNdScv), each with its own controls; no primary label changed")
  tagit((pA | pB) / pC + plot_layout(heights = c(1, 2.3)))
})

## =================================================================================================================
## F02 Landscape of sequon-changing mutations (Q1)
## =================================================================================================================
add_fig("F02", "landscape", "Sequon-changing somatic mutations in TCGA (Q1)", 180, 175,
        "P2_09_Q1_counts.tsv; P2_09_Q1_loss_position_in_sequon.tsv; P2_09_Q1_recurrent_variants.tsv", function() {
  q <- rt("P2_09_Q1_counts.tsv"); qa <- q[q$set == "main" & q$scope == "pooled", ]
  cats <- c(loss_T1 = "T1 confirmed", loss_T2 = "T2 predicted", loss_T3 = "T3 motif only", gain_plausible = "Plausible", gain_unknown = "Unknown",
            gain_implausible = "Implausible")
  a <- qa[qa$category %in% names(cats), ]
  a$side <- factor(ifelse(grepl("^loss", a$category), "Losses", "Gains"),
                   levels = c("Losses", "Gains"))
  a$lab <- factor(cats[a$category], levels = rev(cats))
  fills <- c(loss_T1 = COL$T1, loss_T2 = COL$T2, loss_T3 = COL$T3, gain_plausible = SEQ_HI, gain_unknown = "#8A7FE0", gain_implausible = "#CFCAF3")
  pA <- ggplot(a, aes(n_mutations, lab, fill = category)) + geom_col(width = 0.65) +
    gt(aes(label = paste0(fmtn(n_mutations), " (", fmtn(n_cases), " tumours)")), hjust = -0.06, size = 1.95) +
    scale_fill_manual(values = fills, guide = "none") + facet_grid(side ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.5)), labels = scales::label_comma()) +
    labs(x = "Missense mutations", y = NULL, title = "Losses and gains",
         subtitle = sprintf("%s missense mutations in %s tumours (33 types)", fmtn(qa$n_mutations[qa$category == "all_missense"]), fmtn(qa$scope_cases_total[1])))
  ps <- rt("P2_09_Q1_loss_position_in_sequon.tsv"); ps <- ps[ps$scope == "pooled", ]
  pos_map <- c("1 (N)" = "N (position 1)", "2 (X to P)" = "X to P (position 2)", "3 (S/T)" = "S/T (position 3)")
  ps$pos <- factor(pos_map[ps$sequon_position], levels = pos_map)
  ps$setl <- factor(ifelse(ps$set == "main", "All tumours", "Without hypermutated"), levels = c("All tumours", "Without hypermutated"))
  ps$tl <- factor(TIER[ps$tier], levels = rev(TIER))
  ps$share <- ave(ps$n_sequon_losses, ps$set, ps$tier, FUN = function(v) v / sum(v))
  pB <- ggplot(ps, aes(share, tl, fill = pos)) + geom_col(width = 0.65, colour = "white", linewidth = 0.3, position = position_stack(reverse = TRUE)) +
    gt(aes(label = ifelse(share > 0.07, fmtn(n_sequon_losses), "")), position = position_stack(vjust = 0.5, reverse = TRUE), size = 1.9) +
    scale_fill_manual(values = setNames(c(COL$N, COL$XP, COL$ST), pos_map), name = "Lost position") + facet_wrap(~setl, ncol = 1) +
    scale_x_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
    labs(x = "Share of sequon losses", y = NULL, title = "Where in the sequon", subtitle = "Counts inside the bars") + theme(legend.position = "bottom")
  ct <- q[q$set == "main" & q$category == "loss_T1" & q$scope != "pooled", ]
  ct$type <- sub("^TCGA-", "", ct$scope); ct$share <- 100 * ct$n_cases / ct$scope_cases_total
  ct$type <- factor(ct$type, levels = ct$type[order(ct$share)])
  pC <- ggplot(ct, aes(share, type)) + geom_col(fill = COL$glyco, width = 0.7) +
    gt(aes(label = paste0(fmtn(n_cases), "/", fmtn(scope_cases_total))), hjust = -0.1, size = 1.75) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(x = "Tumours with at least one T1 sequon loss (%)", y = NULL, title = "T1 losses by cancer type", subtitle = "Label: tumours with a loss / tumours in the type")
  rv <- rt("P2_09_Q1_recurrent_variants.tsv"); rv <- rv[rv$category == "loss_T1", ]; rv <- rv[order(rv$rank_in_category), ]; rv <- rv[seq_len(min(15, nrow(rv))), ]
  rv$v <- paste0(rv$gene, " ", rv$ref, rv$position, rv$alt); rv$v <- factor(rv$v, levels = rev(unique(rv$v)))
  rl <- rbind(data.frame(v = rv$v, n = rv$n_cases, w = "All tumours"), data.frame(v = rv$v, n = rv$n_cases_without_hypermutated, w = "Without hypermutated"))
  pD <- ggplot(rv, aes(y = v)) + geom_segment(aes(x = n_cases_without_hypermutated, xend = n_cases, yend = v), colour = "grey70", linewidth = 0.45) +
    geom_point(data = rl, aes(x = n, y = v, colour = w, shape = w), size = 1.6, fill = "white", stroke = 0.6) +
    gt(aes(x = n_cases, label = paste0(n_cancer_types, " types")), hjust = -0.35, size = 1.75, colour = "grey30") +
    scale_colour_manual(values = c("All tumours" = COL$glyco, "Without hypermutated" = COL$ref), name = NULL) +
    scale_shape_manual(values = c("All tumours" = 16, "Without hypermutated" = 21), name = NULL) +
    scale_x_continuous(breaks = scales::breaks_pretty(), expand = expansion(mult = c(0.05, 0.35))) +
    labs(x = "Tumours carrying the substitution", y = NULL, title = "Most recurrent T1 losses") + theme(legend.position = "bottom")
  tagit((pA | pB) / (pC | pD) + plot_layout(heights = c(1, 2)))
})

## =================================================================================================================
## F03 Selection test (Q2 primary and A1 context-aware)
## =================================================================================================================
add_fig("F03", "selection_test", "Selection test at sequon positions: Q2 and A1", 180, 210,
        "P2_09_Q2_results.tsv; P2_14_Q2ctx_results.tsv; P2_14_calibration.tsv; P2_08_unit_capability.tsv; P2_14_unit_capability.tsv", function() {
  q2 <- eng(rt("P2_09_Q2_results.tsv")); a1 <- eng(rt("P2_14_Q2ctx_results.tsv"))
  M1 <- "Position-count model (Q2, primary)"; M2 <- "Context-aware model (A1)"
  m <- rbind(cbind(q2[q2$config == "main" & q2$scope == "pooled", ], model = M1), cbind(a1[a1$config == "main" & a1$scope == "pooled", ], model = M2))
  m$model <- factor(m$model, levels = c(M1, M2)); m$y <- ufac(unit_lab(m$tier, m$group)); m$nt <- fmtn(m$nglyco)
  pA <- forest(m, XL_RR, facet = "model", ncol = 2, xlim = c(-0.42, 0.42), ntext = "nt",
               title = "Tumour somatic missense, 33 cancer types pooled", subtitle = paste0(SUB_EQ, "; right: glycosite events"))
  cf <- rbind(cbind(q2[q2$scope == "pooled" & q2$tier == "T1" & q2$group == "both", ], model = M1),
              cbind(a1[a1$scope == "pooled" & a1$tier == "T1" & a1$group == "both", ], model = M2))
  cf$model <- factor(ifelse(cf$model == M1, "Q2 position-count", "A1 context-aware"), levels = c("Q2 position-count", "A1 context-aware")); cf$y <- cfac(cf$config)
  pB <- forest(cf, XL_RR, facet = "model", ncol = 2, xlim = c(-0.25, 0.25), lines = NULL, title = "T1 N + S/T losses across set-ups")
  cb <- rt("P2_14_calibration.tsv"); cb <- cb[cb$scope == "pooled" & cb$tier == "T1" & cb$run == min(cb$run), ]
  cmod <- c(calib_ctx = "Context model", calib_uniform = "Uniform rate")
  cb$est <- cb$log2rr; cb$label <- vclass(cb$raw_label)
  cb$y <- factor(paste0(cmod[cb$config], " ", U_DOT, " ", GRP[cb$group]), levels = rev(as.vector(t(outer(cmod, GRP, function(a, b) paste0(a, " ", U_DOT, " ", b))))))
  pC <- forest(cb, "log2 ratio, high- vs low-expected background", lines = NULL, xlim = c(-0.15, 1.35),
               title = "Calibration of the expected rate (A1)",
               subtitle = "Background split at the median expected rate in each\nprotein and residue; a calibrated model gives about 0")
  cp <- rbind(cap_long(rt("P2_08_unit_capability.tsv"), M1, tiers = "T1"), cap_long(rt("P2_14_unit_capability.tsv"), M2, tiers = "T1"))
  pD <- captile(cp, "What the test can detect (planted effects, T1 units)", "Number of 20 planted data sets labelled as expected; run 1")
  tagit((pA / ((pB | pC) + plot_layout(widths = c(1.2, 1))) / pD + plot_layout(heights = c(1.25, 0.85, 0.55), guides = "collect")) & theme(legend.position = "bottom", legend.box = "vertical"))
})

## =================================================================================================================
## F04 / F05 Per cancer type
## =================================================================================================================
pertype <- function(file, model_name) {
  e <- eng(rt(file)); e <- e[e$config == "main" & e$scope != "pooled" & e$tier == "T1", ]
  e$type <- sub("^TCGA-", "", e$scope); e$y <- factor(e$type, levels = rev(sort(unique(e$type))))
  e$grp <- factor(GRP[e$group], levels = GRP); e$nt <- fmtn(e$nglyco)
  forest(e, XL_RR, facet = "grp", ncol = 3, xlim = c(-1.7, 1.7), ntext = "nt", title = paste0("T1 sequon losses by cancer type: ", model_name),
         subtitle = paste0(SUB_EQ, "\nBH within each column of 33 types; right: glycosite events; types with too few data are Not tested"))
}
add_fig("F04", "per_type_context", "T1 losses by cancer type, context-aware model (A1)", 180, 160, "P2_14_Q2ctx_results.tsv",
        function() pertype("P2_14_Q2ctx_results.tsv", "context-aware model (A1)"))
add_fig("F05", "per_type_position", "T1 losses by cancer type, position-count model (Q2)", 180, 160, "P2_09_Q2_results.tsv",
        function() pertype("P2_09_Q2_results.tsv", "position-count model (Q2)"))

## =================================================================================================================
## F06 Germline control (A2) and constrained-residue check
## =================================================================================================================
add_fig("F06", "germline_control", "Germline control (A2) and constrained Cys/His positive control", 180, 195,
        "P2_14_Q2ctx_results.tsv; P2_15_germline_results.tsv; P2_17_constrained_results.tsv", function() {
  a1 <- eng(rt("P2_14_Q2ctx_results.tsv")); gm <- eng(rt("P2_15_germline_results.tsv")); cs <- eng(rt("P2_17_constrained_results.tsv"))
  S1 <- "Tumour somatic, glycosites (A1)"; S2 <- "Germline, glycosites (A2)"
  x <- rbind(cbind(a1[a1$config == "main" & a1$scope == "pooled", ], src = S1), cbind(gm[gm$config == "main" & gm$scope == "pooled", ], src = S2))
  x$src <- factor(x$src, levels = c(S1, S2)); x$y <- ufac(unit_lab(x$tier, x$group)); x$nt <- fmtn(x$nglyco)
  pA <- forest(x, XL_RR, facet = "src", ncol = 2, xlim = c(-0.36, 0.36), ntext = "nt", title = "Same test code on tumour and germline variants",
               subtitle = paste0(SUB_EQ, "; right: glycosite events\nGermline variation is not depleted at glycosite N, so the tumour Equivalent reads Equivalent_biology_control_FAILED"))
  c1 <- cs[cs$config == "main" & cs$scope == "pooled" & !is.na(cs$est), ]
  c1$y <- ufac(unit_lab(c1$tier, c1$group, CS_TIER, CS_GRP), unit_levels(CS_TIER, CS_GRP)); c1$nt <- fmtn(c1$nglyco)
  pB <- forest(c1, "log2 rate ratio, constrained residue against other Cys/His", xlim = c(-1, 0.35), ntext = "nt",
               title = "Positive control: germline, constrained Cys/His",
               subtitle = "UniProt metal-binding, active-site and disulfide residues;\nthe test finds the expected depletion (decision unit: metal-binding Cys + His)")
  c2 <- cs[cs$scope == "pooled" & cs$tier == "T1" & !is.na(cs$est), ]; c2$y <- cfac(c2$config); c2$grp <- factor(CS_GRP[c2$group], levels = CS_GRP)
  pD <- forest(c2, "log2 rate ratio", facet = "grp", ncol = 3, xlim = c(-1.05, 0.2), title = "Metal-binding Cys/His across set-ups")
  g2 <- gm[gm$scope == "pooled" & gm$tier == "T1", ]; g2$y <- cfac(g2$config); g2$grp <- factor(GRP[g2$group], levels = GRP)
  pC <- forest(g2, XL_RR, facet = "grp", ncol = 3, xlim = c(-0.2, 0.2), lines = NULL, title = "Germline T1 glycosites across set-ups")
  tagit((pA / (pB | pD) / pC + plot_layout(heights = c(1.3, 0.9, 0.6), guides = "collect")) & theme(legend.position = "bottom", legend.box = "vertical"))
})

## =================================================================================================================
## F07 / F08 Mutational signatures (A3)
## =================================================================================================================
SIG <- c("SBS1", "SBS5", "APOBEC", "Tobacco", "UV", "MMR", "POLE", "other")
SIGL <- c(SBS1 = "SBS1 (clock-like)", SBS5 = "SBS5 (clock-like)", APOBEC = "APOBEC", Tobacco = "Tobacco", UV = "UV", MMR = "MMR deficiency",
          POLE = "POLE", other = "Other", CpG_CtoT = "CpG C>T")
SIGCOL <- setNames(c(CAT8[1:7], "#C8C8C8"), SIGL[SIG])
SETL <- c(T1_loss_N = "T1 loss, N", T1_loss_ST = "T1 loss, S/T", T1_loss_XP = "T1 loss, X to P", bg_N = "Background N", bg_ST = "Background S/T",
          loss_any = "All sequon losses", gain_plausible = "Gains, plausible", gain_unknown = "Gains, unknown", all_missense = "All missense")
add_fig("F07", "signatures", "Mutational signatures behind sequon losses and gains (A3)", 180, 160,
        "P2_16_shares_pooled.tsv; P2_16_differences.tsv", function() {
  sp <- rt("P2_16_shares_pooled.tsv"); sp <- sp[sp$set %in% names(SETL), ]
  a <- sp[sp$category %in% SIG, ]; a$setl <- factor(SETL[a$set], levels = rev(SETL)); a$cat <- factor(SIGL[a$category], levels = SIGL[SIG])
  nn <- unique(sp[, c("set", "n")]); nn$setl <- factor(SETL[nn$set], levels = rev(SETL))
  pA <- ggplot(a, aes(share, setl, fill = cat)) + geom_col(width = 0.72, colour = "white", linewidth = 0.25, position = position_stack(reverse = TRUE)) +
    gt(data = nn, aes(x = 1.01, y = setl, label = paste0("n = ", fmtn(n))), inherit.aes = FALSE, hjust = 0, size = 1.85) +
    scale_fill_manual(values = SIGCOL, name = NULL) + scale_x_continuous(labels = scales::label_percent(), expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    labs(x = "Share of mutations, pooled over cancer types", y = NULL, title = "Signature groups",
         subtitle = "Per-type signature fits; each mutation shared by its channel probabilities") +
    theme(legend.position = "bottom", plot.margin = margin(3, 32, 3, 3)) + guides(fill = guide_legend(nrow = 2))
  cg <- sp[sp$category == "CpG_CtoT", ]; cg$setl <- factor(SETL[cg$set], levels = rev(SETL))
  pB <- ggplot(cg, aes(share, setl)) + geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.5, colour = COL$ref, na.rm = TRUE) +
    geom_point(size = 1.6, colour = COL$ref) + scale_x_continuous(labels = scales::label_percent()) +
    labs(x = "Share at CpG C>T (95% CI)", y = NULL, title = "CpG C>T") + theme(axis.text.y = element_blank())
  df <- rt("P2_16_differences.tsv")
  CON <- c("T1_loss_ST minus T1_loss_N" = "T1 S/T loss minus T1 N loss", "T1_loss_N minus bg_N" = "T1 N loss minus background N",
           "T1_loss_ST minus bg_ST" = "T1 S/T loss minus background S/T", "loss_any minus all_missense" = "All losses minus all missense",
           "gain_plausible minus all_missense" = "Plausible gains minus all missense")
  df <- df[df$contrast %in% names(CON), ]; df$con <- factor(CON[df$contrast], levels = CON)
  df$cat <- factor(SIGL[df$category], levels = rev(SIGL[c(SIG, "CpG_CtoT")]))
  df$ex <- ifelse(df$lo95 > 0 | df$hi95 < 0, "95% CI excludes 0", "95% CI includes 0")
  pC <- ggplot(df, aes(difference, cat)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = COL$ref) +
    geom_linerange(aes(xmin = lo95, xmax = hi95, colour = ex), linewidth = 0.55) + geom_point(aes(colour = ex), size = 1.3) +
    scale_colour_manual(values = c("95% CI excludes 0" = "#1A1A1A", "95% CI includes 0" = "#A6A6A6"), name = NULL) +
    facet_wrap(~con, nrow = 1, labeller = label_wrap_gen(20)) + scale_x_continuous(labels = function(v) formatC(100 * v, format = "f", digits = 0)) +
    labs(x = "Difference in share, percentage points (95% CI)", y = NULL, title = "Differences between sets (descriptive)",
         subtitle = "The position-3 excess of losses follows CpG C>T, UV and APOBEC exposure, the same processes that act on background S/T") +
    theme(legend.position = "bottom", panel.grid.major.y = element_line(colour = COL$grid, linewidth = 0.2))
  tagit(((pA | pB) + plot_layout(widths = c(3.3, 1))) / pC + plot_layout(heights = c(1, 1.05)))
})
add_fig("F08", "signatures_by_type", "Signature shares of T1 losses by cancer type (A3)", 180, 175, "P2_16_shares_by_type.tsv", function() {
  sb <- rt("P2_16_shares_by_type.tsv"); sb <- sb[sb$set %in% c("T1_loss_N", "T1_loss_ST"), ]
  cats <- c(SIG, "CpG_CtoT")
  L <- do.call(rbind, lapply(cats, function(k) data.frame(set = sb$set, type = sub("^TCGA-", "", sb$project), n = sb$n, cat = SIGL[k], share = as.numeric(sb[[k]]))))
  L$setl <- factor(SETL[L$set], levels = SETL[c("T1_loss_N", "T1_loss_ST")])
  L$cat <- factor(L$cat, levels = SIGL[cats]); L$yl <- paste0(L$type, " (", L$n, ")")
  L$yl <- factor(L$yl, levels = rev(sort(unique(L$yl))))
  ggplot(L, aes(cat, yl, fill = share)) + geom_tile(colour = "white", linewidth = 0.4) +
    gt(aes(label = ifelse(share >= 0.1, formatC(100 * share, format = "f", digits = 0), ""), colour = share >= 0.45), size = 1.7) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = SEQ2_LO, high = SEQ2_HI, limits = c(0, 1), labels = scales::label_percent(), name = "Share") +
    facet_wrap(~setl, nrow = 1, scales = "free_y") +
    labs(x = NULL, y = NULL, title = "Signature shares of T1 sequon losses by cancer type",
         subtitle = "Cancer type (number of losses); numbers in tiles are percentages of 10% or more; CpG C>T is counted beside the signature groups") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
})

## =================================================================================================================
## F09 Substitution spectra
## =================================================================================================================
add_fig("F09", "spectra", "Substitution spectra used by the context model", 180, 110, "P2_14_spectra.tsv; P2_15_germline_spectra.tsv", function() {
  SUBS <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  prep <- function(ch, val, what) {
    d <- data.frame(channel = ch, v = val, stringsAsFactors = FALSE)
    d$sub <- substr(d$channel, 3, 5); d$tri <- paste0(substr(d$channel, 1, 1), substr(d$channel, 3, 3), substr(d$channel, 7, 7))
    d$ord <- match(d$sub, SUBS) * 100 + match(substr(d$channel, 1, 1), c("A", "C", "G", "T")) * 10 + match(substr(d$channel, 7, 7), c("A", "C", "G", "T"))
    d <- d[order(d$ord), ]; d$share <- d$v / sum(d$v); d$what <- what; d$sub <- factor(d$sub, levels = SUBS)
    d$x <- factor(paste(d$sub, d$tri), levels = unique(paste(d$sub, d$tri))); d
  }
  sp <- rt("P2_14_spectra.tsv"); sp <- sp[sp$set == "main", ]; ag <- aggregate(count ~ channel, data = sp, FUN = sum)
  gs <- rt("P2_15_germline_spectra.tsv")
  d <- rbind(prep(ag$channel, ag$count, sprintf("Tumour somatic SNVs\n%s, 33 types", fmtn(sum(ag$count)))),
             prep(gs$channel, gs$variants_A, sprintf("Germline missense\n%s variants", fmtn(sum(gs$variants_A)))))
  d$what <- factor(d$what, levels = unique(d$what))
  ggplot(d, aes(x, share, fill = sub)) + geom_col(width = 0.75) + facet_grid(what ~ sub, scales = "free", space = "free_x") +
    scale_fill_manual(values = setNames(CAT8[1:6], SUBS), guide = "none") + scale_y_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(labels = function(v) sub("^\\S+ ", "", v)) +
    labs(x = "Trinucleotide context (pyrimidine-centred)", y = "Share of substitutions", title = "Substitution spectra behind the expected rates",
         subtitle = "Tumour: chosen-file SNVs per channel (counts pooled over types; the model uses per-type rates); germline: option A variant set") +
    theme(axis.text.x = element_text(angle = 90, size = 3.6, vjust = 0.5, hjust = 1), panel.spacing.x = unit(0.6, "mm"),
          strip.text.y = element_text(size = 5.6))
})

## =================================================================================================================
## F10 Fan et al. 2018 cross-check
## =================================================================================================================
add_fig("F10", "fan_crosscheck", "Cross-check with Fan et al. 2018", 180, 85, "P2_13_fan_counts.tsv; P2_13_fan_positions.tsv", function() {
  fc <- rt("P2_13_fan_counts.tsv"); fc <- fc[!is.na(fc$fan_2018) & !is.na(fc$ours), ]
  fc$q <- factor(trimws(fc$quantity), levels = rev(unique(trimws(fc$quantity))))
  fl <- rbind(data.frame(q = fc$q, v = fc$fan_2018, w = "Fan et al. 2018"), data.frame(q = fc$q, v = fc$ours, w = "This study (TCGA only)"))
  pA <- ggplot(fc, aes(y = q)) + geom_segment(aes(x = fan_2018, xend = ours, yend = q), colour = "grey70", linewidth = 0.45) +
    geom_point(data = fl, aes(x = v, y = q, colour = w, shape = w), size = 1.6) +
    scale_colour_manual(values = c("Fan et al. 2018" = COL$ref, "This study (TCGA only)" = COL$glyco), name = NULL) +
    scale_shape_manual(values = c("Fan et al. 2018" = 17, "This study (TCGA only)" = 16), name = NULL) +
    scale_x_log10(labels = scales::label_comma()) + scale_y_discrete(labels = scales::label_wrap(40)) + labs(x = "Count (log scale)", y = NULL, title = "Counts beside Fan et al. 2018",
         subtitle = "Different data sources and sequon definitions; not a test") + theme(legend.position = "bottom")
  fp <- rt("P2_13_fan_positions.tsv")
  L <- rbind(data.frame(tier = fp$tier, v = fp$main_ratio_3_to_1, w = "Observed, all tumours"),
             data.frame(tier = fp$tier, v = fp$nohyper_ratio_3_to_1, w = "Observed, without hypermutated"),
             data.frame(tier = fp$tier, v = fp$expected_ratio_ST_to_N, w = "Expected from background rates"))
  L$tl <- factor(TIER[L$tier], levels = rev(TIER)); L$w <- factor(L$w, levels = unique(L$w))
  L$yn <- as.numeric(L$tl) + (2 - as.numeric(L$w)) * 0.22
  pB <- ggplot(L, aes(v, yn, colour = w, shape = w)) + geom_point(size = 1.9, stroke = 0.6) +
    scale_y_continuous(breaks = seq_along(levels(L$tl)), labels = levels(L$tl), expand = expansion(add = 0.45)) +
    scale_colour_manual(values = setNames(c(COL$glyco, COL$N, COL$ref), levels(L$w)), name = NULL) +
    scale_shape_manual(values = setNames(c(16, 17, 23), levels(L$w)), name = NULL) +
    scale_x_continuous(limits = c(1, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Ratio of S/T to N losses", y = NULL, title = "The position-3 excess",
         subtitle = "More S/T than N losses (Fan et al.);\nthe background rates expect it") +
    theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
  tagit((pA | pB) + plot_layout(widths = c(1.25, 1)))
})

## =================================================================================================================
## F11 Structure (Q3, FoldX)
## =================================================================================================================
UNIT6 <- c("1" = "N + S/T losses (primary)", "2" = "N losses", "3" = "S/T losses", "4" = "Without fragment models",
           "5" = "Stricter matching (same mutant residue)", "6" = "Second control draw")
add_fig("F11", "structure_foldx", "Structural impact of T1 sequon losses (Q3, FoldX)", 180, 150, "P2_11_foldx_ddg.tsv; P2_11_comparison_results.tsv", function() {
  fx <- rt("P2_11_foldx_ddg.tsv")
  lo <- fx[fx$set == "loss" & fx$role %in% c("N", "ST"), ]; cl <- fx[fx$set == "ctrl_d1" & fx$role %in% c("N", "ST"), ]
  RL <- c(N = "N losses (position 1)", ST = "S/T losses (position 3)")
  d <- rbind(data.frame(grp = "T1 sequon loss", role = lo$role, ddg = lo$ddg_mean, plddt = lo$plddt, rsa = lo$rsa),
             data.frame(grp = "Matched control substitution", role = cl$role, ddg = cl$ddg_mean, plddt = cl$plddt, rsa = cl$rsa))
  d$grp <- factor(d$grp, levels = c("T1 sequon loss", "Matched control substitution")); d$role <- factor(RL[d$role], levels = RL)
  GC <- c("T1 sequon loss" = COL$glyco, "Matched control substitution" = COL$control)
  nn <- aggregate(ddg ~ grp + role, data = d, FUN = length)
  pA <- ggplot(d, aes(ddg, colour = grp)) + geom_vline(xintercept = c(1, 3), linetype = "22", linewidth = 0.25, colour = "grey60") +
    stat_ecdf(linewidth = 0.55, pad = FALSE, na.rm = TRUE) + facet_wrap(~role, nrow = 1) + scale_colour_manual(values = GC, name = NULL) +
    coord_cartesian(xlim = c(-2, 6)) + scale_y_continuous(labels = scales::label_percent()) +
    labs(x = paste0("FoldX ", U_DD, " (kcal/mol; > 0 destabilising)"), y = "Cumulative share", title = paste0("FoldX ", U_DD, " of losses and matched controls"),
         subtitle = sprintf("Losses %s (N) and %s (S/T); controls from draw 1; dashed: 1 and 3 kcal/mol",
                            fmtn(sum(d$grp == "T1 sequon loss" & d$role == RL["N"])), fmtn(sum(d$grp == "T1 sequon loss" & d$role == RL["ST"])))) +
    theme(legend.position = "bottom")
  L <- rbind(data.frame(grp = d$grp, m = "pLDDT at the residue", v = d$plddt), data.frame(grp = d$grp, m = "Relative solvent exposure", v = d$rsa))
  pB <- ggplot(L, aes(grp, v, fill = grp)) + geom_violin(colour = NA, alpha = 0.55, scale = "width", na.rm = TRUE) +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3, fill = "white", na.rm = TRUE) + facet_wrap(~m, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = GC, guide = "none") + scale_x_discrete(labels = c("Loss", "Control")) +
    labs(x = NULL, y = NULL, title = "Where the residues lie", subtitle = "AlphaFold models; only residues with pLDDT >= 70 were modelled")
  cr <- rt("P2_11_comparison_results.tsv")
  cr$est <- cr$median_d; cr$lo95 <- cr$ci95_lo; cr$hi95 <- cr$ci95_hi; cr$lo90 <- cr$ci90_lo; cr$hi90 <- cr$ci90_hi; cr$label <- vclass(cr$label)
  cr$y <- factor(UNIT6[as.character(cr$unit)], levels = rev(UNIT6)); cr$nt <- paste0(fmtn(cr$n_pairs), " pairs")
  pC <- forest(cr, paste0("Median difference in ", U_DD, ", loss minus matched controls (kcal/mol)"), band = c(-0.5, 0.5), lines = NULL, xlim = c(-0.62, 0.62),
               ntext = "nt", title = "Losses against matched controls", subtitle = paste0("Grey band: pre-specified margin ", U_PM, "0.5 kcal/mol;\nthick 90% CI, thin 95% CI"))
  tagit((pA / (pB | pC) + plot_layout(heights = c(1, 1), guides = "keep")))
})

## =================================================================================================================
## F12 Variant-effect scores (B1) and F13 aromatic split (C1)
## =================================================================================================================
MEAS <- c(FoldX = paste0("FoldX ", U_DD, " (kcal/mol)"), AM = "AlphaMissense pathogenicity", ESM = "ESM1b damage (-LLR)")
add_fig("F12", "variant_effect", "AlphaMissense and ESM1b beside FoldX (B1)", 180, 175, "P2_18_A_results.tsv; P2_18_A_correlation.tsv; P2_18_B_results.tsv", function() {
  ar <- rt("P2_18_A_results.tsv"); ar <- ar[ar$run == 1, ]
  ar$est <- ar$median_d; ar$lo95 <- ar$ci95_lo; ar$hi95 <- ar$ci95_hi; ar$lo90 <- ar$ci90_lo; ar$hi90 <- ar$ci90_hi; ar$label <- vclass(ar$label)
  ar$meas <- factor(MEAS[ar$measure], levels = MEAS); ar$y <- factor(UNIT6[as.character(ar$unit)], levels = rev(UNIT6))
  bd <- unique(data.frame(meas = ar$meas, xmin = -ar$margin, xmax = ar$margin))
  pA <- forest(ar, "Median difference, loss minus matched controls (> 0: loss predicted more damaging)", bands = bd, lines = NULL,
               facet = "meas", ncol = 3, scales = "free_x", title = "TCGA T1 losses against matched controls",
               subtitle = "Grey band: pre-specified margin of each score (FoldX 0.5, AlphaMissense 0.082, ESM1b 0.71); thick 90% CI, thin 95% CI")
  cr <- rt("P2_18_A_correlation.tsv")
  cl <- rbind(data.frame(pair = paste0(cr$measure_x, " and ", cr$measure_y), v = cr$spearman_d, w = "Paired differences"),
              data.frame(pair = paste0(cr$measure_x, " and ", cr$measure_y), v = cr$spearman_loss_value, w = "Loss values"))
  cl$pair <- gsub("AM", "AlphaMissense", gsub("ESM", "ESM1b", cl$pair)); cl$pair <- factor(cl$pair, levels = rev(unique(cl$pair)))
  pB <- ggplot(cl, aes(v, pair, colour = w, shape = w)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = COL$ref) + geom_point(size = 1.9) +
    scale_colour_manual(values = c("Paired differences" = COL$glyco, "Loss values" = COL$ref), name = NULL) +
    scale_shape_manual(values = c("Paired differences" = 16, "Loss values" = 17), name = NULL) + scale_x_continuous(limits = c(0, 1)) +
    labs(x = "Spearman correlation", y = NULL, title = "Agreement between scores") + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 1))
  br <- rt("P2_18_B_results.tsv"); br <- br[br$run == 1, ]
  VAR <- c(all19 = "all 19 substitutions", snv7 = "single-nucleotide substitutions", strict19 = "all 19, structure-matched")
  br$est <- br$median_d; br$lo95 <- br$ci95_lo; br$hi95 <- br$ci95_hi; br$lo90 <- br$ci90_lo; br$hi90 <- br$ci90_hi; br$label <- vclass(br$label)
  br$meas <- factor(MEAS[br$score], levels = MEAS[c("AM", "ESM")])
  lev <- as.vector(t(outer(TIER, VAR, function(a, b) paste0(a, " ", U_DOT, " ", b))))
  br$y <- factor(paste0(TIER[br$tier], " ", U_DOT, " ", VAR[br$variant]), levels = rev(lev))
  bd2 <- unique(data.frame(meas = br$meas, xmin = -br$margin, xmax = br$margin))
  nullr <- geom_linerange(data = br, aes(xmin = null_min, xmax = null_max, y = y), inherit.aes = FALSE, colour = "#C9C3EE", linewidth = 2.6, na.rm = TRUE)
  pC <- forest(br, "Median difference, sequon N against background N positions", bands = bd2, lines = NULL, facet = "meas", ncol = 2, scales = "free_x",
               pre = nullr, title = "Proteome-wide: every possible substitution at sequon N",
               subtitle = "Violet bar: range of 50 null medians (random background sets); > 0: sequon N scored more damaging than background N")
  tagit((pA / ((pC | pB) + plot_layout(widths = c(2.4, 1))) + plot_layout(heights = c(1, 1.1), guides = "collect")) & theme(legend.position = "bottom", legend.box = "vertical"))
})
add_fig("F13", "aromatic_split", "Aromatic residue before the sequon (C1)", 180, 95, "P2_18_C1_results.tsv", function() {
  c1 <- rt("P2_18_C1_results.tsv")
  c1$sp <- ifelse(c1$unit == 5, paste0(c1$split, ", stricter matching"), c1$split)
  c1$est <- c1$diff_group_minus_rest; c1$lo95 <- c1$diff_ci95_lo; c1$hi95 <- c1$diff_ci95_hi; c1$lo90 <- c1$diff_ci90_lo; c1$hi90 <- c1$diff_ci90_hi
  c1$label <- desc_lab(nrow(c1)); c1$meas <- factor(MEAS[c1$measure], levels = MEAS)
  c1$sp <- paste0(c1$sp, " (n = ", c1$n_group_pairs, ")"); c1$y <- factor(c1$sp, levels = rev(unique(c1$sp)))
  forest(c1, "Difference in median loss-minus-control score, aromatic group minus the rest", band = NULL, lines = NULL, facet = "meas", ncol = 3, scales = "free_x",
         title = "Losses with an aromatic residue before the sequon (descriptive)",
         subtitle = "F/Y/W at N-2 enhances occupancy (aromatic sequon); thick 90% CI, thin 95% CI; n = loss pairs in the group")
})

## =================================================================================================================
## F14 / F15 Clinical links (Q4)
## =================================================================================================================
SPEC <- c(P_loss_only_vs_control_only = "Loss only vs control only (primary)", R_both_vs_control_only = "Loss and control vs control only",
          Q_any_loss_vs_control_only = "Any loss vs control only", M_any_loss_adj_control_carrier = "Any loss, adjusted for control carrier")
XL_HR <- "Hazard ratio for progression-free interval (log scale)"
add_fig("F14", "clinical", "Clinical links of T1 sequon losses (Q4)", 180, 170, "P2_12_primary_secondary.tsv; P2_12_per_type.tsv; P2_12_cohort.tsv", function() {
  ps <- rt("P2_12_primary_secondary.tsv")
  ps$est <- ps$HR; ps$lo95 <- ps$HR_lo; ps$hi95 <- ps$HR_hi; ps$lo90 <- exp(ps$b - qnorm(0.95) * ps$se); ps$hi90 <- exp(ps$b + qnorm(0.95) * ps$se)
  ps$label <- vclass(ps$label); ps$y <- factor(SPEC[ps$spec], levels = rev(SPEC)); ps$nt <- paste0(fmtn(ps$n_group), " (", fmtn(ps$events_group), " events)")
  pA <- forest(ps, XL_HR, band = c(0.8, 1.25), lines = NULL, ref = 1, log_x = TRUE, xlim = c(0.55, 1.8), ntext = "nt",
               title = "Cox model, 29 cancer types (strata), age, sex and log missense count",
               subtitle = "Grey band: TOST margin 0.80 to 1.25; thick 90% CI, thin 95% CI; right: carriers (events) in the exposure group")
  pt <- rt("P2_12_per_type.tsv"); pt <- pt[istrue(pt$tested), ]
  pt$est <- pt$HR; pt$lo95 <- pt$HR_lo; pt$hi95 <- pt$HR_hi; pt$label <- desc_lab(nrow(pt))
  pt$y <- factor(pt$type, levels = rev(sort(pt$type))); pt$nt <- paste0(pt$n_any_loss, " (", pt$events_any_loss, ")")
  pB <- forest(pt, XL_HR, band = c(0.8, 1.25), lines = NULL, ref = 1, log_x = TRUE, xlim = c(0.3, 4), ntext = "nt", title = "Any loss vs control only, per cancer type",
               subtitle = "Descriptive; types with enough events; right: loss carriers (events)")
  co <- rt("P2_12_cohort.tsv"); allp <- rt("P2_12_per_type.tsv")
  co <- co[co$type %in% allp$type & co$grp %in% c("control_only", "loss_only", "both") & !is.na(co$PFI.time) & co$PFI.time > 0 & !is.na(co$PFI), ]
  GL <- c(control_only = "Control only", loss_only = "Loss only", both = "Loss and control")
  fit <- survival::survfit(survival::Surv(PFI.time / 365.25, PFI) ~ grp, data = co)
  st <- rep(names(fit$strata), fit$strata)
  km <- data.frame(time = fit$time, surv = fit$surv, g = sub("^grp=", "", st), stringsAsFactors = FALSE)
  km <- rbind(data.frame(time = 0, surv = 1, g = unique(km$g)), km)
  ng <- table(co$grp); km$gl <- factor(paste0(GL[km$g], " (n = ", fmtn(ng[km$g]), ")"), levels = paste0(GL, " (n = ", fmtn(ng[names(GL)]), ")"))
  pC <- ggplot(km, aes(time, surv, colour = gl)) + geom_step(linewidth = 0.5) +
    scale_colour_manual(values = setNames(c(COL$control, COL$glyco, COL$N), levels(km$gl)), name = NULL) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) + scale_y_continuous(labels = scales::label_percent()) + scale_x_continuous(breaks = seq(0, 10, 2)) +
    labs(x = "Years", y = "Progression-free", title = "Kaplan-Meier curves (descriptive)",
         subtitle = "Pooled over the 29 types, not adjusted;\nthe test is the stratified Cox model in A") +
    theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1), legend.background = element_blank())
  tagit((pA / (pB | pC) + plot_layout(heights = c(0.55, 1), guides = "keep")))
})
add_fig("F15", "clinical_sensitivity", "Clinical analysis: sensitivity analyses, leave-one-type-out and power (Q4)", 180, 215,
        "P2_12_results_long.tsv; P2_12_leave_one_type_out.tsv; P2_12_power_run1.tsv", function() {
  rule <- function(r) vclass(ifelse(r$p < 0.05, ifelse(r$HR > 1, "Associated_worse", "Associated_better"), ifelse(r$p_tost < 0.05, "Equivalent", "Not_clear")))
  CT <- c(x_lo = "Loss only vs control only (primary)", x_any = "Any loss vs control only")
  pick <- function(r) r[(r$model == "M1_3level" & r$term == "x_lo") | (r$model == "M2_any_loss" & r$term == "x_any"), ]
  AN <- c(main = "Main", S1_strict_controls = "S1 strict controls", S2_no_hypermutated = "S2 without hypermutated", S3_no_rescued_rows = "S3 without rescued rows",
          S4_N_loss_only = "S4 N losses only", S5_ST_loss_only = "S5 S/T losses only", S6_OS_endpoint = "S6 overall survival",
          S7_plus_DLBC_KICH = "S7 plus DLBC and KICH", S8_unadjusted = "S8 unadjusted")
  rl <- pick(rt("P2_12_results_long.tsv")); rl <- rl[rl$analysis %in% names(AN), ]
  rl$est <- rl$HR; rl$lo95 <- rl$HR_lo; rl$hi95 <- rl$HR_hi; rl$lo90 <- exp(rl$b - qnorm(0.95) * rl$se); rl$hi90 <- exp(rl$b + qnorm(0.95) * rl$se)
  rl$label <- rule(rl); rl$y <- factor(AN[rl$analysis], levels = rev(AN)); rl$ct <- factor(CT[rl$term], levels = CT)
  pA <- forest(rl, XL_HR, band = c(0.8, 1.25), lines = NULL, ref = 1, log_x = TRUE, xlim = c(0.45, 1.9), facet = "ct", ncol = 2,
               title = "Sensitivity analyses", subtitle = "Labels by the pre-specified rule (Wald p < 0.05, else TOST p < 0.05 for 0.80 to 1.25); descriptive")
  lo <- pick(rt("P2_12_leave_one_type_out.tsv"))
  lo$est <- lo$HR; lo$lo95 <- lo$HR_lo; lo$hi95 <- lo$HR_hi; lo$lo90 <- exp(lo$b - qnorm(0.95) * lo$se); lo$hi90 <- exp(lo$b + qnorm(0.95) * lo$se)
  lo$label <- rule(lo); lo$y <- factor(paste("Without", lo$left_out), levels = rev(sort(unique(paste("Without", lo$left_out))))); lo$ct <- factor(CT[lo$term], levels = CT)
  pB <- forest(lo, XL_HR, band = c(0.8, 1.25), lines = NULL, ref = 1, log_x = TRUE, xlim = c(0.55, 1.6), facet = "ct", ncol = 2, title = "Leave one cancer type out")
  pw <- rt("P2_12_power_run1.tsv"); pw$sp <- factor(SPEC[pw$spec], levels = SPEC); pw$set <- paste("Simulation set", pw$set)
  pC <- ggplot(pw, aes(HR_planted, power_05, colour = sp, linetype = set)) + geom_hline(yintercept = 0.8, linetype = "22", linewidth = 0.25, colour = "grey55") +
    geom_line(linewidth = 0.5) + geom_point(size = 1.3) + scale_colour_manual(values = setNames(c(COL$glyco, COL$N, COL$ST, COL$ref), SPEC), name = NULL) +
    scale_linetype_manual(values = c("solid", "22"), name = NULL) + scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
    scale_x_continuous(breaks = c(1, 1.25, 1.5)) + labs(x = "Planted hazard ratio", y = "Share with p < 0.05", title = "Power by simulation",
         subtitle = "Weibull model fitted on the real cohort; dashed line: 80% rule at HR 1.5") +
    theme(legend.position = "right") + guides(colour = guide_legend(ncol = 1))
  tagit((pA / pC / pB + plot_layout(heights = c(0.9, 0.6, 2.1), guides = "keep")))
})

## =================================================================================================================
## F16 / F17 Gain side (B3)
## =================================================================================================================
add_fig("F16", "gain_side", "Sequon-creating substitutions (B3)", 180, 160, "P2_19_gain_results.tsv; P2_19_calibration.tsv; P2_11_foldx_gains_ddg.tsv; P2_11_foldx_ddg.tsv", function() {
  gr <- rt("P2_19_gain_results.tsv"); g <- eng(gr)
  g$y <- ufac(unit_lab(g$tier, g$group, G_TIER, G_GRP), unit_levels(G_TIER, G_GRP)); g$nt <- fmtn(g$nglyco)
  m <- g[g$config == "main" & g$scope == "pooled", ]
  pA <- forest(m, "log2 rate ratio, sequon-creating against other substitutions", xlim = c(-0.36, 0.36), ntext = "nt",
               title = "Gains pooled over 33 types (context model)", subtitle = paste0(SUB_EQ2, "\nPlausible: T1 protein, new N outside the cytoplasm; right: events"))
  cf <- g[g$scope == "pooled" & g$group == "both" & g$tier %in% c("T1", "T3"), ]; cf$y <- cfac(cf$config)
  cf$tl <- factor(c(T1 = "Plausible", T3 = "Implausible")[cf$tier], levels = c("Plausible", "Implausible"))
  pB <- forest(cf, "log2 rate ratio", facet = "tl", ncol = 2, xlim = c(-0.36, 0.2), title = "Across set-ups",
               subtitle = "New N + new S/T; position counts (last row) read\nthe implausible comparison as Inconclusive")
  cb <- rt("P2_19_calibration.tsv"); cb <- cb[cb$scope == "pooled" & cb$tier == "T1" & cb$run == min(cb$run), ]
  cmod <- c(calib_ctx = "Context model", calib_uniform = "Uniform rate")
  cb$est <- cb$log2rr; cb$lo90 <- cb$lo90; cb$label <- vclass(cb$raw_label)
  cb$y <- factor(paste0(cmod[cb$config], " ", U_DOT, " ", G_GRP[cb$group]), levels = rev(as.vector(t(outer(cmod, G_GRP, function(a, b) paste0(a, " ", U_DOT, " ", b))))))
  pC <- forest(cb, "log2 ratio, high- against low-expected background", lines = NULL, xlim = c(-0.2, 2.7), title = "Calibration of the expected rate",
               subtitle = "Passed on the pooled unit; new S/T alone 0.15\n(expected rate less exact for S/T-creating changes)")
  gd <- rt("P2_11_foldx_gains_ddg.tsv"); ld <- rt("P2_11_foldx_ddg.tsv"); ld <- ld[ld$set == "loss" & ld$role %in% c("N", "ST"), ]
  GS <- c("T1 sequon losses", "Gains, plausible", "Gains, unknown")
  dd <- rbind(data.frame(w = GS[1], v = ld$ddg_mean), data.frame(w = GS[2], v = gd$ddg_mean[gd$set == "gain_plausible"]),
              data.frame(w = GS[3], v = gd$ddg_mean[gd$set == "gain_unknown"]))
  nn <- table(dd$w); dd$w <- factor(paste0(dd$w, " (", fmtn(nn[dd$w]), ")"), levels = paste0(GS, " (", fmtn(nn[GS]), ")"))
  pD <- ggplot(dd, aes(v, colour = w)) + geom_vline(xintercept = c(1, 3), linetype = "22", linewidth = 0.25, colour = "grey60") +
    stat_ecdf(linewidth = 0.55, pad = FALSE, na.rm = TRUE) + scale_colour_manual(values = setNames(c(COL$glyco, COL$N, COL$ST), levels(dd$w)), name = NULL) +
    coord_cartesian(xlim = c(-2, 7)) + scale_y_continuous(labels = scales::label_percent()) +
    labs(x = paste0("FoldX ", U_DD, " (kcal/mol)"), y = "Cumulative share", title = paste0("FoldX ", U_DD, " of gains"),
         subtitle = "Descriptive; T1 losses for reference") + guides(colour = guide_legend(ncol = 1))
  tagit(((pA | pB) / (pC | pD) + plot_layout(heights = c(1.15, 1), guides = "collect")) & theme(legend.position = "bottom", legend.box = "vertical"))
})
add_fig("F17", "gain_shortlist", "Shortlist of sequon-creating substitutions (B3)", 180, 160,
        "P2_19_shortlist.tsv; P2_19_netnglyc_scores.tsv; P2_19_reglyco_filled.tsv; funnel steps 1 to 5 from the SqnPAN_19 log", function() {
  sl <- rt("P2_19_shortlist.tsv"); ng <- rt("P2_19_netnglyc_scores.tsv"); rg <- rt("P2_19_reglyco_filled.tsv")
  fn <- data.frame(step = c("Plausible gains (unique substitutions)", "SURFY surface protein", "AlphaFold model in one fragment",
                            paste0("pLDDT ", U_GE, " 70 and RSA ", U_GE, " 0.25 at the new N"), "No proline at N+3",
                            paste0("FoldX ", U_DD, " present and ", U_LE, " 3 kcal/mol"), "NetNGlyc potential > 0.5"),
                   n = c(3237, 1828, 1485, 714, 673, nrow(sl), sum(ng$netnglyc_potential > 0.5, na.rm = TRUE)))
  fn$step <- factor(fn$step, levels = rev(fn$step))
  pA <- ggplot(fn, aes(n, step)) + geom_col(width = 0.65, fill = SEQ_HI) + gt(aes(label = fmtn(n)), hjust = -0.12, size = TXT) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18)), labels = scales::label_comma()) +
    labs(x = "Substitutions", y = NULL, title = "Filters to the shortlist", subtitle = sprintf("%s substitutions in %s proteins after the FoldX filter", fmtn(nrow(sl)), fmtn(length(unique(sl$accession)))))
  ng$gt <- factor(ifelse(ng$gain_type == "new N", "New N", "New S/T"), levels = c("New N", "New S/T"))
  pB <- ggplot(ng, aes(netnglyc_potential, fill = gt)) + geom_histogram(binwidth = 0.025, boundary = 0, colour = "white", linewidth = 0.15, position = "stack", na.rm = TRUE) +
    geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.3) + scale_fill_manual(values = c("New N" = COL$N, "New S/T" = COL$ST), name = NULL) +
    labs(x = "NetNGlyc 1.0 potential (mutant sequence)", y = "Substitutions", title = "Predicted sequon use", subtitle = "Dashed: default threshold 0.5") +
    theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1), legend.background = element_blank())
  top <- ng[ng$id %in% rg$id, ]; top$lab <- paste(top$gene, top$mutation)
  pC <- ggplot(ng, aes(foldx_ddg, netnglyc_potential, colour = gt)) + geom_hline(yintercept = 0.5, linetype = "22", linewidth = 0.25, colour = "grey55") +
    geom_point(size = 0.7, alpha = 0.6, na.rm = TRUE) + geom_point(data = top, size = 1.6, shape = 21, fill = NA, colour = "black", stroke = 0.5) +
    ggrepel::geom_text_repel(data = top, aes(label = lab), family = FONT, size = 1.8, colour = "black", seed = 33700002, min.segment.length = 0,
                             segment.size = 0.2, max.overlaps = 50, box.padding = 0.5, force = 8, force_pull = 0.3, max.time = 3, max.iter = 50000) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.25))) +
    scale_colour_manual(values = c("New N" = COL$N, "New S/T" = COL$ST), guide = "none") +
    labs(x = paste0("FoldX ", U_DD, " of the substitution (kcal/mol)"), y = "NetNGlyc potential", title = "Shortlist", subtitle = "Black rings: top 10, one per protein")
  RS <- function(x) ifelse(x == "built_clash_free", "Clash-free", sub("^unresolved_([0-9]+)_clash_partners$", "Unresolved (\\1)", x))
  rg <- rg[order(rg$rank), ]; rg$row <- paste0(rg$rank, ". ", rg$gene, " ", rg$mutation)
  rg$row <- factor(rg$row, levels = rev(rg$row))
  CC <- c("Site", "NetNGlyc", "pLDDT", "RSA", U_DD, "Tumours", "Re-Glyco")
  tl <- do.call(rbind, list(data.frame(row = rg$row, col = CC[2], txt = f2(rg$netnglyc_potential, 3)), data.frame(row = rg$row, col = CC[1], txt = paste0("N", rg$glycosite_N, " ", rg$sequon)),
                            data.frame(row = rg$row, col = CC[3], txt = f2(rg$plddt_at_N, 1)), data.frame(row = rg$row, col = CC[4], txt = f2(rg$rsa_at_N, 2)),
                            data.frame(row = rg$row, col = CC[5], txt = f2(rg$foldx_ddg, 2)), data.frame(row = rg$row, col = CC[6], txt = as.character(rg$n_cases)),
                            data.frame(row = rg$row, col = CC[7], txt = RS(rg$clash_status))))
  tl$col <- factor(tl$col, levels = CC); tl$ok <- ifelse(tl$col == CC[7] & grepl("^Unresolved", tl$txt), "grey45", "black")
  pD <- ggplot(tl, aes(col, row, label = txt, colour = ok)) + gt(size = 1.9) + scale_colour_identity() + scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = "Top 10: mutant models and Re-Glyco", subtitle = "Re-Glyco 0.2.6, one complex N-glycan (G54258NG);\nNetNGlyc jury 9/9 for all ten; a clash-free glycan\nis not evidence of glycosylation") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(face = "bold", size = 5.6))
  tagit((pA | pB) / ((pC | pD) + plot_layout(widths = c(1, 1.55))) + plot_layout(heights = c(0.8, 1)))
})

## =================================================================================================================
## F18 / F19 Hotspots and ActiveDriver (B4)
## =================================================================================================================
add_fig("F18", "hotspots", "Single-position hotspots at T1 glycosites (B4)", 180, 155,
        "P2_20_T1_hotspots.tsv; P2_20_summary.tsv; P2_20_C3_planted.tsv; P2_20_positive_control.tsv; P2_20_background_hotspots.tsv; P2_20_C2_null.tsv", function() {
  h <- rt("P2_20_T1_hotspots.tsv"); h$hot <- grepl("^Hotspot", h$label_main); h$lab <- paste0(h$gene, " ", h$residue, h$position)
  pA <- ggplot(h, aes(expected_main, k_main)) + geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.3, colour = "grey55") +
    geom_point(aes(colour = hot), size = 0.75, alpha = 0.55) + geom_point(data = h[h$hot, ], colour = COL$up, size = 1.8) +
    ggrepel::geom_text_repel(data = h[h$hot, ], aes(label = paste0(lab, "\nFDR ", formatC(fdr_main, format = "g", digits = 2))), family = FONT, size = 1.9,
                             seed = 33700002, nudge_x = -0.6, nudge_y = 0.2, segment.size = 0.2, lineheight = 0.9) +
    scale_colour_manual(values = c(`TRUE` = COL$up, `FALSE` = "grey55"), guide = "none") + scale_x_log10() +
    scale_y_continuous(breaks = scales::breaks_pretty()) +
    labs(x = "Expected events at the position (context model, log scale)", y = "Observed events", title = "T1 glycosite positions with events",
         subtitle = sprintf("%s positions with at least one event; red: called hotspot (label Hotspot_calibration_FAILED)", fmtn(nrow(h))))
  sm <- rt("P2_20_summary.tsv")
  L <- rbind(data.frame(setup = sm$setup, m = "Hotspot (FDR < 0.05)", est = sm$OR_hotspot, lo95 = sm$OR_lo, hi95 = sm$OR_hi),
             data.frame(setup = sm$setup, m = paste0("Nominal p < 0.001 and k ", U_GE, " 3"), est = sm$OR_nominal, lo95 = sm$ORn_lo, hi95 = sm$ORn_hi))
  L$y <- factor(paste0(L$m, ", ", ifelse(L$setup == "main", "all tumours", "without hypermutated")),
                levels = rev(paste0(rep(unique(L$m), each = 2), ", ", c("all tumours", "without hypermutated")))); L$label <- desc_lab(nrow(L))
  pB <- forest(L, "Odds ratio, T1 glycosites vs background (95% CI)", legend = "none", band = NULL, lines = NULL, ref = 1, log_x = TRUE, xlim = c(0.1, 10),
               title = "T1 against background (same model)", subtitle = "Descriptive; the calibration control failed (ratio 3.11),\nso single calls are anticonservative")
  c3 <- rt("P2_20_C3_planted.tsv"); a3 <- do.call(rbind, lapply(split(c3, c3$d), function(z) data.frame(d = z$d[1], m = mean(z$recovery), lo = min(z$recovery), hi = max(z$recovery))))
  pC <- ggplot(a3, aes(factor(d), m)) + geom_col(width = 0.55, fill = SEQ_HI) + geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.4) +
    gt(aes(label = f2(m)), vjust = -0.9, size = 1.9) + scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
    labs(x = "Planted extra events at one position", y = "Recovered", title = "Planted hotspots", subtitle = "Mean and range of 20 repeats")
  pc <- rt("P2_20_positive_control.tsv"); pc$lab <- paste(pc$gene, pc$residue); pc$lab <- factor(pc$lab, levels = rev(pc$lab))
  pc$called <- ifelse(istrue(pc$called), "Called", "Not called")
  pD <- ggplot(pc, aes(k, lab, fill = called)) + geom_col(width = 0.6) + geom_point(aes(x = expected), shape = 124, size = 2.2, colour = "black") +
    scale_fill_manual(values = c(Called = COL$equiv, "Not called" = COL$nt), name = NULL) +
    labs(x = "Events (bar) and expected (tick)", y = NULL, title = "Known hotspots", subtitle = "U2AF1 S34 has no events in the masked MAFs") +
    theme(legend.position = "bottom")
  bg <- rt("P2_20_background_hotspots.tsv"); bg <- bg[order(-bg$k), ][1:min(20, nrow(bg)), ]
  bg$lab <- factor(paste0(bg$gene, " ", bg$residue, bg$position), levels = rev(unique(paste0(bg$gene, " ", bg$residue, bg$position))))
  bg$drv <- ifelse(istrue(bg$bailey_driver), "Bailey driver gene", "Other gene")
  pE <- ggplot(bg, aes(k, lab, fill = drv)) + geom_col(width = 0.65) + scale_fill_manual(values = c("Bailey driver gene" = COL$N, "Other gene" = COL$nt), name = NULL) +
    labs(x = "Events", y = NULL, title = "Top background hotspots", subtitle = "Non-sequon N, S, T positions") + theme(legend.position = "bottom")
  c2 <- rt("P2_20_C2_null.tsv"); n0 <- tapply(c2$T1_hotspots == 0, c2$setup, mean)
  pC <- pC + labs(caption = sprintf("Null control: no T1 hotspot in %s (main) and %s (without hypermutated) of the null repeats", pct(n0["main"], 0), pct(n0["nohyper"], 0))) +
    theme(plot.caption = element_text(size = 5.4, hjust = 0, colour = "grey25"))
  tagit((pA | pB) / (pC | pD | pE) + plot_layout(heights = c(1, 1.15)))
})
add_fig("F19", "activedriver", "ActiveDriver test of glycosite regions (B4)", 180, 115,
        "P2_21_activedriver_results.tsv; P2_21_controls.tsv; P2_21_concordance.tsv", function() {
  ad <- rt("P2_21_activedriver_results.tsv"); ad$enr <- grepl("^Enriched", ad$label)
  ctl <- rt("P2_21_controls.tsv")
  sh <- rbind(data.frame(w = sprintf("Real T1 sites (%s proteins)", fmtn(nrow(ad))), v = mean(ad$enr), k = "real"),
              data.frame(w = paste0("Null, random non-sequon N, repeat ", ctl$rep[ctl$control == "null"]), v = ctl$share[ctl$control == "null"], k = "null"),
              data.frame(w = "Positive control (5 proteins)", v = ctl$share[ctl$control == "positive"], k = "positive"))
  sh$w <- factor(sh$w, levels = rev(sh$w))
  pA <- ggplot(sh, aes(v, w, fill = k)) + geom_col(width = 0.6) + geom_vline(xintercept = 0.05, linetype = "22", linewidth = 0.3) +
    gt(aes(label = pct(v)), hjust = -0.15, size = 1.9) + scale_fill_manual(values = c(real = COL$glyco, null = COL$control, positive = COL$equiv), guide = "none") +
    scale_x_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.2))) +
    labs(x = "Share of proteins called (FDR < 0.05)", y = NULL, title = "Call rates", subtitle = "Dashed: null rule 5%; the real rate is below the null rate")
  ad$fdrc <- pmax(ad$fdr, 1e-12); ad$drv <- ifelse(istrue(ad$bailey_driver), "Bailey driver gene", "Other gene")
  ad$st <- factor(ifelse(ad$enr, "Enriched glycosite region", "Not called"), levels = c("Enriched glycosite region", "Not called"))
  lab <- ad[ad$enr, ]; lab <- lab[order(lab$fdr), ][1:min(12, sum(ad$enr)), ]
  if ("ERBB3" %in% ad$symbol & !("ERBB3" %in% lab$symbol)) lab <- rbind(lab, ad[ad$symbol == "ERBB3", ][1, ])
  pB <- ggplot(ad[!is.na(ad$fdr), ], aes(n_mutations, -log10(fdrc), colour = st, shape = drv)) + geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.3) +
    geom_point(size = 0.9, alpha = 0.7) + ggrepel::geom_text_repel(data = lab, aes(label = symbol), family = FONT, size = 1.8, colour = "black", seed = 33700002,
                                                                     segment.size = 0.2, max.overlaps = 40, show.legend = FALSE) +
    scale_colour_manual(values = c("Enriched glycosite region" = COL$glyco, "Not called" = "grey70"), name = NULL) +
    scale_shape_manual(values = c("Bailey driver gene" = 17, "Other gene" = 16), name = NULL) + scale_x_log10() +
    labs(x = "Missense mutations in the protein (log scale)", y = "-log10 FDR", title = "Proteins tested",
         subtitle = "Individual calls are hypotheses, not glycosite-region drivers") + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
  cc <- rt("P2_21_concordance.tsv")
  cc$a <- ifelse(istrue(cc$ActiveDriver_enriched), "ActiveDriver: called", "ActiveDriver: not called"); cc$b <- ifelse(istrue(cc$P2_20_hotspot), "Hotspot: yes", "Hotspot: no")
  pC <- ggplot(cc, aes(b, a, fill = log10(Freq + 1))) + geom_tile(colour = "white", linewidth = 0.8) + gt(aes(label = fmtn(Freq), colour = Freq > 1000), size = 2.4) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") + scale_fill_gradient(low = SEQ_LO, high = SEQ_HI, guide = "none") +
    labs(x = NULL, y = NULL, title = "Agreement with B4 hotspots", subtitle = "Proteins; the shared one is ERBB3") + theme(axis.line = element_blank(), axis.ticks = element_blank())
  tagit(((pA / pC) | pB) + plot_layout(widths = c(1, 1.35)))
})

## =================================================================================================================
## F20 CPTAC protein amount in loss carriers (B2)
## =================================================================================================================
add_fig("F20", "cptac_protein", "Protein amount in carriers of T1 sequon losses, CPTAC (B2)", 180, 215,
        "P2_22_loss_case_level.tsv; P2_22_pairs.tsv.gz; P2_22_primary.tsv; P2_22_beside.tsv; P2_22_controls.tsv; P2_22_summary.tsv; P2_22_glycopeptide.tsv", function() {
  lc <- rt("P2_22_loss_case_level.tsv")
  why <- ifelse(istrue(lc$excluded_trunc_same_gene), "Truncating mutation in the same gene (left out)",
                ifelse(!is.na(lc$why_not_usable), lc$why_not_usable, ifelse(istrue(lc$in_loss_set), "Usable: protein and RNA values", "Other")))
  wt <- as.data.frame(table(why), stringsAsFactors = FALSE)
  wt$why <- sub("^no gene key in the proteome table$", "Protein not in the global proteome", sub("^carrier has no protein value$", "Carrier has no protein value", wt$why))
  wt$why <- factor(wt$why, levels = wt$why[order(wt$Freq)])
  su <- rt("P2_22_summary.tsv"); sv <- setNames(su$value, su$item)
  pA <- ggplot(wt, aes(Freq, why)) + geom_col(width = 0.6, fill = COL$glyco) + gt(aes(label = Freq), hjust = -0.2, size = TXT) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "LOSS carrier-protein pairs", y = NULL, title = sprintf("%s LOSS pairs in CPTAC tumours", fmtn(nrow(lc))),
         subtitle = sprintf("Usable %s pairs in %s proteins;\nwithin-protein route %s pairs in %s proteins", sv["n_loss"], sv["n_loss_proteins"], sv["n_within"], sv["n_within_proteins"]))
  pr <- rt("P2_22_pairs.tsv.gz", PRO); ok <- !istrue(pr$excluded_trunc_same_gene) & !is.na(pr$y_adj)
  GL <- c(loss = "T1 sequon loss", control = "Control missense", trunc = "Truncating")
  dd <- pr[ok & ((pr$group %in% c("loss", "control") & istrue(pr$in_loss_set)) | pr$group == "trunc"), c("group", "y_adj")]
  nn <- table(dd$group); dd$g <- factor(paste0(GL[dd$group], "\n(n = ", fmtn(nn[dd$group]), ")"), levels = paste0(GL, "\n(n = ", fmtn(nn[names(GL)]), ")"))
  GC <- setNames(c(COL$glyco, COL$control, COL$trunc), levels(dd$g))
  pB <- ggplot(dd, aes(g, y_adj, fill = g)) + geom_hline(yintercept = 0, linewidth = 0.3, colour = COL$ref) +
    geom_violin(colour = NA, alpha = 0.5, scale = "width") + geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.3, fill = "white") +
    geom_point(data = dd[dd$group != "trunc", ], aes(colour = g), size = 0.45, alpha = 0.5, position = position_jitter(width = 0.12, height = 0, seed = 33700001)) +
    scale_fill_manual(values = GC, guide = "none") + scale_colour_manual(values = GC, guide = "none") + coord_cartesian(ylim = c(-3, 3)) +
    labs(x = NULL, y = "Protein level, RNA-adjusted (log2)", title = "Carrier protein levels",
         subtitle = "Analysis-set pairs; truncating = positive control;\nlog2 vs reference tumours; axis cut at -3 and 3")
  pm <- rt("P2_22_primary.tsv"); bs <- rt("P2_22_beside.tsv"); pc <- rt("P2_22_controls.tsv")
  r1 <- data.frame(y = sprintf("Primary: within-protein, RNA-adjusted (%s proteins)", pm$n_units), est = pm$D_log2, lo95 = pm$lo95, hi95 = pm$hi95, lo90 = pm$lo90, hi90 = pm$hi90,
                   lab = pm$final_label, stringsAsFactors = FALSE)
  BL <- c("other estimator" = "Pooled estimator", "unadjusted outcome" = "Not RNA-adjusted", "nohyper" = "Without hypermutated",
          "mutations in all used files only" = "Mutations in all used files", "mechanism N" = "N losses only", "mechanism ST" = "S/T losses only", "mechanism XP" = "X to P losses only")
  b1 <- bs[bs$analysis %in% names(BL), ]
  r2 <- data.frame(y = paste0(BL[b1$analysis], " (", b1$n_loss, " loss pairs)"), est = b1$D_log2, lo95 = b1$lo95, hi95 = b1$hi95, lo90 = b1$lo90, hi90 = b1$hi90,
                   lab = b1$descriptive_pattern_not_a_result_label, stringsAsFactors = FALSE)
  p1 <- pc[grepl("^positive control", pc$control), ]
  r3 <- data.frame(y = ifelse(grepl("y_unadj", p1$control), "Positive control: truncating, not adjusted (pre-specified)", "Positive control: truncating, RNA-adjusted"),
                   est = as.numeric(p1$share_Lower_or_D), lo95 = as.numeric(p1$share_Equivalent_or_lo95), hi95 = as.numeric(p1$share_Not_clear_or_hi95),
                   lo90 = NA_real_, hi90 = NA_real_, lab = p1$median_D_or_label, stringsAsFactors = FALSE)
  ee <- rbind(r1, r2, r3); ee$label <- vclass(ee$lab); ee$y <- factor(ee$y, levels = rev(ee$y))
  pC <- forest(ee, "Difference in protein level, loss (or truncating) minus control carriers (log2)", xlim = c(-0.45, 0.45),
               title = "Protein amount: estimates", subtitle = paste0(SUB_EQ, "\nPrimary label Equivalent_positive_control_FAILED: truncating mutations lower protein (ratio 0.91)\nbut by less than the 0.32 threshold; rows beside the primary: descriptive patterns"))
  pl <- pc[!grepl("^positive control", pc$control), ]
  PL <- rbind(data.frame(sc = pl$control, w = "Lower", v = as.numeric(pl$share_Lower_or_D)), data.frame(sc = pl$control, w = "Equivalent", v = as.numeric(pl$share_Equivalent_or_lo95)),
              data.frame(sc = pl$control, w = "Not clear", v = as.numeric(pl$share_Not_clear_or_hi95)))
  PL$sc <- factor(sub("^planted ", "Planted ", sub("^null$", "Null", PL$sc)), levels = rev(sub("^planted ", "Planted ", sub("^null$", "Null", pl$control))))
  PL$w <- factor(PL$w, levels = c("Lower", "Equivalent", "Not clear"))
  pD <- ggplot(PL, aes(v, sc, fill = w)) + geom_col(width = 0.6, colour = "white", linewidth = 0.25, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = c(Lower = COL$down, Equivalent = COL$equiv, "Not clear" = COL$inconc), name = NULL) +
    scale_x_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
    labs(x = "Share of pseudo data sets", y = NULL, title = "Null and planted controls", subtitle = "500 null and 200 planted sets (log2 shifts)") + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1))
  b2 <- bs[grepl("^cancer type", bs$analysis), ]
  ct <- data.frame(y = paste0(sub("^cancer type ", "", b2$analysis), " (", b2$n_loss, ")"), est = b2$D_log2, lo95 = b2$lo95, hi95 = b2$hi95, lo90 = b2$lo90, hi90 = b2$hi90)
  ct$label <- desc_lab(nrow(ct)); ct$y <- factor(ct$y, levels = rev(sort(ct$y)))
  pE <- forest(ct, "Difference (log2)", xlim = c(-0.9, 0.9), lines = NULL, legend = "none", title = "Per CPTAC cancer type", subtitle = "Descriptive; (loss pairs)")
  gp <- rt("P2_22_glycopeptide.tsv"); gp <- gp[!is.na(gp$observed_log2_vs_reference_median), ]
  pF <- ggplot(gp, aes(expected_log2_by_design, observed_log2_vs_reference_median)) + geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.3, colour = "grey55") +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = COL$ref) + geom_point(colour = COL$glyco, size = 1.3, alpha = 0.8) +
    labs(x = "Expected by design (log2)", y = "Observed glycopeptide (log2)", title = "Glycopeptides at lost sites",
         subtitle = sprintf("%s rows at %s sites; descriptive,\nnever read as glycosylation\nloss or retention", nrow(gp), sv["glyco_sites_with_carrier_value"]))
  tagit((pA | pB) / pC / (pD | pE | pF) + plot_layout(heights = c(1, 1.15, 0.95)))
})

## =================================================================================================================
## F21 dNdScv
## =================================================================================================================
add_fig("F21", "dndscv", "dNdScv check of selection at T1 glycosites", 180, 170,
        "P2_23_results.tsv; P2_23_leave_one_out.tsv; P2_23_absolute.tsv; P2_23_planted.tsv; P2_23_capability.tsv", function() {
  dn <- rt("P2_23_results.tsv")
  m <- dn[dn$setup == "main" & dn$mapping == "strict" & dn$run == 1, ]
  m$est <- m$log2rr; m$label <- vclass(m$final_label); m$y <- ufac(unit_lab(m$tier, m$group)); m$nt <- paste0(fmtn(m$ev_g), " vs ", fmtn(m$ev_b))
  pA <- forest(m, XL_RR, xlim = c(-0.36, 0.36), ntext = "nt", title = "dNdScv expected rates, main run",
               subtitle = paste0(SUB_EQ2, "\nRight: glycosite vs background events; T1 N and S/T alone:\nInconclusive, test not sensitive"))
  SU <- c("main strict 1" = "Main (strict mapping)", "main strict 3" = "Main, second seed", "main window 1" = "Window mapping", "nocaps strict 1" = "No mutation caps",
          "nohyper strict 1" = "Without hypermutated", "nodrivers strict 1" = "Without driver genes")
  s1 <- dn[dn$tier == "T1" & dn$group == "both", ]; s1$k <- paste(s1$setup, s1$mapping, s1$run); s1 <- s1[s1$k %in% names(SU), ]
  s1$y <- SU[s1$k]; s1$lab <- ifelse(is.na(s1$final_label), s1$raw_label, s1$final_label)
  lo <- rt("P2_23_leave_one_out.tsv"); lo <- lo[lo$tier == "T1" & lo$group == "both" & !grepl("^\\(none", lo$left_out), ]
  lo$y <- paste("Without", sub("^TCGA-", "", lo$left_out)); lo$lab <- lo$raw_label
  cols <- c("y", "log2rr", "lo95", "hi95", "lo90", "hi90", "lab")
  sb <- rbind(s1[, cols], lo[, cols]); sb$est <- sb$log2rr; sb$label <- vclass(sb$lab); sb$y <- factor(sb$y, levels = rev(unique(sb$y)))
  pB <- forest(sb, XL_RR, xlim = c(-0.22, 0.22), lines = NULL, legend = "none", title = "T1 N + S/T across set-ups", subtitle = "Leave-one-out of the five largest projects")
  ab <- rt("P2_23_absolute.tsv"); ab$g <- factor(GRP[ab$group], levels = rev(GRP)); ab$cl <- factor(ifelse(grepl("^T1", ab$class), "T1 glycosite", "Background, same genes"),
                                                                                                   levels = c("T1 glycosite", "Background, same genes"))
  ab$yn <- as.numeric(ab$g) + ifelse(ab$cl == "T1 glycosite", 0.14, -0.14)
  pC <- ggplot(ab, aes(dnds_cal, yn, colour = cl)) + geom_vline(xintercept = 1, linewidth = 0.3, colour = COL$ref) +
    geom_linerange(aes(xmin = lo95_cal, xmax = hi95_cal), linewidth = 0.6) + geom_point(size = 1.6) +
    scale_y_continuous(breaks = seq_along(levels(ab$g)), labels = levels(ab$g), expand = expansion(add = 0.4)) +
    scale_colour_manual(values = c("T1 glycosite" = COL$glyco, "Background, same genes" = COL$control), name = NULL) +
    labs(x = "dN/dS, divided by the synonymous calibration 1.059 (95% CI)", y = NULL, title = "Absolute dN/dS", subtitle = "Descriptive, no label") +
    theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 1))
  pl <- rt("P2_23_planted.tsv"); pl$g <- factor(c(N = "N", ST = "S/T", both = "N + S/T")[pl$group], levels = c("N", "S/T", "N + S/T")); pl$label <- vclass(pl$raw_label)
  pl$r <- factor(formatC(pl$ratio, format = "fg"), levels = formatC(sort(unique(pl$ratio)), format = "fg"))
  pD <- ggplot(pl, aes(r, fill = label)) + geom_bar(width = 0.7, colour = "white", linewidth = 0.2) + facet_wrap(~g, nrow = 1) +
    scale_fill_manual(values = VCOL, name = NULL) +
    labs(x = "Planted rate ratio", y = "Planted sets (of 20)", title = "Planted effects, T1") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(fill = guide_legend(ncol = 2))
  cp <- rt("P2_23_capability.tsv"); cp$y <- ufac(unit_lab(cp$tier, cp$group)); cp$al <- ifelse(istrue(cp$equivalent_allowed), "Equivalent allowed", "Equivalent not allowed")
  pE <- ggplot(cp, aes(n_equivalent, y, fill = al)) + geom_col(width = 0.6) + geom_vline(xintercept = 16, linetype = "22", linewidth = 0.3) +
    scale_fill_manual(values = c("Equivalent allowed" = COL$equiv, "Equivalent not allowed" = COL$nt), name = NULL) + scale_x_continuous(limits = c(0, 20)) +
    labs(x = "Null sets labelled Equivalent (of 20)", y = NULL, title = "Null control", subtitle = paste0("Rule: ", U_GE, " 16 of 20;\nno Changed call in any null set")) +
    theme(legend.position = "bottom") + guides(fill = guide_legend(ncol = 1))
  tagit((pA | pB) / (pC | pD | pE) + plot_layout(heights = c(1, 1)))
})

## =================================================================================================================
## F22 Synthesis on one scale
## =================================================================================================================
add_fig("F22", "synthesis", "All headline estimates scaled to their equivalence margins", 180, 125,
        "P2_09, P2_14, P2_15, P2_17, P2_19, P2_23, P2_11, P2_18, P2_22, P2_12 result tables", function() {
  one <- function(d) d[d$config == "main" & d$scope == "pooled" & d$tier == "T1" & d$group == "both", ]
  rows <- list()
  addr <- function(block, name, est, lo, hi, margin, lab, unit) rows[[length(rows) + 1]] <<- data.frame(block = block, y = name, est = est / margin, lo90 = lo / margin, hi90 = hi / margin,
      lab = lab, txt = sprintf("%s %s (%s to %s); margin %s", unit, f2(est, 3), f2(lo, 3), f2(hi, 3), trimws(formatC(margin, format = "fg", digits = 3))), stringsAsFactors = FALSE)
  B1 <- "Selection (log2 rate ratio)"; B2 <- "Structure and predicted effect"; B3 <- "Protein amount and survival"
  for (z in list(list("P2_09_Q2_results.tsv", "Q2 position-count model (primary), T1 losses"), list("P2_14_Q2ctx_results.tsv", "A1 context-aware model, T1 losses"))) {
    r <- one(rt(z[[1]])); addr(B1, z[[2]], r$log2rr, r$lo90, r$hi90, EQ, r$final_label, "log2RR") }
  d <- rt("P2_23_results.tsv"); r <- d[d$setup == "main" & d$mapping == "strict" & d$run == 1 & d$tier == "T1" & d$group == "both", ]
  addr(B1, "dNdScv rates, T1 losses", r$log2rr, r$lo90, r$hi90, EQ, r$final_label, "log2RR")
  r <- one(rt("P2_19_gain_results.tsv")); addr(B1, "B3 sequon gains, plausible", r$log2rr, r$lo90, r$hi90, EQ, r$final_label, "log2RR")
  r <- one(rt("P2_15_germline_results.tsv")); addr(B1, "A2 germline, T1 glycosites (depletion expected)", r$log2rr, r$lo90, r$hi90, EQ, r$final_label, "log2RR")
  r <- one(rt("P2_17_constrained_results.tsv")); addr(B1, "Germline metal-binding Cys + His (positive control)", r$log2rr, r$lo90, r$hi90, EQ, r$final_label, "log2RR")
  r <- rt("P2_11_comparison_results.tsv"); r <- r[r$unit == 1, ]; addr(B2, paste0("Q3 FoldX ", U_DD, ", loss minus control"), r$median_d, r$ci90_lo, r$ci90_hi, 0.5, r$label, "kcal/mol")
  a <- rt("P2_18_A_results.tsv"); a <- a[a$run == 1 & a$unit == 1, ]
  for (k in c("AM", "ESM")) { r <- a[a$measure == k, ]; addr(B2, paste0("B1 ", c(AM = "AlphaMissense", ESM = "ESM1b")[k], ", loss minus control"), r$median_d, r$ci90_lo, r$ci90_hi, r$margin, r$label, "difference") }
  r <- rt("P2_22_primary.tsv"); addr(B3, "B2 CPTAC protein, loss minus control carriers", r$D_log2, r$lo90, r$hi90, EQ, r$final_label, "log2")
  ps <- rt("P2_12_primary_secondary.tsv")
  for (k in c("P_loss_only_vs_control_only", "Q_any_loss_vs_control_only")) {
    r <- ps[ps$spec == k, ]; z <- qnorm(0.95)
    rows[[length(rows) + 1]] <- data.frame(block = B3, y = paste0("Q4 PFI, ", c(P_loss_only_vs_control_only = "loss only vs control only (primary)", Q_any_loss_vs_control_only = "any loss vs control only")[k]),
      est = r$b / log(1.25), lo90 = (r$b - z * r$se) / log(1.25), hi90 = (r$b + z * r$se) / log(1.25), lab = r$label,
      txt = sprintf("HR %s (%s to %s); margin 0.80 to 1.25", f2(r$HR), f2(exp(r$b - z * r$se)), f2(exp(r$b + z * r$se))), stringsAsFactors = FALSE)
  }
  s <- do.call(rbind, rows); s$label <- vclass(s$lab); s$block <- factor(s$block, levels = c(B1, B2, B3)); s$y <- factor(s$y, levels = rev(unique(s$y)))
  p1 <- forest(s, "Estimate divided by its equivalence margin (thick bar: 90% CI)", band = c(-1, 1), lines = NULL, xlim = c(-4.3, 1.6),
               title = "Every headline estimate on one scale", subtitle = "Grey band: inside the pre-specified equivalence margin of each analysis (labels as written by the analysis scripts)") +
    facet_grid(block ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(16)) + theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 5.4))
  p2 <- ggplot(s, aes(0, y, label = txt)) + gt(hjust = 0, size = 1.7) + facet_grid(block ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(16)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + theme_void(base_family = FONT) + theme(strip.text = element_blank(), plot.margin = margin(3, 2, 3, 2))
  (p1 | p2) + plot_layout(widths = c(2, 1.15))
})

## =================================================================================================================
## F23 Label map of every set-up; F24 capability map
## =================================================================================================================
add_fig("F23", "sensitivity_map", "Labels of every pooled unit in every set-up", 180, 200,
        "P2_09_Q2_results.tsv; P2_14_Q2ctx_results.tsv; P2_15_germline_results.tsv; P2_17_constrained_results.tsv; P2_19_gain_results.tsv; P2_23_results.tsv", function() {
  short <- c(T1 = "T1", T2 = "T2", T3 = "T3"); sg <- c(N = "N", ST = "S/T", both = "N+S/T")
  mk <- function(file, an) { e <- eng(rt(file)); e <- e[e$scope == "pooled", ]
    data.frame(analysis = an, x = CFG[e$config], y = paste(short[e$tier], sg[e$group]), est = e$est, label = e$label, stringsAsFactors = FALSE) }
  L <- rbind(mk("P2_09_Q2_results.tsv", "Q2 position-count model (tumour)"), mk("P2_14_Q2ctx_results.tsv", "A1 context-aware model (tumour)"),
             mk("P2_15_germline_results.tsv", "A2 germline glycosites"), mk("P2_17_constrained_results.tsv", "Germline constrained residues (N = Cys, S/T = His)"),
             mk("P2_19_gain_results.tsv", "B3 sequon gains (T1 to T3 = plausibility)"))
  dn <- rt("P2_23_results.tsv"); SU <- c("main strict 1" = "Main (strict)", "main strict 3" = "Second seed", "main window 1" = "Window mapping", "nocaps strict 1" = "No caps",
                                           "nohyper strict 1" = "Without hypermutated", "nodrivers strict 1" = "Without drivers")
  dn$k <- paste(dn$setup, dn$mapping, dn$run); dn <- dn[dn$k %in% names(SU), ]
  L <- rbind(L, data.frame(analysis = "dNdScv", x = SU[dn$k], y = paste(short[dn$tier], sg[dn$group]), est = dn$log2rr,
                           label = vclass(ifelse(is.na(dn$final_label), dn$raw_label, dn$final_label)), stringsAsFactors = FALSE))
  L$analysis <- factor(L$analysis, levels = unique(L$analysis))
  L$x <- factor(L$x, levels = unique(c(CFG, SU)))
  L$y <- factor(L$y, levels = rev(as.vector(t(outer(short, sg, paste)))))
  ggplot(L, aes(x, y, fill = label)) + geom_tile(colour = "white", linewidth = 0.5) +
    gt(aes(label = f2(est), colour = label %in% VLEV[1:3]), size = 1.65) + scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_manual(values = VCOL, limits = VLEV, drop = FALSE, name = NULL) + facet_wrap(~analysis, ncol = 2, scales = "free_x") +
    labs(x = NULL, y = NULL, title = "Labels of every pooled unit in every set-up", subtitle = "Numbers: log2 rate ratio; colour: label after the control rules (set-ups beside the main run carry raw labels where no final label exists)") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1, size = 5.4), legend.position = "bottom")
})
add_fig("F24", "capability_map", "What each test can detect: planted effects in all pooled units", 180, 175,
        "P2_08_unit_capability.tsv; P2_14_unit_capability.tsv; P2_15_unit_capability.tsv; P2_17_unit_capability.tsv; P2_19_unit_capability.tsv; P2_23_planted.tsv", function() {
  L <- rbind(cap_long(rt("P2_08_unit_capability.tsv"), "Q2 position-count"), cap_long(rt("P2_14_unit_capability.tsv"), "A1 context-aware"),
             cap_long(rt("P2_15_unit_capability.tsv"), "A2 germline glycosites"), cap_long(rt("P2_17_unit_capability.tsv"), "Cys/His (germline)", CS_TIER, CS_GRP),
             cap_long(rt("P2_19_unit_capability.tsv"), "B3 sequon gains", G_TIER, G_GRP))
  pl <- rt("P2_23_planted.tsv"); want <- c("0.67" = "Depleted", "0.8" = "Depleted", "1.25" = "Enriched", "1.5" = "Enriched", "1" = "Equivalent")
  key <- formatC(pl$ratio, format = "fg"); hit <- pl$raw_label == want[key]
  sc <- setNames(SCEN, c("0.67", "0.8", "1.25", "1.5", "1"))
  ag <- aggregate(hit ~ group + key, data = data.frame(group = pl$group, key = key, hit = hit), FUN = sum)
  L <- rbind(L, data.frame(analysis = "dNdScv (T1 only)", y = unit_lab("T1", ag$group), scen = sc[ag$key], n = ag$hit, stringsAsFactors = FALSE))
  L$analysis <- factor(L$analysis, levels = unique(L$analysis))
  captile(L, "What each test can detect", "Number of 20 planted data sets called as expected (run 1); main set-up, pooled units", ncol = 3) + theme(legend.position = "bottom")
})

## =================================================================================================================
## Output
## =================================================================================================================
dev_print <- function(open_fun, p) { open_fun(); tryCatch(print(p), finally = grDevices::dev.off()); invisible(TRUE) }
save_fig <- function(p, base, w, h) {
  if (MODE == "draft") {
    f <- file.path(OUT, paste0(base, ".png"))
    dev_print(function() ragg::agg_png(f, width = w, height = h, units = "mm", res = 300, background = "white"), p); return(f)
  }
  f1 <- file.path(OUT, "pdf", paste0(base, ".pdf")); f2p <- file.path(OUT, "tiff", paste0(base, ".tiff")); f3 <- file.path(OUT, "png", paste0(base, ".png"))
  dev_print(function() grDevices::cairo_pdf(f1, width = w / 25.4, height = h / 25.4, family = FONT), p)
  dev_print(function() ragg::agg_tiff(f2p, width = w, height = h, units = "mm", res = 600, compression = "lzw", background = "white"), p)
  dev_print(function() ragg::agg_png(f3, width = w, height = h, units = "mm", res = 300, background = "white"), p)
  paste(f1, f2p, f3, sep = " | ")
}

idx <- list(); n_ok <- 0L
for (f in FIGS) {
  if (length(ONLY) && !(f$id %in% ONLY)) next
  t0 <- Sys.time(); nw <- 0L; base <- paste0("P2_", f$id, "_", f$name)
  st <- tryCatch(withCallingHandlers({ p <- f$fun(); save_fig(p, base, f$w, f$h); "made" },
                                     warning = function(w) { nw <<- nw + 1L; lg("  warning in ", f$id, ": ", conditionMessage(w)); invokeRestart("muffleWarning") }),
                 error = function(e) { while (grDevices::dev.cur() > 1) grDevices::dev.off(); lg("  ERROR in ", f$id, ": ", conditionMessage(e)); paste0("ERROR: ", conditionMessage(e)) })
  if (identical(st, "made")) n_ok <- n_ok + 1L
  lg(f$id, " ", f$name, ": ", st, ", ", nw, " warnings, ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
  idx[[length(idx) + 1]] <- data.frame(id = f$id, file = base, title = f$title, width_mm = f$w, height_mm = f$h, status = st, warnings = nw, sources = f$src, stringsAsFactors = FALSE)
}
if (length(idx)) write.table(do.call(rbind, idx), file.path(OUT, "P2_24_figure_index.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
lg("P2_24 finished: ", n_ok, " of ", length(idx), " figures made in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min; folder ", OUT)
