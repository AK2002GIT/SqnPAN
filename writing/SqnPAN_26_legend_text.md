# SqnPAN legends, alt text, Table 1 and Table 2 texts (Tables SM2 and SM3 of the Supplementary Methods) and supplementary table titles
# Read by scripts/SqnPAN_26_tables_legends_workbook.py. Line prefixes: TITLE, LEGEND, ALT, HEADER, ROW, FOOT, T2ROW, STAB.

### Figure_1
TITLE Study design, sequon evidence tiers and the landscape of sequon-changing mutations
LEGEND (A) The N-glycosylation sequon N-X-S/T (X is not proline) and the three classes of sequon loss, with the number of Tier 1 events in each (N to another residue, 1,009; X to proline, 28; S/T to a residue other than S or T, 1,721). A gain is a substitution that creates a new N or a new S/T and completes a new sequon. (B) Sequons in the reviewed human proteome (UniProt release 2026_03) by evidence tier. Tier 1 (experimental) sites carry experimental UniProt evidence or are N-GlycositeAtlas sites (13,247), Tier 2 (annotated) sites carry other UniProt N-glycosylation annotation (8,087) and Tier 3 sequons are motif only (38,218). (C) Data sources and analyses. The primary analyses were pre-specified. Each analysis added after the primary analysis has its own controls and changed no primary label. (D) Sequon losses by tier and gains by region class among 1,460,792 missense mutations in 10,190 TCGA tumours of 33 cancer types, with the number of tumours carrying each class. (E) Position of the lost residue within the sequon, in all tumours and without hypermutated tumours (more than 1,000 missense mutations); counts are printed inside bars that hold more than 7% of losses. TCGA, The Cancer Genome Atlas; CPTAC, Clinical Proteomic Tumor Analysis Consortium; GDC, Genomic Data Commons; MAF, mutation annotation format; PFI, progression-free interval; TCGA-CDR, TCGA Clinical Data Resource; ΔΔG, change in folding free energy.
ALT Five-panel overview. A sequon diagram shows the three ways a missense substitution can destroy an N-X-S/T sequon, a bar chart gives the number of sequons in three evidence tiers, a flow chart lists the data and analyses, and bar charts count losses and gains and show where in the sequon the lost residue lies.

### Figure_2
TITLE Selection at sequon positions estimated by three frameworks
LEGEND (A) log2 rate ratio of sequon-destroying substitutions at glycosite positions against other N, S and T positions of the same proteins, by tier and lost position, pooled over 33 TCGA cancer types, in the position-count model (primary) and the context-aware trinucleotide model. (B) The same ratio with expected rates from dNdScv (main run); numbers on the right are glycosite against background events, and Tier 1 N and S/T alone are Inconclusive because the test is not sensitive there. (C) dNdScv estimates for Tier 1 N + S/T losses across set-ups and with each of the five largest TCGA projects left out. (D) Selection estimates divided by their equivalence margins, with germline Tier 1 glycosites and the constrained metal-binding Cys + His positive control beside them; all headline estimates are shown in Supplementary Fig. S17. The grey band is the equivalence region (±0.14 log2, about 10%), dashed lines mark the change threshold (±0.32 log2), thick bars are 90% and thin bars 95% confidence intervals (CIs), and colour gives the label by the pre-specified rule. The context-aware model was added after the primary analysis, and dNdScv was run after the first two results were known (not blind).
ALT Forest plots of log2 rate ratios for sequon losses at glycosites against other positions in three selection models. All Tier 1 estimates for N plus S/T lie inside the grey equivalence band of ±0.14, and the germline Cys and His positive control lies well below it.

### Figure_3
TITLE Germline comparison and positive control at constrained Cys/His
LEGEND (A) The same test code applied to tumour somatic missense mutations and to germline missense variants (gnomAD, TOPMed, ExAC, 1000 Genomes and ESP) at glycosite positions, by tier and lost position. Germline variation is not depleted at Tier 1 sequons (log2 rate ratio 0.003, 90% CI −0.014 to 0.021), so the tumour label is read with the note that the germline biology control failed. (B) Positive control of the method, with germline variants at UniProt metal-binding, active-site and disulfide Cys and His residues tested against other Cys and His positions; the decision unit, metal-binding Cys + His, shows the expected depletion (−0.526, −0.589 to −0.465). (C) Germline Tier 1 glycosites across set-ups (main run, structure-matched strata, position counts in the context-model protein set and variants reported by two or more resources). (D) Metal-binding Cys and His across the same set-ups without structure-matched strata. Bands, lines and bars as in Fig. 2. This analysis was added after the primary analysis. gnomAD, Genome Aggregation Database; TOPMed, Trans-Omics for Precision Medicine; ExAC, Exome Aggregation Consortium; ESP, NHLBI Exome Sequencing Project.
ALT Forest plots comparing tumour and germline variants with the same test. Germline Tier 1 glycosites sit inside the equivalence band, while germline metal-binding Cys and His residues are clearly depleted in every set-up.

### Figure_4
TITLE Mutational processes behind sequon losses and the cross-check with Fan et al. 2018
LEGEND (A) Shares of signature groups among Tier 1 losses by lost position, all sequon losses, background N and S/T positions, gains and all missense mutations, pooled over cancer types; each mutation is shared among signatures by its per-type channel probabilities, and n is the number of mutations. (B) Share of mutations at CpG C>T, with 95% CI. (C) Differences in share between sets in percentage points, with 95% CI (descriptive); the position-3 excess of losses follows CpG C>T, UV and APOBEC exposure, the processes that act on background S/T. (D) Counts of this study beside those reported by Fan et al. 2018 {Fan, 2018}; LOG and GOG are the terms of Fan et al. for loss and gain of N-glycosylation sequons, and data sources and sequon definitions differ, so this is not a test. (E) Ratio of S/T to N losses observed in all tumours and without hypermutated tumours, and expected from the background rates, by tier. The signature analysis was added after the primary analysis, and the cross-check was done after the main results were known (not blind). SBS, single-base substitution signature; MMR, mismatch repair; POLE, DNA polymerase epsilon; APOBEC, apolipoprotein B mRNA-editing catalytic polypeptide-like cytidine deaminases; UV, ultraviolet light.
ALT Stacked bars of mutational signature shares, a CpG C>T share chart, a forest plot of share differences, a log-scale bar chart of counts beside Fan et al. 2018 and a dot plot of S/T-to-N loss ratios. S/T losses carry more CpG C>T, UV and APOBEC mutations, and the observed S/T excess matches the ratio expected from background rates.

### Figure_5
TITLE Predicted structural and functional impact of sequon losses
LEGEND (A) Cumulative distributions of FoldX ΔΔG for Tier 1 N (732) and S/T (1,134) losses and their matched control substitutions (first control draw); dashed lines mark 1 and 3 kcal/mol, and ΔΔG above 0 is destabilising. (B) Median difference, loss minus matched controls, for FoldX ΔΔG, AlphaMissense pathogenicity and ESM1b damage (negative log-likelihood ratio), in the primary unit (N + S/T losses, 1,443 pairs) and five set-ups; the grey band is the pre-specified margin of each score (FoldX 0.5 kcal/mol, AlphaMissense 0.082, ESM1b 0.71), and values above 0 mean that the loss is predicted more damaging. (C) Proteome-wide median score difference between every possible substitution at sequon N and at background N positions, by tier, for all 19 substitutions, single-nucleotide substitutions only and structure-matched sets; the violet bar is the range of 50 null medians from random background sets. (D) Spearman correlation between scores, for loss values and for paired differences. Thick bars are 90% and thin bars 95% CIs. FoldX is the primary structural comparison; AlphaMissense and ESM1b were added after the primary analysis. LLR, log-likelihood ratio.
ALT Cumulative ΔΔG curves for losses and matched controls overlap, forest plots of loss-minus-control differences for FoldX, AlphaMissense and ESM1b lie inside their margins, a proteome-wide forest plot compares sequon N with background N, and a small bar chart shows the correlation between scores.

### Figure_6
TITLE Protein amount in CPTAC carriers of sequon losses and clinical links
LEGEND (A) The 263 loss carrier-protein pairs in CPTAC tumours and the reasons pairs were left out (truncating mutation in the same gene, no protein value for the carrier, protein not in the global proteome); 141 pairs in 133 proteins were usable and 86 pairs in 81 proteins entered the within-protein analysis. (B) RNA-adjusted protein levels (log2 against reference tumours) of loss carriers (141), control missense carriers (245) and truncating-mutation carriers (14,026, positive control); the axis is cut at −3 and 3. (C) Difference in protein level, loss (or truncating) carriers minus control missense carriers (log2). The primary within-protein, RNA-adjusted estimate is Equivalent with the positive control failed, because truncating mutations lowered protein (ratio 0.91) by less than the 0.32 change threshold; rows beside the primary are descriptive. (D) Share of 500 null and 200 planted pseudo data sets (log2 shifts of −0.32, −0.5 and −1.0) given each label. (E) Kaplan-Meier curves of progression-free interval pooled over 29 TCGA cancer types for control-only (4,711), loss-only (183) and loss-and-control (1,207) carriers, not adjusted (descriptive). (F) Hazard ratios from the Cox model stratified by cancer type and adjusted for age, sex and log missense count; the primary contrast, loss only against control only (183 carriers, 60 events), is Not clear, the grey band is the equivalence margin of 0.80 to 1.25, and numbers on the right are carriers (events) in the exposure group. Thick bars are 90% and thin bars 95% CIs. The protein analysis was added after the primary analysis. HR, hazard ratio.
ALT Bar chart of carrier-protein pairs, box plots of protein levels by carrier group, a forest plot of protein differences inside the equivalence band with the truncating control below zero, bars of control outcomes, overlapping Kaplan-Meier curves and a forest plot of hazard ratios near 1.

### Figure_7
TITLE Sequon-creating substitutions, shortlist and glycan modelling
LEGEND (A) log2 rate ratio of sequon-creating against other substitutions, pooled over 33 cancer types (context-aware model), for gains that create a new N, a new S/T or either, in plausible regions (a new N outside the cytoplasm in a protein with a Tier 1 site), plausible or unknown regions and implausible regions (comparison); numbers on the right are events. (B) New N + new S/T gains across set-ups; position counts in the context-model protein set read the implausible comparison as Inconclusive. (C) Filters from plausible gains (unique substitutions) to the shortlist of 614 substitutions in 426 proteins. (D) NetNGlyc 1.0 potential of the new sequon in the mutant sequence, for new N and new S/T substitutions; the dashed line is the default threshold of 0.5. (E) FoldX ΔΔG against NetNGlyc potential for the shortlist; black rings mark the top 10, one per protein. (F) The top 10 with mutant models and Re-Glyco glycan modelling (Re-Glyco 0.2.6, one complex N-glycan, GlyTouCan G54258NG), with site and sequon, NetNGlyc potential, pLDDT, RSA, ΔΔG, tumours carrying the substitution and the Re-Glyco result; eight glycans were placed without clashes and two were unresolved (number of clashing protein residues in brackets), and the NetNGlyc jury agreement was 9 of 9 for all ten. A clash-free glycan is not evidence of glycosylation, and shortlist entries are candidates, not predictions of glycosylation. Bands, lines and bars in A and B as in Fig. 2. This analysis was added after the primary analysis. pLDDT, predicted local distance difference test; RSA, relative solvent accessibility; SURFY, in silico human surfaceome.
ALT Forest plots show gains in plausible regions inside the equivalence band, a bar chart filters gains to 614 shortlisted substitutions, a histogram shows NetNGlyc potentials, a scatter plot compares FoldX ΔΔG with NetNGlyc potential, and a table summarises the top 10, eight of them with clash-free glycan models.

### Figure_8
TITLE Single-position hotspots at glycosites and ActiveDriver test of glycosite regions
LEGEND (A) Observed against expected events (context model, log scale) at the 2,407 Tier 1 glycosite positions with at least one event; red marks the single position called a hotspot, ERBB3 T355 (FDR 0.0014; label Hotspot, calibration failed). (B) Odds ratio of hotspot calls at Tier 1 glycosites against background N, S and T positions with the same model, for FDR below 0.05 and for nominal p below 0.001 with at least three events, in all tumours and without hypermutated tumours (descriptive); the calibration control failed (ratio 3.11), so single calls are anticonservative. (C) Recovery of planted hotspots (extra events at one position; mean and range of 20 repeats); in the null control no Tier 1 hotspot was called in 100% (main) and 99% (without hypermutated tumours) of the repeats. (D) Events and expected events at known hotspots; U2AF1 S34 has no events in the masked MAF files. (E) ActiveDriver call rates (FDR below 0.05) for the 5,246 proteins with Tier 1 sites, three null sets of random non-sequon N positions and a positive control of five proteins; the dashed line is the 5% null rule, and the real rate (2.7%) is below the null rate. (F) Agreement between ActiveDriver calls and single-position hotspots; the one shared protein is ERBB3. (G) Proteins tested, missense mutations in the protein against −log10 FDR; single calls are hypotheses, not glycosite-region drivers. Both analyses were added after the primary analysis. FDR, false discovery rate (Benjamini-Hochberg).
ALT A scatter of observed against expected events at Tier 1 positions marks one hotspot, ERBB3 T355. Odds ratio, planted-control and known-hotspot panels follow, then bars show the ActiveDriver call rate below the null rate, a tile chart of agreement with hotspots and a scatter of the proteins tested.

### Figure_S1
TITLE Tier 1 sequon losses by cancer type and the most recurrent losses
LEGEND (A) Share of TCGA tumours with at least one Tier 1 sequon loss by cancer type; labels give tumours with a loss over tumours in the type. (B) The 15 most recurrent Tier 1 losses, with tumours carrying each substitution in all tumours and without hypermutated tumours, and the number of cancer types. Cancer type codes are listed in Supplementary Table S1.
ALT Horizontal bar chart of the share of tumours with a Tier 1 sequon loss in each of 33 cancer types, and a bar chart of the 15 most recurrent Tier 1 losses.

### Figure_S2
TITLE Substitution spectra behind the expected rates
LEGEND Share of substitutions in the 96 pyrimidine-centred trinucleotide channels for tumour somatic single-nucleotide variants (SNVs; 2,149,932 SNVs from one selected file per tumour, 33 cancer types, counts pooled over types for display, while the model uses per-type rates) and for germline missense variants (8,356,004 variants, primary variant set).
ALT Bar charts of the 96 trinucleotide substitution channels for tumour somatic SNVs and germline missense variants, the spectra used to build the expected rates.

### Figure_S3
TITLE Selection test across set-ups, calibration of the expected rate and planted effects
LEGEND (A) log2 rate ratio for Tier 1 N + S/T losses across set-ups (main, structure-matched strata, without hypermutated tumours, all primary-tumour MAF files, without rescued rows and position counts in the context-model protein set) in the position-count and context-aware models. (B) Calibration of the expected rate, as the log2 ratio of observed events at background positions of high against low expected rate (split at the median expected rate in each protein and residue); a calibrated model gives about 0. (C) Number of 20 planted data sets labelled as expected (first run) at planted rate ratios of 0.67 and 0.8 (expected Depleted), none (expected Equivalent) and 1.25 and 1.5 (expected Enriched), Tier 1 units; cells with 12 or more are marked. Rescued rows are mutations mapped to UniProt through another annotated transcript (Supplementary Methods, SM1). Bands, lines and bars as in Fig. 2.
ALT Forest plots of Tier 1 estimates across set-ups and of calibration ratios, and a tile chart of how many planted data sets each model labelled as expected.

### Figure_S4
TITLE Tier 1 sequon losses by cancer type, position-count model
LEGEND log2 rate ratio of Tier 1 sequon losses at glycosites against other positions in each of 33 TCGA cancer types, for N, S/T and N + S/T, in the position-count model. Benjamini-Hochberg correction within each column of 33 types; numbers on the right are glycosite events, and types with too few data are Not tested. Bands, lines and bars as in Fig. 2. Cancer type codes are listed in Supplementary Table S1.
ALT Forest plot of Tier 1 sequon loss estimates for each of 33 cancer types in the position-count model, in three columns for N, S/T and N plus S/T, coloured by label.

### Figure_S5
TITLE Tier 1 sequon losses by cancer type, context-aware model
LEGEND log2 rate ratio of Tier 1 sequon losses at glycosites against other positions in each of 33 TCGA cancer types, for N, S/T and N + S/T, in the context-aware model. Benjamini-Hochberg correction within each column of 33 types; numbers on the right are glycosite events, and types with too few data are Not tested. Bands, lines and bars as in Fig. 2.
ALT Forest plot of Tier 1 sequon loss estimates for each of 33 cancer types in the context-aware model, in three columns for N, S/T and N plus S/T, coloured by label.

### Figure_S6
TITLE dNdScv absolute dN/dS, planted effects and null control
LEGEND (A) Absolute dN/dS at Tier 1 glycosites and at background positions in the same genes, divided by the synonymous calibration of 1.059, with 95% CI (descriptive, no label). (B) Labels of 20 planted data sets at each planted rate ratio, Tier 1 units. (C) Null sets labelled Equivalent (of 20) for each unit; the rule requires at least 16 of 20 and no Changed call in any null set, and units below the rule are not allowed an Equivalent label.
ALT Three dNdScv panels show absolute dN/dS for glycosites and background, a bar chart of planted data set labels, and a bar chart of null sets labelled Equivalent against the rule of 16 of 20.

### Figure_S7
TITLE What each test can detect: planted effects in all pooled units
LEGEND Number of 20 planted data sets called as expected (first run, main set-up, pooled units) for the position-count and context-aware models, germline glycosites, germline Cys/His, sequon gains and dNdScv (Tier 1 only), at planted rate ratios as in Supplementary Fig. S3C; cells with 12 or more are marked.
ALT Tile chart of the number of planted data sets each test and pooled unit called as expected, at five planted rate ratios.

### Figure_S8
TITLE Labels of every pooled unit in every set-up
LEGEND Numbers are log2 rate ratios and colour gives the label after the control rules; set-ups beside the main run carry the label before the control rules where no final label exists. Panels show the position-count and context-aware models (tumour), germline glycosites, germline Cys (N row) and His (S/T row) by constrained class (metal-binding, active site, disulfide), sequon gains by region class (plausible, unknown, implausible) and dNdScv.
ALT Tile chart of the label and log2 rate ratio of every pooled unit in every set-up, grouped by test.

### Figure_S9
TITLE Signature shares of Tier 1 sequon losses by cancer type
LEGEND Shares of signature groups among Tier 1 N and S/T losses in each cancer type with enough losses (number of losses in brackets); numbers in tiles are percentages of 10% or more, and CpG C>T is counted beside the signature groups. Abbreviations as in Fig. 4.
ALT Heat map of mutational signature shares for Tier 1 N and S/T losses in each cancer type with enough losses.

### Figure_S10
TITLE Structural location of loss and matched control residues
LEGEND pLDDT and relative solvent exposure at the residue in AlphaFold models for Tier 1 sequon losses and matched control substitutions; only residues with pLDDT of 70 or more were modelled with FoldX.
ALT Plots of pLDDT and relative solvent exposure at loss residues and matched control residues.

### Figure_S11
TITLE Aromatic residue before the sequon
LEGEND Difference in median loss-minus-control score between losses with an aromatic residue before the sequon and the other losses, for FoldX ΔΔG, AlphaMissense and ESM1b (descriptive). Groups are F, Y or W at N-2 (primary, and with stricter matching), F only at N-2, F/Y/W at N-2 with T at +2, F/Y/W at N-3, F/Y/W at N-2 or N-3, and F at N-1 {Price, 2011}; n is the number of loss pairs in the group. An aromatic residue at N-2 raises sequon occupancy (enhanced aromatic sequon). Thick bars are 90% and thin bars 95% CIs.
ALT Forest plot of differences between losses with an aromatic residue before the sequon and the other losses, for three scores, all descriptive.

### Figure_S12
TITLE Protein amount by CPTAC cancer type and glycopeptides at lost sites
LEGEND (A) Difference in protein level, loss minus control carriers (log2), by CPTAC cancer type, with the number of loss pairs in brackets (descriptive). (B) Observed glycopeptide level at the lost site against the level expected by design, log2(1 − VAF) for a heterozygous carrier with equal allele expression, because a glycopeptide carrying the substitution cannot be identified against reference sequences; 20 rows at 20 sites, descriptive, never read as glycosylation loss or retention. BRCA, breast cancer; CCRCC, clear cell renal cell carcinoma; CRC, colorectal cancer; GBM, glioblastoma; HNSCC, head and neck squamous cell carcinoma; LSCC, lung squamous cell carcinoma; LUAD, lung adenocarcinoma; OV, ovarian cancer; UCEC, uterine corpus endometrial carcinoma; VAF, variant allele fraction.
ALT Forest plot of protein differences by CPTAC cancer type and a scatter plot of observed against expected glycopeptide levels at 20 lost sites.

### Figure_S13
TITLE Progression-free interval by cancer type
LEGEND Hazard ratio for progression-free interval, any loss against control only, in each TCGA cancer type with enough events (descriptive); numbers on the right are loss carriers (events). Thick bars are 90% and thin bars 95% CIs.
ALT Forest plot of hazard ratios for progression-free interval by cancer type, descriptive.

### Figure_S14
TITLE Clinical analysis: sensitivity analyses, leave-one-type-out and power
LEGEND (A) Hazard ratios of the primary and any-loss contrasts in sensitivity analyses (main, strict controls, without hypermutated tumours, without rescued rows, N losses only, S/T losses only, overall survival, with DLBC and KICH, unadjusted); labels by the pre-specified rule (Wald p below 0.05, otherwise TOST p below 0.05 for 0.80 to 1.25), descriptive. (B) Each cancer type left out in turn. (C) Power by simulation from a Weibull model fitted on the real cohort, as the share of simulated cohorts with p below 0.05 at each planted hazard ratio, for four contrasts and two simulation sets; the dashed line marks the 80% rule at a hazard ratio of 1.5. TOST, two one-sided tests.
ALT Forest plots of hazard ratios in sensitivity analyses and leave-one-type-out runs, and power curves from simulation.

### Figure_S15
TITLE Calibration of the gain model and FoldX of sequon-creating substitutions
LEGEND (A) Calibration of the expected rate for gains (as in Supplementary Fig. S3B) for new N, new S/T and both; the pooled unit passed, and new S/T alone gave 0.15 because the expected rate is less exact for S/T-creating changes. (B) Cumulative FoldX ΔΔG of plausible (2,494) and unknown (1,944) gains, with Tier 1 losses (1,866) for reference (descriptive).
ALT Forest plot of calibration ratios for the gain model and cumulative ΔΔG curves for gains beside Tier 1 losses.

### Figure_S16
TITLE Top background hotspots
LEGEND The 20 non-sequon N, S and T positions with the most events in TCGA tumours, coloured by whether the gene is a known driver gene {Bailey, 2018}.
ALT Bar chart of the 20 background N, S and T positions with the most events, coloured by known driver gene.

### Figure_S17
TITLE All headline estimates scaled to their equivalence margins
LEGEND Each headline estimate and its 90% CI divided by its pre-specified equivalence margin, grouped into selection, structure and predicted effect, and protein amount and survival; the text gives each estimate, its CI and its margin. The grey band is the region inside the margin, and labels are as pre-specified.
ALT Forest plot of every headline estimate divided by its margin. Selection, structural and protein estimates lie inside the band, the primary progression-free interval contrast crosses its lower margin, and the germline Cys and His positive control lies far below it.

### Table_1
TITLE Data sources, cohorts, evidence tiers and set-ups
HEADER Resource | Release or version | What was used | Size | Used for
ROW TCGA somatic mutations | GDC masked MAF files (GRCh38, GENCODE v36 transcripts) | Single-nucleotide variants of one selected primary-tumour file per tumour; missense mutations mapped to UniProt positions | 10,190 tumours in 33 cancer types; 1,460,792 missense mutations; 2,149,932 SNVs | Catalogue, selection tests, signatures, hotspots, ActiveDriver
ROW Human proteome and glycosites | UniProt 2026_03 (reviewed human entries); N-GlycositeAtlas | Sequons N-X-S/T (X not P) and their evidence tiers | Tier 1, 13,247; Tier 2, 8,087; Tier 3, 38,218 sequons | Every analysis
ROW Germline missense variants | UniProt 2026_03 human variation file (gnomAD, TOPMed, ExAC, 1000 Genomes, ESP) | Missense variants at sequon and comparison positions | 8,356,004 variants | Germline comparison; positive control at constrained Cys/His
ROW Constrained residues | UniProt 2026_03 features | Metal-binding, active-site and disulfide Cys and His | See Supplementary Table S14 | Positive control at constrained Cys/His
ROW Protein structures | AlphaFold Protein Structure Database v6, human proteome UP000005640 | pLDDT, relative solvent accessibility and models | One model or fragment set per protein | Structure-matched strata, FoldX, sequon creation
ROW Predicted effects | FoldX 5; AlphaMissense amino-acid substitution scores; ESM1b log-likelihood ratios (all human isoforms) | ΔΔG, pathogenicity and damage scores of losses, matched controls and gains | 1,443 loss-control pairs | Predicted structural and functional impact
ROW Clinical data | TCGA Clinical Data Resource | Progression-free interval, overall survival, age and sex | 9,433 patients; 6,101 patients in 29 cancer types in the Cox model | Clinical analysis
ROW CPTAC proteomes | Glycoprotein-Notebook global proteome and glycoproteome tables; GDC CPTAC masked MAF files | Protein and RNA levels of loss, control missense and truncating carriers | 1,017 tumours with a MAF in 10 cancer types | Protein amount in carriers
ROW Driver genes | Bailey et al. 2018, Table S1 | Known cancer driver genes | 299 genes | Set-up without driver genes; hotspot and ActiveDriver annotation
ROW Published variants | Fan et al. 2018 | Recurrent loss and gain variants | 13 loss and 39 gain variants | Cross-check with Fan et al. 2018
ROW Prediction and modelling of new sequons | NetNGlyc 1.0; Re-Glyco 0.2.6 with GlycoShape; SURFY | Sequon potential, glycan placement and surface annotation | 614 shortlisted substitutions; top 10 modelled | Sequon creation and shortlist
HEADER Set-up | Definition | Used in
ROW Main | One selected primary-tumour MAF file per tumour, all 33 cancer types, rescued rows included | All analyses
ROW Without hypermutated tumours | Tumours with more than 1,000 missense rows in the selected file left out | Selection tests, gains, hotspots, protein amount, clinical analysis
ROW All primary-tumour MAF files | Union of the primary-tumour files of a tumour, each variant counted once per tumour | Selection tests, gains
ROW Without rescued rows | Rows mapped to UniProt through another annotated transcript left out | Selection tests, gains, clinical analysis
ROW Structure-matched strata | Strata also matched by AlphaFold pLDDT (70 or more, below 70) and relative exposure (0.25 or more, below 0.25); proteins with a usable model only | Selection tests, germline comparison, gains
ROW Position counts, context-model protein set | Position-count model restricted to the proteins of the context-aware model | Selection tests, germline comparison, Cys/His control, gains
ROW Germline variants in two or more resources | Variants reported by at least two germline resources | Germline comparison, Cys/His control
ROW Window mapping (dNdScv) | Glycosite positions mapped onto the dNdScv reference coding sequences by sequence window, beside strict mapping [confirm: Supplementary Methods SM6] | dNdScv
ROW Second seed (dNdScv) | Main run repeated with a second random seed for the bootstrap | dNdScv
ROW No mutation caps (dNdScv) | Default dNdScv caps on mutations per gene and per sample removed | dNdScv
ROW Without driver genes (dNdScv) | The 299 driver genes of Bailey et al. 2018 left out | dNdScv
ROW Leave-one-out (dNdScv) | Each of the five largest TCGA projects (UCEC, SKCM, COAD, LUAD, STAD) left out in turn | dNdScv
FOOT Cancer type codes are listed in Supplementary Table S1. Tier 1, experimental UniProt evidence or N-GlycositeAtlas site; Tier 2, other UniProt N-glycosylation annotation; Tier 3, motif only.
FOOT CPTAC, Clinical Proteomic Tumor Analysis Consortium; ESP, NHLBI Exome Sequencing Project; ExAC, Exome Aggregation Consortium; GDC, Genomic Data Commons; gnomAD, Genome Aggregation Database; MAF, mutation annotation format; pLDDT, predicted local distance difference test; SNV, single-nucleotide variant; SURFY, in silico human surfaceome; TCGA, The Cancer Genome Atlas; TOPMed, Trans-Omics for Precision Medicine; ΔΔG, change in folding free energy.

### Table_2
TITLE Headline estimates with their margins, labels, controls and analysis status
HEADER Analysis | Unit and cohort | Estimate (90% CI) | Margin | Label | Controls | Analysis status | Figure
T2ROW q2 | Position-count selection test | Tier 1 N + S/T losses, TCGA tumours, 33 types | ±0.14 log2 | Planted-effect and null controls met; germline biology control failed | Primary (pre-specified) | Fig. 2A
T2ROW ctx | Context-aware selection test | Tier 1 N + S/T losses, TCGA tumours, 33 types | ±0.14 log2 | Calibration passed; planted-effect and null controls met; germline biology control failed | Added after the primary analysis | Fig. 2A
T2ROW dnds | dNdScv expected rates | Tier 1 N + S/T losses, TCGA tumours, 33 types | ±0.14 log2 | Null control met (16 of 20); planted-effect, synonymous calibration and driver positive controls passed | Post hoc, not blind | Fig. 2B, C
T2ROW germ | Germline comparison | Tier 1 N + S/T, germline missense variants | ±0.14 log2 | Positive control at constrained Cys/His passed | Added after the primary analysis | Fig. 3A
T2ROW cys | Positive control at constrained Cys/His | Metal-binding Cys + His, germline missense variants | ±0.14 log2 | This row is the positive control of the method | Added after the primary analysis | Fig. 3B
T2ROW sig | Mutational signatures | CpG C>T share, Tier 1 S/T losses minus Tier 1 N losses (percentage points, 95% CI) | None | Descriptive | Added after the primary analysis | Fig. 4C
T2ROW fan | Cross-check with Fan et al. 2018 | Position-3 excess, log2 observed over expected S/T-to-N loss ratio, Tier 1 | None | Descriptive | Post hoc, not blind | Fig. 4E
T2ROW foldx | FoldX ΔΔG | Tier 1 N + S/T losses minus matched controls (kcal/mol), 1,443 pairs | ±0.5 kcal/mol | Second control draw and stricter matching give the same label | Primary (pre-specified) | Fig. 5A, B
T2ROW am | AlphaMissense | Tier 1 N + S/T losses minus matched controls, 1,443 pairs | ±0.082 | Second control draw gives the same label | Added after the primary analysis | Fig. 5B
T2ROW esm | ESM1b | Tier 1 N + S/T losses minus matched controls, 1,443 pairs | ±0.71 | Second control draw gives the same label | Added after the primary analysis | Fig. 5B
T2ROW cptac | Protein amount in CPTAC carriers | Loss minus control missense carriers, within protein, RNA-adjusted (log2), 81 proteins | ±0.14 log2 | Null and planted controls passed; truncating positive control failed | Added after the primary analysis | Fig. 6C
T2ROW pfi1 | Progression-free interval | Loss only against control only, Cox model stratified by 29 types | HR 0.80 to 1.25 | Power by simulation (Supplementary Fig. S14C) | Primary (pre-specified) | Fig. 6F
T2ROW pfi2 | Progression-free interval | Any loss against control only, Cox model stratified by 29 types | HR 0.80 to 1.25 | Holm correction over three secondary contrasts | Primary (pre-specified), secondary contrast | Fig. 6F
T2ROW gain | Sequon creation | New N + new S/T in plausible regions, TCGA tumours, 33 types | ±0.14 log2 | Calibration passed; planted-effect and null controls met | Added after the primary analysis | Fig. 7A
T2ROW hot | Single-position hotspots | Tier 1 glycosites against background N, S and T positions, odds ratio of calls (95% CI) | None | Planted and null controls passed; calibration failed (ratio 3.11) | Added after the primary analysis | Fig. 8B
T2ROW ad | ActiveDriver | Proteins with an enriched glycosite region, real Tier 1 sites against null sets | 5% null rule | Positive control called in 4 of 5 proteins | Added after the primary analysis | Fig. 8E
FOOT Estimates are given with 90% confidence intervals, the basis of the equivalence test (two one-sided tests at the 5% level), except where 95% is stated. Labels follow the pre-specified rules (Supplementary Methods, SM20). Rate ratios compare glycosite positions with other N, S and T positions of the same proteins; for constrained residues, other Cys and His positions; for gains, other substitutions.
FOOT CI, confidence interval; CPTAC, Clinical Proteomic Tumor Analysis Consortium; HR, hazard ratio; TCGA, The Cancer Genome Atlas; ΔΔG, change in folding free energy.

### Supplementary_tables
STAB S1 | Cancer types and tumour counts | TCGA cancer type codes, full names, tumours and missense mutations
STAB S2 | Sequon losses and gains by cancer type | Mutations and tumours in each class, all tumours and without hypermutated tumours
STAB S3 | Position of the lost residue | Sequon losses by tier and position in the sequon, by cancer type
STAB S4 | Most recurrent sequon-changing substitutions | Tumours and cancer types for the most recurrent losses and gains
STAB S5 | Position-count selection test | Every unit, set-up and cancer type
STAB S6 | Context-aware selection test | Every unit, set-up and cancer type, with calibration status
STAB S7 | Calibration of the expected rates | High against low expected background, context model and uniform rate
STAB S8 | dNdScv selection test by set-up | Every tier and position in every dNdScv set-up
STAB S9 | dNdScv leave-one-out | Tier 1 units with each of the five largest TCGA projects left out
STAB S10 | dNdScv absolute dN/dS | Observed and expected events, raw and calibrated dN/dS
STAB S11 | What each test can detect | Planted-effect and null results for every unit of every test
STAB S12 | Substitution spectra | Tumour SNVs by cancer type and channel, and germline variants by channel
STAB S13 | Germline comparison | Germline missense variants at glycosite positions, every unit and set-up
STAB S14 | Positive control at constrained Cys/His | Germline variants at metal-binding, active-site and disulfide Cys and His
STAB S15 | Signature shares, pooled | Signature group and CpG C>T shares of each mutation set
STAB S16 | Signature shares by cancer type | Signature group and CpG C>T shares of each mutation set in each cancer type
STAB S17 | Differences in signature shares | Differences between mutation sets with 95% CI
STAB S18 | Counts beside Fan et al. 2018 | Counts reported by Fan et al. 2018 and the closest counts of this study
STAB S19 | Variants of Fan et al. 2018 in TCGA | Each published loss and gain variant and whether it is seen in TCGA tumours
STAB S20 | Position-3 excess | Observed and expected S/T-to-N loss ratios by tier
STAB S21 | FoldX ΔΔG per substitution | Every modelled loss, matched control and method control
STAB S22 | Losses against matched controls, three scores | FoldX, AlphaMissense and ESM1b differences in every set-up
STAB S23 | Proteome-wide scores at sequon N | Every possible substitution at sequon N against background N, with null ranges
STAB S24 | Aromatic residue before the sequon | Group-minus-rest differences for every aromatic split and score
STAB S25 | CPTAC loss carrier-protein pairs | Every loss carrier-protein pair, protein values and reasons for leaving pairs out
STAB S26 | Protein amount estimates | Primary estimate, controls and descriptive patterns
STAB S27 | Glycopeptides at lost sites | Observed and expected glycopeptide levels at lost Tier 1 sites
STAB S28 | Progression-free interval models | Every Cox model and contrast, with 90% CIs, equivalence tests and labels
STAB S29 | Progression-free interval by cancer type | Any loss against control only in each cancer type
STAB S30 | Progression-free interval, leave one type out | Primary and secondary contrasts with each cancer type left out
STAB S31 | Sequon creation test | Every gain unit and set-up
STAB S32 | Shortlist of sequon-creating substitutions | The 614 shortlisted substitutions with structure, FoldX and NetNGlyc values
STAB S33 | Top 10 with Re-Glyco | Glycan modelling of the top 10 shortlisted substitutions
STAB S34 | FoldX of sequon-creating substitutions | Summary of ΔΔG for plausible and unknown gains
STAB S35 | Single-position test at Tier 1 glycosites | Every Tier 1 position with at least one event
STAB S36 | Background hotspots | Non-sequon N, S and T positions called hotspots
STAB S37 | Hotspot summary and controls | Odds ratios of Tier 1 against background, and the known-hotspot positive control
STAB S38 | ActiveDriver test of glycosite regions | Every protein tested
STAB S39 | ActiveDriver controls and agreement | Positive and null controls, and agreement with single-position hotspots
STAB S40 | Software and versions | Programs, data releases and package versions
