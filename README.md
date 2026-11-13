# SqnPAN analysis code

Code, supplementary tables and glycan models for the article

> Balasundar A. **Somatic loss of N-glycosylation sequons occurs at the rate expected from mutational processes across 33 tumour types.** (submitted, 2026)

The code reproduces every result, table and figure of the article from public files. Somatic missense mutations of 10,190 TCGA tumours (33 cancer types) are mapped to N-glycosylation sequons (N-X-S/T) with experimental evidence, and sequon-destroying substitutions are compared with substitutions at other N, S and T positions of the same proteins by a position-count test, a trinucleotide-context model and dNdScv, each against an equivalence margin of ±0.14 log2 (about 10%), with null and planted-effect controls. Added analyses cover germline variation, mutational signatures, FoldX, AlphaMissense and ESM1b scores, CPTAC protein levels, sequon gains, hotspots, clinical links and a cross-check with Fan et al. (2018).

Author: Anishkiran Balasundar, Independent Researcher, Kanyakumari district, Tamil Nadu, India. ORCID [0009-0007-7056-5356](https://orcid.org/0009-0007-7056-5356)

## Contents

```
Codework_SqnPAN/
  README.md                 this file
  LICENSE                   MIT licence (code)
  CITATION.cff              how to cite this archive
  .zenodo.json              Zenodo metadata
  requirements.txt          Python packages (versions used)
  install_R_packages.R      R packages (versions used are listed below)
  scripts/                  analysis scripts SqnPAN_01 to SqnPAN_26 (run in number order)
  writing/                  SqnPAN_26_legend_text.md, text read by SqnPAN_26
  data/README_data.md       every input file: where it goes, where it comes from, its md5
  data/raw/                 GDC manifests and case maps, NetNGlyc 1.0 result pages
  results/reglyco_gain/     mutant models and Re-Glyco glycan models of the top 10 sequon gains (Figure 7F)
  results/tables/           P2_19_reglyco_filled.tsv (Re-Glyco results, filled by hand from the web runs)
  supplementary/            SqnPAN_SupplementaryTables_S01-S40.xlsx (Supplementary Tables S1 to S40)
```

## Software

Analyses were run on Windows with the versions below (also in Supplementary Table S40).

| Software | Version |
|---|---|
| Python | 3.14.0; numpy 2.3.5, pandas 3.1.0, openpyxl 3.1.5, python-docx 1.2.0, psutil 7.2.2 |
| R | 4.6.0; survival 3.8.6, MutationalPatterns 3.22.0, ActiveDriver 1.0.0, dndscv 0.0.1.0 (GitHub im3sanger/dndscv), ggplot2 4.0.3, patchwork 1.3.2, ggrepel 0.9.8, scales 1.4.0, ragg 1.5.2, systemfonts 1.3.2, readxl 1.5.0, BiocManager 1.30.27 |
| FoldX | 5.1, Windows build `foldx_1_20270131.exe` (academic licence from https://foldxsuite.crg.eu, not included) |
| Web tools | NetNGlyc 1.0 (DTU Health Tech), Re-Glyco 0.2.6 (GlycoShape) |

```
python -m pip install -r requirements.txt
Rscript install_R_packages.R
```

Figures use the Times New Roman font.

## Folder set-up

All paths are relative to this folder. Python scripts find it from their own location. R scripts use the folder above `scripts/` when run with `Rscript`; when run with `source()` (RStudio) they use the working directory, so set the working directory to this folder first, or set the environment variable `SQNPAN_ROOT` to its full path. Place the input files as listed in `data/README_data.md`. The scripts create `data/processed/`, `results/logs/`, `results/tables/` and `results/figures/`.

## Run order

Run every command from this folder. Each script writes a log per section to `results/logs/`, checks the md5 of its inputs and stops when a check fails. A later step refuses to run until the gate file of the step before it exists.

| Step | Command | What it does |
|---|---|---|
| 1 | `python scripts/SqnPAN_01_verify_maf_download.py data/raw results/logs/SqnPAN_01_log.txt` | TCGA masked MAF download against the manifest (presence, size, md5, header) |
| 2 | `python scripts/SqnPAN_02_verify_maf_chunks.py chunk ...` / `summary ...` | Same checks in restartable chunks (usage in the header) |
| 3 | `python scripts/SqnPAN_03_check_maf_uniprot_mapping.py chunk ...` / `summary ...` | MAF protein columns against UniProt 2026_03 sequences |
| 4 | `python scripts/SqnPAN_04_check_alphafold_archive.py all` | AlphaFold v6 human archive: index, model reading, UniProt comparison, use table (helpers in `SqnPAN_04_alphafold_helpers.py`) |
| 5 | `python scripts/SqnPAN_05_structure_check.py all` | Input guards, table shapes, glycosite tiers, one MAF file per case |
| 6 | `python scripts/SqnPAN_06_build_tables.py all` | Control C1, site tables, mutation table with loss and gain classes |
| 7 | `python scripts/SqnPAN_07_structure_features.py 1`, then `2`, then `3` | Per-residue pLDDT and relative solvent exposure |
| 8 | `python scripts/SqnPAN_08_selection_controls.py 1`, then `2`, then `3` | Q2 controls C2 (null) and C3 (planted effects) |
| 9 | `python scripts/SqnPAN_09_selection_test.py all` | Q2 primary result (position-count model) and final labels |
| 10 | `python scripts/SqnPAN_10_foldx_setup.py 1`, then `2` | FoldX format test and timing |
| 11 | `python scripts/SqnPAN_11_foldx_runs.py all` | FoldX ΔΔG of losses against matched controls, and of gains |
| 12 | `Rscript scripts/SqnPAN_12_clinical.R` | Clinical links (progression-free interval), control C5, driver genes |
| 13 | `python scripts/SqnPAN_13_fan_crosscheck.py all` | Cross-check with Fan et al. (2018) |
| 14 | `python scripts/SqnPAN_14_context_q2.py all` | A1: Q2 with trinucleotide-context expected rates |
| 15 | `python scripts/SqnPAN_15_germline_q2.py all` | A2: germline positive control |
| 16 | `Rscript scripts/SqnPAN_16_signatures.R` | A3: mutational signatures of sequon losses |
| 17 | `python scripts/SqnPAN_17_constrained_check.py all` | Germline test at metal-binding, active-site and disulfide Cys and His |
| 18 | `python scripts/SqnPAN_18_variant_effect_aromatic.py all` | B1: AlphaMissense and ESM1b; C1: aromatic sequon split |
| 19 | `python scripts/SqnPAN_19_gain_side.py all`, then `6` | B3: sequon gains, shortlist, NetNGlyc and Re-Glyco inputs (manual steps below) |
| 20 | `python scripts/SqnPAN_20_hotspots.py all` | B4: single-site hotspots |
| 21 | `Rscript scripts/SqnPAN_21_activedriver.R` | B4: ActiveDriver |
| 22 | `python scripts/SqnPAN_22_cptac_protein.py 1`, then `2` | B2: protein level in CPTAC carriers (confirmation step below) |
| 23 | `python scripts/SqnPAN_23_dndscv_prep.py`, then `Rscript scripts/SqnPAN_23_dndscv.R` | dNdScv at Tier 1 glycosite residues |
| 24 | `Rscript scripts/SqnPAN_24_figures.R draft` | Theme, palette and panel code; draft renditions |
| 25 | `Rscript scripts/SqnPAN_25_figures_final.R final` | Figures 1 to 8 and Supplementary Figures S1 to S17 (PDF, 600 dpi TIFF, PNG) |
| 26 | `python scripts/SqnPAN_26_tables_legends_workbook.py` | Supplementary Tables S1 to S40 workbook, Tables SM2 and SM3, legends and alt text |

Files ending in `_helpers.py` or `_helpers.R` hold shared functions and are not run on their own. Scripts 15 and 17 need step 14; 19, 20 and 22 need 14; 21 needs 20; 23 needs 8, 14 and 20; 25 reads the panel code of 24; 26 reads the figure index of 25.

## Manual steps

1. **FoldX** (steps 5, 10, 11, 19). Put the FoldX 5.1 Windows build in `foldx5_Windows/` in this folder. FoldX builds expire, and the scripts check the md5 of the build that was used (`FOLDX_EXE` and `FOLDX_MD5` in `SqnPAN_10_foldx_helpers.py`); with another build, change both values and note that ΔΔG values can differ.
2. **NetNGlyc 1.0** (step 19). Section 3 writes `data/processed/P2_19_netnglyc_input_<k>.fasta`. Submit each file to NetNGlyc 1.0 with default settings and save the result pages in `data/raw/netnglyc/` (the pages used are included). Then run `python scripts/SqnPAN_19_gain_side.py 6`.
3. **Re-Glyco** (step 19, Figure 7F). Section 6 writes the mutant models of the top 10 substitutions to `results/reglyco_gain/<gene_mutation>/` and `results/tables/P2_19_reglyco_input.tsv`. Each model was run in Re-Glyco 0.2.6 (GlycoShape) with one complex N-glycan (GlyTouCan G54258NG) at the new sequon; the outputs are in `results/reglyco_gain/` and the clash results in `results/tables/P2_19_reglyco_filled.tsv`.
4. **CPTAC counts** (step 22). Section 1 prints the counts and the analysis route. Compare them with the COUNTS block in the header of `SqnPAN_22_cptac_protein.py`, then create `data/processed/P2_22_counts_confirmed.txt` with one line that starts with `CONFIRMED` and contains `route=within` (the route printed by section 1). Section 2 refuses to run without this file.

## Notes for reading the code

- **Pre-specification.** The header of each analysis script holds its dated pre-specification (question, rules, controls, seeds, labels), written before the result was computed. Changes made after a first run are written as dated amendments in the same header and are called post hoc. "Analysis plan" refers to the study plan; its sections, rules (for example rule 2A) and Amendments 1 to 4 are described in the Supplementary Methods.
- **Labels used in the code.** Q1 classification of losses and gains; Q2 rate of sequon-destroying substitutions against background substitutions; Q3 structural effect (FoldX). T1, T2, T3 are Tier 1 (experimental evidence), Tier 2 (other UniProt annotation) and Tier 3 (motif only). Added analyses (Amendment 4): A1 context-aware Q2, A2 germline positive control, A3 mutational signatures, B1 AlphaMissense and ESM1b, B2 CPTAC protein level, B3 gain side and shortlist, B4 hotspots and ActiveDriver, C1 aromatic sequon split. C-xx are the controls of each script. "Equivalent" means that the 90% confidence interval lies inside ±0.14 log2.
- **File names.** Scripts carry the SqnPAN prefix and a serial number. The intermediate files they write keep the prefix `P2_<serial>_` of the step that made them (for example `P2_06_mutations.tsv.gz` from SqnPAN_06), because later scripts and the md5 checks refer to these names.
- **Seeds.** Every random step uses a fixed seed from a reserved block named in the script header; repeated runs give identical tables.
- **md5 checks.** Scripts check the md5 of their inputs, and some check the md5 of other scripts (SqnPAN_09 checks SqnPAN_08; SqnPAN_11 checks both SqnPAN_10 files; SqnPAN_12 checks its helpers; SqnPAN_13 checks SqnPAN_06_helpers.py). Editing those scripts, or a change of line endings, makes the check stop the run; `.gitattributes` keeps the files byte for byte in git.
- **Test options.** Options such as `limit=N`, `test=1` and the environment variables `P2_LOGDIR`, `P2_TABLEDIR`, `P2_PROCDIR`, `P2_WORKERS`, `P2_NBOOT` are for quick test runs only and were not used for the reported results.
- **Run time.** The FoldX runs (step 11) and the resampling steps (8, 9, 14, 15) take many hours on a laptop; sections that stop can be restarted and continue from saved parts.

## Licence

The code is released under the MIT licence (`LICENSE`). The supplementary workbook and the model files are released under CC BY 4.0. Values derived from third-party resources keep the terms of their providers; AlphaMissense (CC BY-NC-SA 4.0) and ESM1b (CC BY-NC 4.0) scores are for non-commercial research only.

## Citation

Please cite the article and this archive (see `CITATION.cff`).
