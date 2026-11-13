# Input files

All inputs are public. Put them at the paths below (relative to the main folder). The md5 values are the ones the scripts check; a different file makes the check print DIFFERENT or stop the run. Exact sources, releases and download dates are in Table SM1 of the Supplementary Methods, and software and data releases in Supplementary Table S40.

Files marked **included** are in this archive. The others are too large or are third-party files that must be downloaded from their providers.

## Mutations (GDC, open access)

| Path | Content | md5 |
|---|---|---|
| `data/raw/TCGA_masked_maf_manifest_2026-10-21.txt` | **included**. GDC manifest of the 10,640 TCGA masked somatic MAF files (33 projects) | |
| `data/raw/maf_tcga/<id>/<filename>` | the MAF files, downloaded with the GDC Data Transfer Tool 2.3 (`gdc-client download -m <manifest> -d data/raw/maf_tcga`); size and md5 of each file are in the manifest | |
| `data/raw/GDC_TCGA_masked_maf_case_project_2026-10-22.tsv` | **included**. File, case and project table of the TCGA MAF files | b6bd436e630584d4738d9bdacf4d5b7a |
| `data/raw/CPTAC_masked_maf_manifest.txt` | **included**. GDC manifest of the CPTAC masked somatic MAF files | b10990250fd4d89d77977ec4a1dba003 |
| `data/raw/maf_cptac/<id>/<filename>` | the CPTAC MAF files, downloaded the same way | |
| `data/raw/GDC_CPTAC_masked_maf_case_map_2026-11-09.tsv` | **included**. GDC file to case map of the CPTAC MAF files | cf1cd301ad6ed2746fb91c531479801e |

## Proteins, glycosites and structures

| Path | Content | md5 |
|---|---|---|
| `data/raw/Uniprot/stream.gz` | UniProtKB/Swiss-Prot human reviewed, release 2026_03, TSV with Entry, Entry Name, Gene Names (primary), Length, Sequence, Glycosylation, Ensembl, AlphaFoldDB | af9af22f1e6eea80e96751c93c9da2b3 |
| `data/raw/Uniprot/uniprot_2026_03_features_stream.gz` | same entries, TSV with Entry, Length, Sequence, Topological domain, Transmembrane, Intramembrane, Signal peptide, Subcellular location [CC], Domain [FT], Disulfide bond, Active site, Binding site | bf26668a884779609e7c1fe46bb6aec4 |
| `data/raw/Uniprot_variation/homo_sapiens_variation.txt.gz` | UniProt human variation file, release 2026_03 | ddfc80d965ebafb900d64fe8c51d7cec |
| `data/raw/Uniprot_variation/humsavar.txt` | UniProt humsavar, release 2026_03 | |
| `data/raw/nglycositeatlas/HumanAll.xlsx` | N-GlycositeAtlas, human table | 23d97d828400e59f1ca0837b5d59cf8a |
| `data/raw/alphafold/UP000005640_9606_HUMAN_v6.tar` | AlphaFold Protein Structure Database, human proteome v6 (5,177,506,304 bytes; never unpacked) | 5bc6bf3812832ceb6fe8f7c9d64fd6be |
| `data/raw/gencode/gencode.v36.pc_transcripts.fa.gz` | GENCODE v36 protein-coding transcripts (the GDC MAF annotation) | 61eaea20b1275e2e19def3e5882ed7ae |
| `data/raw/gencode.v42.pc_translations.fa.gz` | GENCODE v42 protein-coding translations (CPTAC Glycoprotein-Notebook reference) | 2b83b9c1751f60a7cc3ad499973a80af |
| `data/raw/table_S3_surfaceome.xlsx` | SURFY in silico human surfaceome, Table S3 (Bausch-Fluck et al., 2018) | bf95ed25ffc933fbe0c0383ab19b8652 |

## Variant effect scores

| Path | Content | md5 |
|---|---|---|
| `data/raw/alphamissense/AlphaMissense_aa_substitutions.tsv.gz` | AlphaMissense amino-acid substitution scores (Zenodo, version 3) | |
| `data/raw/esm1b/ALL_hum_isoforms_ESM1b_LLR.zip` | ESM1b log-likelihood ratio scores, all human isoforms | 12c4a5d2e8d2466196748da038cb0617 |

## Clinical and driver tables

| Path | Content | md5 |
|---|---|---|
| `data/raw/TCGA-CDR-SupplementalTableS1.xlsx` | TCGA Clinical Data Resource, Table S1 (Liu et al., 2018) | a4591b2dcee39591f59e5e25a6ce75fa |
| `data/raw/Bailey2018_TableS1_driver_genes.xlsx` | Cancer driver genes, Table S1 (Bailey et al., 2018) | 144a9ccbbf038236decc1226df24620b |

## CPTAC protein and glycopeptide tables (Glycoprotein-Notebook, Zhang et al., 2025)

One zip per cohort in `data/raw/`, read by SqnPAN_22.

| File | Cohort | md5 |
|---|---|---|
| `BRC_P_PDC000121.zip` | BRCA | 9266f7f1dbaab7272cf417577a9d9c00 |
| `CCRCC_P_PDC000128.zip` | CCRCC | 0e6efc6795a575c877629cdfa4da6c97 |
| `CRC_P_PDC000117.zip` | CRC | 42bbe35ac99c3c03836b077954b7609b |
| `GBM_P_PDC000205.zip` | GBM | 9a91371d09d91859961e8e6cf6884969 |
| `HNSCC_P_PDC000222.zip` | HNSCC | 5ffc664b61c35dd4974cc6ce3f70bf37 |
| `LSCC_P_PDC000232.zip` | LSCC | 13a137f91ac9c7e58749910da1a6fd25 |
| `LUAD_P_PDC000149.zip` | LUAD | 8d29ffce58d8dacd8307df31c931b61c |
| `OVC_G_PDC000250.zip` | OV | 3ca47a71c659bb24b1bec5a2ac71fa82 |
| `PDAC_PandG_271_272.zip` | PDAC | a4e25b62509f40f5d3334070b64e4c7d |
| `UCEC_P_PDC000126.zip` | UCEC | 39be02a0ab90544cfcdbd7088ac93f7a |

`CCRCC_G_PDC000471.zip`, `UCEC_G_PDC000445.zip` and `UCEC_P_PDC000441.zip` were also downloaded; their protein and RNA files are byte-identical to the files used (checked in SqnPAN_22).

## Web-tool outputs

| Path | Content |
|---|---|
| `data/raw/netnglyc/` | **included**. NetNGlyc 1.0 result pages for the shortlist FASTA files of SqnPAN_19 section 3 |
| `results/reglyco_gain/` | **included**. Mutant models (FoldX) and Re-Glyco 0.2.6 outputs of the top 10 sequon gains |
| `results/tables/P2_19_reglyco_filled.tsv` | **included**. Re-Glyco clash results of the top 10, entered by hand |

## Program

| Path | Content | md5 |
|---|---|---|
| `foldx5_Windows/foldx_1_20270131.exe` | FoldX 5.1 Windows build (academic licence, https://foldxsuite.crg.eu; not included) | c66ffd383d4d1d274e940bf99979bf00 |

## Values typed in the code

The Fan et al. (2018) variant tables used by SqnPAN_13 are typed into the script as `FAN_` constants, with their source in the script header.
