#!/usr/bin/env python3
# SqnPAN_26_tables_legends_workbook.py - supplementary tables workbook and the table and legend documents
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Builds, from the finished result tables (read only) and the text file writing/SqnPAN_26_legend_text.md:
#   writing/P2_Tables_1_2.docx                 Table 1 (two parts) and Table 2 (landscape), Times New Roman
#   writing/P2_Legends_AltText.docx            legends of 8 main and 17 supplementary figures, alt text, table titles
#   results/tables_final/P2_Supplementary_Tables.xlsx   README + Supplementary Tables S1 to S40
#   results/tables_final/P2_Table_2.tsv        Table 2 numbers as computed (for checking)
#   results/tables_final/P2_26_scan.tsv        working-term scan of every sheet name, header and text cell
#   results/tables_final/P2_26_log.txt         log (also printed)
# No random numbers are used (no seed). Estimates in Table 2 are read from the result tables, never typed.
# Run (from the repository folder):  python scripts/SqnPAN_26_tables_legends_workbook.py
# Needs: pandas, openpyxl, python-docx  (pip install pandas openpyxl python-docx)
import os, sys, re, math, shutil, importlib, importlib.util, subprocess, traceback, datetime
from pathlib import Path

missing = [m for m in ("pandas", "openpyxl", "docx") if importlib.util.find_spec(m) is None]
if missing:
    print("Missing Python packages: " + ", ".join(missing))
    print("Install with:  python -m pip install pandas openpyxl python-docx")
    sys.exit(1)
import pandas as pd
import numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter
import docx
from docx.shared import Pt, Mm
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
OUTT = ROOT / "results" / "tables_final"; OUTT.mkdir(parents=True, exist_ok=True)
OUTW = ROOT / "writing"
MD = OUTW / "SqnPAN_26_legend_text.md"
FIGIDX = ROOT / "results" / "figures" / "P2_25_final" / "P2_25_figure_index.tsv"
RSCRIPT = Path(os.environ.get("RSCRIPT") or shutil.which("Rscript") or r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe")
LOGF = OUTT / "P2_26_log.txt"; LOGF.write_text("", encoding="utf-8")
ERRORS = []
MINUS = "\u2212"; PM = "\u00b1"; DDG = "\u0394\u0394G"

def lg(*a):
    m = datetime.datetime.now().strftime("%H:%M:%S") + "  " + " ".join(str(x) for x in a)
    print(m)
    with open(LOGF, "a", encoding="utf-8") as f: f.write(m + "\n")

def rt(name):
    return pd.read_csv(TAB / (name + ".tsv"), sep="\t", low_memory=False, keep_default_na=False, na_values=["", "NA", "NaN", "nan"])

lg("SqnPAN_26 start; root", ROOT)

# ---------------------------------------------------------------- 1. legend text
def read_md(path):
    blocks, cur = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("### "):
            cur = line[4:].strip(); blocks[cur] = []
        elif cur is not None and line.strip():
            blocks[cur].append(line.rstrip())
    return blocks

B = read_md(MD)
def field(block, key):
    for l in B[block]:
        if l.startswith(key + " "): return l[len(key) + 1:]
    raise KeyError(block + " " + key)

FIGS = ["Figure_%d" % i for i in range(1, 9)] + ["Figure_S%d" % i for i in range(1, 18)]
for f in FIGS:
    for k in ("TITLE", "LEGEND", "ALT"): field(f, k)
lg("text file read:", len(FIGS), "figures, Table_1, Table_2, supplementary table titles")
try:
    idx = pd.read_csv(FIGIDX, sep="\t")
    mism = [r.file for r in idx.itertuples() if r.file in B and field(r.file, "TITLE") != r.title]
    lg("figure titles that differ from the figure index:", ", ".join(mism) if mism else "none")
except Exception as e:
    lg("figure index not read:", e)

STAB = []
for l in B["Supplementary_tables"]:
    if l.startswith("STAB "):
        s, t, d = [x.strip() for x in l[5:].split("|")]
        STAB.append((s, t, d))
assert len(STAB) == 40 and [s for s, _, _ in STAB] == ["S%d" % i for i in range(1, 41)], "supplementary table list must be S1 to S40"

# ---------------------------------------------------------------- 2. plain wording
TIER = {"T1": "Tier 1", "T2": "Tier 2", "T3": "Tier 3"}
GROUP = {"N": "N (position 1)", "ST": "S/T (position 3)", "both": "N + S/T", "XP": "X to P (position 2)"}
CFG = {"main": "Main", "strict": "Structure-matched strata", "nohyper": "Without hypermutated tumours",
       "union": "All primary-tumour MAF files", "norescue": "Without rescued rows",
       "uniform_restricted": "Position counts, context-model protein set", "sensB": "Germline variants in two or more resources",
       "calib_ctx": "Context model", "calib_uniform": "Uniform rate"}
FAMILY = {"primary": "Primary unit", "pooled_other": "Other pooled units", "cancer_T1_N": "Per cancer type, Tier 1 N",
          "cancer_T1_ST": "Per cancer type, Tier 1 S/T", "cancer_T1_both": "Per cancer type, Tier 1 N + S/T"}
SETQ1 = {"main": "All tumours", "nohyper": "Without hypermutated tumours"}
CAT = {"all_missense": "All missense mutations", "loss_T1": "Tier 1 sequon losses", "loss_T2": "Tier 2 sequon losses",
       "loss_T3": "Tier 3 sequon losses", "loss_any": "All sequon losses", "gain_plausible": "Gains, plausible regions",
       "gain_unknown": "Gains, unknown regions", "gain_implausible": "Gains, implausible regions", "gain_any": "All gains",
       "gain_plausible_or_unknown": "Gains, plausible or unknown regions"}
SIGSET = {"T1_loss_N": "Tier 1 loss, N", "T1_loss_XP": "Tier 1 loss, X to P", "T1_loss_ST": "Tier 1 loss, S/T",
          "loss_any": "All sequon losses", "gain_plausible": "Gains, plausible", "gain_unknown": "Gains, unknown",
          "bg_N": "Background N", "bg_ST": "Background S/T", "all_missense": "All missense"}
SIGCAT = {"SBS1": "SBS1 (clock-like)", "SBS5": "SBS5 (clock-like)", "APOBEC": "APOBEC", "Tobacco": "Tobacco", "UV": "UV",
          "MMR": "MMR deficiency", "POLE": "POLE", "other": "Other", "CpG_CtoT": "CpG C>T"}
LABEL = {"Enriched": "Changed, higher", "Depleted": "Changed, lower", "Higher": "Changed, higher", "Lower": "Changed, lower",
         "Not_clear": "Not clear", "Hotspot_calibration_FAILED": "Hotspot (calibration failed)",
         "Not hotspot_calibration_FAILED": "Not hotspot (calibration failed)",
         "Equivalent_positive_control_FAILED": "Equivalent (positive control failed)",
         "Equivalent_biology_control_FAILED": "Equivalent (germline biology control failed)",
         "Inconclusive, test not sensitive (C2 failed)": "Inconclusive, test not sensitive (null control failed)"}
LABCOLS = {"raw_label", "final_label", "label", "label_main", "label_nohyper", "label_raw", "descriptive_pattern_not_a_result_label",
           "q2_N_final_label", "q2_ST_final_label"}
MODEL = {"M1_3level": "Three groups (loss only, loss and control, control only)", "M2_any_loss": "Any loss against control only",
         "M3_adj_control": "Any loss, adjusted for control carrier"}
TERM = {"x_lo": "Loss only vs control only", "x_both": "Loss and control vs control only", "x_any": "Any loss vs control only",
        "x_L": "Any loss (adjusted for control carrier)", "x_C": "Control carrier (adjusted for loss)"}
ANALYSIS = {"main": "Main", "S1_strict_controls": "Strict controls", "S2_no_hypermutated": "Without hypermutated tumours",
            "S3_no_rescued_rows": "Without rescued rows", "S4_N_loss_only": "N losses only", "S5_ST_loss_only": "S/T losses only",
            "S6_OS_endpoint": "Overall survival", "S7_plus_DLBC_KICH": "Plus DLBC and KICH", "S8_unadjusted": "Unadjusted"}
SCORE = {"FoldX": "FoldX " + DDG + " (kcal/mol)", "AM": "AlphaMissense pathogenicity", "ESM": "ESM1b damage (-LLR)"}
TCGA_NAMES = {"ACC": "Adrenocortical carcinoma", "BLCA": "Bladder urothelial carcinoma", "BRCA": "Breast invasive carcinoma",
    "CESC": "Cervical squamous cell carcinoma and endocervical adenocarcinoma", "CHOL": "Cholangiocarcinoma",
    "COAD": "Colon adenocarcinoma", "DLBC": "Diffuse large B-cell lymphoma", "ESCA": "Oesophageal carcinoma",
    "GBM": "Glioblastoma multiforme", "HNSC": "Head and neck squamous cell carcinoma", "KICH": "Kidney chromophobe",
    "KIRC": "Kidney renal clear cell carcinoma", "KIRP": "Kidney renal papillary cell carcinoma", "LAML": "Acute myeloid leukaemia",
    "LGG": "Brain lower grade glioma", "LIHC": "Liver hepatocellular carcinoma", "LUAD": "Lung adenocarcinoma",
    "LUSC": "Lung squamous cell carcinoma", "MESO": "Mesothelioma", "OV": "Ovarian serous cystadenocarcinoma",
    "PAAD": "Pancreatic adenocarcinoma", "PCPG": "Phaeochromocytoma and paraganglioma", "PRAD": "Prostate adenocarcinoma",
    "READ": "Rectum adenocarcinoma", "SARC": "Sarcoma", "SKCM": "Skin cutaneous melanoma", "STAD": "Stomach adenocarcinoma",
    "TGCT": "Testicular germ cell tumours", "THCA": "Thyroid carcinoma", "THYM": "Thymoma",
    "UCEC": "Uterine corpus endometrial carcinoma", "UCS": "Uterine carcinosarcoma", "UVM": "Uveal melanoma"}

COLS = {
 "config": "Set-up", "tier": "Tier", "group": "Lost position", "scope": "Cancer type", "family": "Multiplicity family",
 "tested": "Tested", "n_informative": "Informative proteins", "expected_glyco": "Expected glycosite events",
 "n_events_glyco": "Glycosite events", "n_events_bg": "Background events", "n_positions_glyco": "Glycosite positions",
 "n_positions_bg": "Background positions", "log2rr": "log2 rate ratio", "rate_ratio": "Rate ratio",
 "lo95": "95% CI lower", "hi95": "95% CI upper", "lo90": "90% CI lower", "hi90": "90% CI upper", "p": "P value", "fdr": "FDR",
 "raw_label": "Label before control rules", "capability": "What the test can detect", "equivalent_allowed": "Equivalent label allowed",
 "final_label": "Label", "calibration": "Calibration of the expected rate", "residues": "Residues", "label": "Label",
 "accession": "UniProt accession", "gene": "Gene", "position": "Position", "ref": "Reference residue", "alt": "Mutant residue",
 "n_cases": "Tumours", "ev_g": "Glycosite events", "ev_b": "Background events",
 "exp_g_h0": "Expected glycosite events without selection", "n_strata": "Strata", "n_genes": "Genes",
 "boot_dropped": "Bootstrap samples dropped", "n_pos_glyco": "Glycosite positions", "n_pos_bg": "Background positions",
 "ci95_lo": "95% CI lower", "ci95_hi": "95% CI upper", "ci90_lo": "90% CI lower", "ci90_hi": "90% CI upper",
 "margin": "Margin", "median_d": "Median difference", "n_proteins": "Proteins", "project": "Cancer type", "n": "Mutations",
 "share": "Share", "category": "Signature group", "note": "Note", "difference": "Difference", "contrast": "Contrast",
 "vaf": "Variant allele fraction", "mechanism": "Lost position", "mutations": "Substitutions", "cancer": "CPTAC cancer type",
 "case": "CPTAC case", "in_all_files": "In all used MAF files", "hypermutated": "Hypermutated tumour", "HR": "Hazard ratio",
 "HR_lo": "HR 95% CI lower", "HR_hi": "HR 95% CI upper", "p_tost": "TOST P value (0.80 to 1.25)", "events": "Events",
 "zph_term": "Proportional hazards test P (term)", "zph_global": "Proportional hazards test P (global)",
 "b": "log hazard ratio", "se": "Standard error", "model": "Model", "term": "Contrast",
 "bailey_driver": "Known driver gene (Bailey et al. 2018)", "k": "Events", "expected": "Expected events",
 "symbol": "Gene", "length": "Protein length", "n_group": "Patients in the exposure group",
 "events_group": "Events in the exposure group", "analysis": "Analysis", "residue": "Residue", "rank": "Rank"}
DROP_ALWAYS = re.compile(r"^(run3_|p209_)|_run3$|^flag$|^tier_meaning$")

TOKENS = [(r"_biology_control_FAILED", " (germline biology control failed)"), (r"_positive_control_FAILED", " (positive control failed)"),
          (r"_calibration_FAILED", " (calibration failed)"), (r"\bNot_clear\b", "Not clear"),
          (r"reference_residue_differs", "reference residue differs"), (r"\(C2 failed\)", "(null control failed)"),
          (r"our flag is Amendment 2 \(plausible or unknown\)", "our flag is the plausible or unknown region class"),
          (r"loss_T12 ", "Tier 1 or Tier 2 losses, "), (r"\(tier T1 or T2\)", ""), (r"loss_T1 ", "Tier 1 losses, "),
          (r"loss_any ", "All sequon losses, "), (r"gain_any ", "All gains, "), (r"gain_pu ", "Plausible or unknown gains, "),
          (r"chosen-file", "selected file"), (r"TCGA projects", "TCGA cancer types"), (r"\btier T([123])\b", r"Tier \1"),
          (r"\bT([123])\b", r"Tier \1"), (r"TCGA-(?!CDR)", ""), (r"\s+$", "")]

def tok(x):
    if not isinstance(x, str): return x
    for a, b in TOKENS: x = re.sub(a, b, x)
    return x

def nice_col(c):
    if c in COLS: return COLS[c]
    s = c.replace("_", " ").strip()
    return s[:1].upper() + s[1:]

def clean(df, cols=None, vals=None, drop=(), keep=None):
    """vals: {original column: dict or function}; cols: extra renames; generic maps for config/tier/group/family/scope/labels."""
    d = df.copy()
    d = d[[c for c in d.columns if not DROP_ALWAYS.search(c) and c not in drop]]
    vals = dict(vals or {})
    gen = {"config": CFG, "tier": TIER, "group": GROUP, "family": FAMILY, "mechanism": GROUP}
    for c in d.columns:
        m = vals.get(c, gen.get(c))
        if m is None and c in LABCOLS: m = LABEL
        if m is None: continue
        if callable(m): d[c] = d[c].map(lambda v: m(v) if pd.notna(v) else v)
        else: d[c] = d[c].map(lambda v: m.get(v, v) if pd.notna(v) else v)
    if "scope" in d.columns: d["scope"] = d["scope"].map(lambda v: "Pooled" if v == "pooled" else v)
    for c in d.columns:
        if not (pd.api.types.is_numeric_dtype(d[c]) or pd.api.types.is_bool_dtype(d[c])): d[c] = d[c].map(tok)
    ren = {c: (cols or {}).get(c, nice_col(c)) for c in d.columns}
    d = d.rename(columns=ren)
    if keep: d = d[[k for k in keep if k in d.columns] + [k for k in d.columns if k not in keep]]
    return d

def ci90_cox(d):
    d = d.copy()
    d["HR 90% CI lower"] = np.exp(d["b"] - 1.6448536 * d["se"])
    d["HR 90% CI upper"] = np.exp(d["b"] + 1.6448536 * d["se"])
    return d

# ---------------------------------------------------------------- 3. supplementary sheets (list of blocks per sheet)
def s1():
    q = rt("P2_09_Q1_counts"); q = q[(q.set == "main") & (q.category == "all_missense")].copy()
    q["code"] = q.scope.str.replace("TCGA-", "", regex=False)
    q = q[q.code != "pooled"].sort_values("code")
    out = pd.DataFrame({"Cancer type": q.code, "Name": q.code.map(TCGA_NAMES), "Tumours": q.scope_cases_total,
                        "Tumours with missense mutations": q.scope_cases_with_mutations, "Missense mutations": q.n_mutations})
    t = rt("P2_09_Q1_counts"); t = t[(t.set == "main") & (t.category == "all_missense") & (t.scope == "pooled")]
    tot = pd.DataFrame({"Cancer type": ["All"], "Name": ["All 33 cancer types"], "Tumours": t.scope_cases_total.values,
                        "Tumours with missense mutations": t.scope_cases_with_mutations.values, "Missense mutations": t.n_mutations.values})
    return [(None, pd.concat([out, tot], ignore_index=True))]

def s2():
    return [(None, clean(rt("P2_09_Q1_counts"), vals={"set": SETQ1, "category": CAT},
        cols={"set": "Tumour set", "category": "Class", "n_mutations": "Mutations", "n_cases": "Tumours with the class",
              "n_unique_variants": "Unique substitutions", "scope_mutations_total": "Missense mutations in the scope",
              "scope_cases_with_mutations": "Tumours with missense mutations", "scope_cases_total": "Tumours in the scope"}))]

def s3():
    pos = {"1 (N)": "N (position 1)", "2 (X to P)": "X to P (position 2)", "3 (S/T)": "S/T (position 3)"}
    return [(None, clean(rt("P2_09_Q1_loss_position_in_sequon"), vals={"set": SETQ1, "sequon_position": pos},
        cols={"set": "Tumour set", "sequon_position": "Lost position", "n_sequon_losses": "Sequon losses"}))]

def s4():
    return [(None, clean(rt("P2_09_Q1_recurrent_variants"), vals={"category": CAT},
        cols={"category": "Class", "n_cases_without_hypermutated": "Tumours, without hypermutated tumours",
              "n_cancer_types": "Cancer types", "cancer_types_with_cases": "Tumours by cancer type", "rank_in_category": "Rank in class"}))]

def s5(): return [(None, clean(rt("P2_09_Q2_results")))]
def s6(): return [(None, clean(rt("P2_14_Q2ctx_results")))]
def s7():
    d = rt("P2_14_calibration"); d = d[d.run == d.run.min()]
    return [(None, clean(d, drop=("run",), cols={"config": "Expected-rate model", "num": "Events, high-expected background",
        "den": "Events, low-expected background", "eligible_proteins": "Eligible proteins", "n_inf": "Informative strata",
        "n_dropped": "Strata dropped"}))]

DN_SU = {("main", "strict", 1): "Main (strict mapping)", ("main", "strict", 3): "Main, second random seed",
         ("main", "window", 1): "Window mapping", ("nocaps", "strict", 1): "No mutation caps",
         ("nohyper", "strict", 1): "Without hypermutated tumours", ("nodrivers", "strict", 1): "Without driver genes"}
def s8():
    d = rt("P2_23_results"); d["su"] = [DN_SU.get((a, b, int(c))) for a, b, c in zip(d.setup, d.mapping, d.run)]
    d = d[d.su.notna()].copy(); d.insert(0, "Set-up", d.pop("su"))
    return [(None, clean(d, drop=("setup", "mapping", "run")))]
def s9():
    lo = {"(none: main run)": "None (main run)"}
    return [(None, clean(rt("P2_23_leave_one_out"), vals={"left_out": lo}, cols={"left_out": "Cancer type left out"}))]
def s10():
    return [(None, clean(rt("P2_23_absolute"), vals={"class": {"T1 glycosite": "Tier 1 glycosite", "background, same genes": "Background, same genes"}},
        cols={"class": "Positions", "observed": "Observed events", "expected": "Expected events", "dnds_raw": "dN/dS",
              "lo95_raw": "dN/dS 95% CI lower", "hi95_raw": "dN/dS 95% CI upper", "dnds_cal": "dN/dS divided by the synonymous calibration",
              "lo95_cal": "Calibrated 95% CI lower", "hi95_cal": "Calibrated 95% CI upper"}))]

CAPCOLS = {"run1_changed_right_1.5": "Planted 1.5, called Changed in the expected direction",
           "run1_changed_right_0.667": "Planted 0.67, called Changed in the expected direction",
           "run1_changed_right_1.25": "Planted 1.25, called Changed in the expected direction",
           "run1_changed_right_0.8": "Planted 0.8, called Changed in the expected direction",
           "run1_equivalent_1.0": "No planted effect, labelled Equivalent", "run1_equivalent_1.5": "Planted 1.5, labelled Equivalent",
           "run1_equivalent_0.667": "Planted 0.67, labelled Equivalent", "run1_median_null_log2rr": "Median log2 rate ratio, null sets",
           "run1_null_cover90": "Null sets, 90% CI covers 0", "run1_null_cover95": "Null sets, 95% CI covers 0",
           "n_tested": "Null sets tested", "n_changed": "Null sets called Changed", "n_equivalent": "Null sets labelled Equivalent",
           "fp_passed": "Null control passed"}
CYS_T = {"T1": "Metal-binding", "T2": "Active site", "T3": "Disulfide"}; CYS_G = {"N": "Cys", "ST": "His", "both": "Cys + His"}
GAIN_T = {"T1": "Plausible", "T2": "Plausible or unknown", "T3": "Implausible (comparison)"}
GAIN_G = {"N": "New N", "ST": "New S/T", "both": "New N + new S/T"}
def s11():
    parts = []
    for name, test, tv, gv in [("P2_08_unit_capability", "Position-count model", TIER, GROUP), ("P2_14_unit_capability", "Context-aware model", TIER, GROUP),
                               ("P2_15_unit_capability", "Germline glycosites", TIER, GROUP), ("P2_17_unit_capability", "Germline constrained Cys/His", CYS_T, CYS_G),
                               ("P2_19_unit_capability", "Sequon gains", GAIN_T, GAIN_G), ("P2_23_capability", "dNdScv", TIER, GROUP)]:
        d = rt(name); d = d[[c for c in d.columns if not c.startswith("run3_")]]
        d = clean(d, vals={"tier": tv, "group": gv}, cols=dict(CAPCOLS, tier="Tier or class", group="Position or residue"))
        d.insert(0, "Test", test); parts.append(d)
    return [(None, pd.concat(parts, ignore_index=True, sort=False))]

def s12():
    a = clean(rt("P2_14_spectra"), vals={"set": SETQ1}, cols={"set": "Tumour set", "n_cases": "Tumours", "channel": "Channel",
              "count": "SNVs", "rate": "Share"})
    g = clean(rt("P2_15_germline_spectra"), cols={"channel": "Channel", "opportunities": "Opportunities",
              "variants_A": "Variants, primary set", "rate_A": "Rate, primary set", "variants_B": "Variants, two or more resources",
              "rate_B": "Rate, two or more resources"})
    return [("Tumour somatic SNVs", a), ("Germline missense variants", g)]

def s13(): return [(None, clean(rt("P2_15_germline_results")))]
def s14():
    return [(None, clean(rt("P2_17_constrained_results"), vals={"tier": CYS_T, "group": CYS_G},
                         cols={"tier": "Constrained class", "group": "Residue"}, drop=("residues",)))]
def s15():
    return [(None, clean(rt("P2_16_shares_pooled"), vals={"set": SIGSET, "category": SIGCAT}, cols={"set": "Mutation set", "n": "Mutations"}))]
def s16():
    d = rt("P2_16_shares_by_type").rename(columns={k: v for k, v in SIGCAT.items()})
    return [(None, clean(d, vals={"set": SIGSET}, cols={"set": "Mutation set", "n": "Mutations"}))]
def s17():
    def con(v):
        for k in sorted(SIGSET, key=len, reverse=True): v = re.sub(r"\b" + k + r"\b", SIGSET[k], v)
        return v
    return [(None, clean(rt("P2_16_differences"), vals={"contrast": con, "category": SIGCAT},
                         cols={"difference": "Difference in share"}))]
def s18():
    return [(None, clean(rt("P2_13_fan_counts"), cols={"quantity": "Quantity (Fan et al. 2018)", "fan_2018": "Fan et al. 2018",
                         "ours": "This study", "our_definition": "Definition in this study"}))]
def s19():
    fs = {"LOG": "Loss of glycosylation (LOG)", "GOG": "Gain of glycosylation (GOG)", "LOG_text": "LOG, named in the text only"}
    st = {"ok": "Reference residue matches", "reference_residue_differs": "Reference residue differs"}
    gs = {"plausible": "Plausible", "unknown": "Unknown", "implausible": "Implausible"}
    return [(None, clean(rt("P2_13_fan_variants"), vals={"fan_set": fs, "status": st, "gain_state_observed": gs},
        cols={"fan_set": "Set of Fan et al.", "gene_fan": "Gene (Fan et al.)", "gene_uniprot": "Gene (UniProt)",
              "n_fan_cancer_types": "Cancer types (Fan et al.)", "fan_cancer_types": "Cancer types named by Fan et al.",
              "status": "Reference check", "residue_in_uniprot": "Residue in UniProt", "sequence_class": "Class",
              "class_agrees_with_fan": "Class agrees with Fan et al.", "lost_sequon_N": "Lost sequon, N position",
              "lost_sequon_tier": "Tier of the lost sequon", "gained_sequon_N": "New sequon, N position",
              "gain_state_observed": "Region class of the gain", "cases_exact_chosen_file": "Tumours with the substitution (selected file)",
              "cases_exact_without_hypermutated": "Tumours with the substitution, without hypermutated tumours",
              "cases_exact_any_primary_file": "Tumours with the substitution (any primary-tumour file)",
              "n_tcga_projects_exact": "Cancer types with the substitution", "tcga_projects_exact": "Tumours by cancer type",
              "other_substitutions_at_residue": "Other substitutions at the residue (tumours)", "label": "Result in TCGA",
              "recurrent_in_tcga": "Recurrent in TCGA"}))]
def s20():
    c = {"main_pos1_N": "N losses, all tumours", "main_pos2_XtoP": "X to P losses, all tumours", "main_pos3_ST": "S/T losses, all tumours",
         "main_ratio_3_to_1": "S/T to N ratio, all tumours", "nohyper_pos1_N": "N losses, without hypermutated tumours",
         "nohyper_pos2_XtoP": "X to P losses, without hypermutated tumours", "nohyper_pos3_ST": "S/T losses, without hypermutated tumours",
         "nohyper_ratio_3_to_1": "S/T to N ratio, without hypermutated tumours", "raw_position3_excess_reproduced": "Position-3 excess seen",
         "q2_N_observed": "N losses observed (position-count model)", "q2_N_expected": "N losses expected",
         "q2_N_final_label": "Label, N", "q2_ST_observed": "S/T losses observed (position-count model)", "q2_ST_expected": "S/T losses expected",
         "q2_ST_final_label": "Label, S/T", "observed_ratio_ST_to_N": "Observed S/T to N ratio", "expected_ratio_ST_to_N": "Expected S/T to N ratio",
         "log2_ratio_of_ratios": "log2 observed over expected ratio", "expected_from_background": "Excess expected from background rates"}
    return [(None, clean(rt("P2_13_fan_positions"), cols=c))]
def s21():
    st = {"loss": "Tier 1 sequon loss", "ctrl_d1": "Matched control, draw 1", "ctrl_d2": "Matched control, draw 2",
          "ctrl_strict": "Matched control, stricter matching", "c4_buried": "Method control, buried residue to D",
          "c4_lysine": "Method control, K to R"}
    return [(None, clean(rt("P2_11_foldx_ddg"), vals={"set": st, "role": {"N": "N (position 1)", "ST": "S/T (position 3)", "odd": "Other"}},
        drop=("model_number",), cols={"set": "Substitution set", "key": "Substitution key", "fragment": "AlphaFold fragment",
        "wt": "Wild-type residue", "mt": "Mutant residue", "role": "Lost position", "plddt": "pLDDT", "rsa": "Relative solvent accessibility",
        "fragment_model": "Fragment model", "loss_key": "Paired loss key", "ddg_mean": DDG + " mean (kcal/mol)", "ddg_sd": DDG + " SD",
        "run0": DDG + " run 1", "run1": DDG + " run 2", "run2": DDG + " run 3"}))]
UNIT6 = {1: "N + S/T losses (primary)", 2: "N losses", 3: "S/T losses", 4: "Without fragment models",
         5: "Stricter matching (same mutant residue)", 6: "Second control draw"}
def s22():
    d = rt("P2_18_A_results"); d = d[d.run == 1].copy(); d.insert(0, "Unit", d.unit.map(UNIT6))
    return [(None, clean(d, vals={"measure": SCORE}, drop=("run", "unit", "description"),
        cols={"measure": "Score", "n_pairs_unit": "Pairs", "n_scored": "Pairs scored", "median_d": "Median difference, loss minus control",
              "median_loss": "Median, losses", "median_controls": "Median, controls", "share_loss_above_cut": "Share of losses above the cut-off",
              "share_controls_above_cut": "Share of controls above the cut-off", "cut": "Cut-off", "n_unique_controls": "Unique controls"}))]
def s23():
    d = rt("P2_18_B_results"); d = d[d.run == 1]
    v = {"all19": "All 19 substitutions", "snv7": "Single-nucleotide substitutions", "strict19": "All 19, structure-matched"}
    return [(None, clean(d, vals={"score": SCORE, "variant": v}, drop=("run",),
        cols={"score": "Score", "variant": "Substitution set", "n_sequons": "Sequons", "median_d": "Median difference, sequon N minus background N",
              "null_min": "Null median, minimum", "null_median": "Null median, median", "null_max": "Null median, maximum",
              "null_share_abs_ge_observed": "Share of null medians at least as large as observed",
              "n_background_positions": "Background positions"}))]
def s24():
    d = rt("P2_18_C1_results").copy()
    d["split"] = [s + ", stricter matching" if u == 5 else s for s, u in zip(d.split, d.unit)]
    return [(None, clean(d, vals={"measure": SCORE}, drop=("unit", "status"),
        cols={"split": "Group", "measure": "Score", "n_group_pairs": "Loss pairs in the group", "n_group_proteins": "Proteins in the group",
              "median_group": "Median, group", "g_ci95_lo": "Group 95% CI lower", "g_ci95_hi": "Group 95% CI upper",
              "g_ci90_lo": "Group 90% CI lower", "g_ci90_hi": "Group 90% CI upper", "n_rest_pairs": "Loss pairs, rest",
              "n_rest_proteins": "Proteins, rest", "median_rest": "Median, rest", "diff_group_minus_rest": "Difference, group minus rest",
              "diff_ci95_lo": "Difference 95% CI lower", "diff_ci95_hi": "Difference 95% CI upper",
              "diff_ci90_lo": "Difference 90% CI lower", "diff_ci90_hi": "Difference 90% CI upper"}))]
WHY = {"no gene key in the proteome table": "Protein not in the global proteome", "carrier has no protein value": "Carrier has no protein value"}
def s25():
    return [(None, clean(rt("P2_22_loss_case_level"), vals={"why_not_usable": WHY}, drop=("group",),
        cols={"lost_T1_sites": "Lost Tier 1 sites", "n_used_files": "MAF files used",
              "excluded_trunc_same_gene": "Truncating mutation in the same gene (left out)", "in_loss_set": "In the analysis set",
              "y_adj": "Protein level, RNA-adjusted (log2)", "y_unadj": "Protein level, not adjusted (log2)",
              "pct_adj_in_reference": "Percentile among reference tumours", "why_not_usable": "Reason not usable"}))]
CPC = {"n_units": "Units (proteins or pairs)", "n_loss": "Loss pairs", "n_control": "Control pairs", "D_log2": "Difference (log2)",
       "ratio": "Ratio", "null": "Null control", "planted": "Planted control", "positive": "Positive control",
       "descriptive_pattern_not_a_result_label": "Pattern (descriptive, not a label)", "estimator": "Estimator", "outcome": "Outcome"}
EST = {"within": "Within protein", "pooled": "Pooled"}; OUTC = {"y_adj": "RNA-adjusted", "y_unadj": "Not adjusted"}
def s26():
    p = rt("P2_22_primary"); p.insert(0, "analysis", "Primary, within protein, RNA-adjusted")
    an = {"other estimator": "Pooled estimator", "unadjusted outcome": "Not RNA-adjusted", "nohyper": "Without hypermutated tumours",
          "mutations in all used files only": "Mutations in all used files", "mechanism N": "N losses only",
          "mechanism ST": "S/T losses only", "mechanism XP": "X to P losses only",
          "BOTH against CONTROL": "Carriers of a loss and a control mutation against control carriers"}
    bs = rt("P2_22_beside")
    out = [("Primary estimate", clean(p, vals={"estimator": EST, "outcome": OUTC}, cols=CPC)),
           ("Descriptive patterns beside the primary", clean(bs, vals={"analysis": lambda v: an.get(v, v.replace("cancer type", "Cancer type")),
                                                                        "estimator": EST, "outcome": OUTC}, cols=CPC))]
    c = rt("P2_22_controls"); isnp = c.control.astype(str).str.match(r"^(null|planted)")
    cv = lambda v: v.replace("null", "Null").replace("planted", "Planted") if isinstance(v, str) else v
    if isnp.any():
        out.append(("Null and planted controls (share of pseudo data sets by label)", clean(c[isnp], vals={"control": cv},
            cols={"control": "Control", "planted_d": "Planted shift (log2)", "repeats": "Pseudo data sets",
                  "share_Lower_or_D": "Share Lower", "share_Higher_or_ratio": "Share Higher", "share_Equivalent_or_lo95": "Share Equivalent",
                  "share_Not_clear_or_hi95": "Share Not clear", "share_not_tested_or_lo90": "Share not tested",
                  "median_D_or_label": "Median difference (log2)"})))
    if (~isnp).any():
        q = c[~isnp].dropna(axis=1, how="all")
        out.append(("Positive control (truncating carriers)", clean(q, vals={"control": {"positive control y_unadj within": "Truncating carriers, not adjusted (pre-specified)", "positive control y_adj within": "Truncating carriers, RNA-adjusted"}},
            drop=("planted_d", "repeats"), cols={"control": "Control", "share_Lower_or_D": "Difference (log2)",
            "share_Higher_or_ratio": "Ratio", "share_Equivalent_or_lo95": "95% CI lower", "share_Not_clear_or_hi95": "95% CI upper",
            "share_not_tested_or_lo90": "90% CI lower", "median_D_or_label": "Label"})))
    return out
def s27():
    return [(None, clean(rt("P2_22_glycopeptide"), cols={"lost_T1_site": "Lost Tier 1 site",
        "expected_log2_by_design": "Expected by design, log2(1 - VAF)", "protein_has_identical_ENSP": "Protein identical to the Ensembl reference",
        "n_glycoform_rows_at_site": "Glycoform rows at the site", "n_glycoforms_with_carrier_value": "Glycoforms with a carrier value",
        "observed_log2_vs_reference_median": "Observed glycopeptide, log2 against reference median",
        "observed_minus_expected": "Observed minus expected (log2)"}))]
PFIV = {"analysis": ANALYSIS, "model": MODEL, "term": TERM}
def s28():
    a = ci90_cox(rt("P2_12_results_long"))
    ps = rt("P2_12_primary_secondary")[["analysis", "model", "term", "family", "p_adj", "tost_adj", "label_raw", "label"]]
    a = a.merge(ps, on=["analysis", "model", "term"], how="left")
    fam = {"primary": "Primary", "secondary_holm3": "Secondary (Holm correction, 3 contrasts)"}
    return [(None, clean(a, vals=dict(PFIV, family=fam), drop=("note",), cols={"n": "Patients", "family": "Contrast role",
        "p": "Wald P value", "p_adj": "Holm-adjusted P", "tost_adj": "Holm-adjusted TOST P", "label_raw": "Label before correction"}))]
def s29():
    return [(None, clean(rt("P2_12_per_type"), cols={"type": "Cancer type", "n_control_only": "Control-only carriers",
        "n_any_loss": "Loss carriers", "events_any_loss": "Events in loss carriers", "p_bh": "FDR"}))]
def s30():
    return [(None, clean(ci90_cox(rt("P2_12_leave_one_type_out")), vals=PFIV, drop=("note",),
                         cols={"left_out": "Cancer type left out", "n": "Patients", "p": "Wald P value"}))]
def s31():
    return [(None, clean(rt("P2_19_gain_results"), vals={"tier": GAIN_T, "group": GAIN_G}, drop=("gain_type",),
                         cols={"tier": "Region class", "group": "Gain type"}))]
SL = {"mutation": "Substitution", "gain_type": "Gain type", "new_sequon_N": "New sequon, N position", "sequon_mutant": "New sequon",
      "projects": "Cancer types", "plddt_at_N": "pLDDT at the new N", "rsa_at_N": "RSA at the new N", "foldx_ddg": "FoldX " + DDG + " (kcal/mol)",
      "nearest_T1_distance": "Distance to the nearest Tier 1 site", "surfy_label_source": "SURFY label source",
      "netnglyc_potential": "NetNGlyc potential", "jury": "NetNGlyc jury agreement", "netnglyc_result": "NetNGlyc result",
      "glycosite_N": "Glycosite", "sequon": "Sequon", "glycan_glytoucan": "Glycan (GlyTouCan)", "glycoshape_id": "GlycoShape ID",
      "clash_status": "Re-Glyco result", "conformer": "Conformer"}
def s32(): return [(None, clean(rt("P2_19_netnglyc_scores"), drop=("id", "netnglyc_file"), cols=SL))]
def s33():
    d = rt("P2_19_reglyco_filled").copy()
    d["conformer"] = d.notes.astype(str).str.extract(r"conformer (model-\d+)", expand=False)
    def cs(v):
        if v == "built_clash_free": return "Clash-free"
        m = re.match(r"unresolved_(\d+)_clash_partners", str(v))
        return "Unresolved (%s clashing protein residues)" % m.group(1) if m else v
    return [(None, clean(d, vals={"clash_status": cs}, drop=("id", "mutant_pdb", "reglyco_seed", "reglyco_folder", "notes"), cols=SL))]
def s34():
    return [(None, clean(rt("P2_11_foldx_gains_summary"), vals={"gain_set": {"plausible": "Plausible regions", "plausible_or_unknown": "Plausible or unknown regions"},
        "models": {"all models": "All models", "without fragment models": "Without fragment models"}},
        cols={"gain_set": "Gains", "models": "Models", "n_planned": "Substitutions planned", "n_with_ddg": "Substitutions with " + DDG,
              "median_ddg": "Median " + DDG + " (kcal/mol)", "q25": "Lower quartile", "q75": "Upper quartile",
              "share_gt1": "Share above 1 kcal/mol", "share_gt3": "Share above 3 kcal/mol"}))]
def s35():
    d = rt("P2_20_T1_hotspots").sort_values("p_main")
    return [(None, clean(d, cols={"k_main": "Events, all tumours", "expected_main": "Expected events, all tumours", "obs_exp": "Observed over expected",
        "p_main": "P value, all tumours", "fdr_main": "FDR, all tumours", "label_main": "Label, all tumours",
        "k_nohyper": "Events, without hypermutated tumours", "p_nohyper": "P value, without hypermutated tumours",
        "fdr_nohyper": "FDR, without hypermutated tumours", "label_nohyper": "Label, without hypermutated tumours",
        "hypermutator_driven": "Driven by hypermutated tumours", "cases_by_cancer": "Tumours by cancer type"}))]
def s36(): return [(None, clean(rt("P2_20_background_hotspots").sort_values("p")))]
def s37():
    sm = rt("P2_20_summary")
    a = clean(sm, vals={"setup": {"main": "All tumours", "nohyper": "Without hypermutated tumours"}},
        cols={"setup": "Tumour set", "T1_positions": "Tier 1 positions", "T1_hotspots": "Tier 1 hotspots", "bg_positions": "Background positions",
              "bg_hotspots": "Background hotspots", "OR_hotspot": "Odds ratio, hotspots", "OR_lo": "95% CI lower", "OR_hi": "95% CI upper",
              "T1_nominal": "Tier 1 positions, nominal call", "bg_nominal": "Background positions, nominal call",
              "OR_nominal": "Odds ratio, nominal calls", "ORn_lo": "Nominal 95% CI lower", "ORn_hi": "Nominal 95% CI upper",
              "T1_obs": "Tier 1 events observed", "T1_exp": "Tier 1 events expected"})
    b = clean(rt("P2_20_positive_control"), cols={"called": "Called"})
    return [("Tier 1 against background", a), ("Known hotspots (positive control)", b)]
def s38():
    d = rt("P2_21_activedriver_results").sort_values("p")
    return [(None, clean(d, cols={"gene": "UniProt accession", "n_T1_sites": "Tier 1 sites", "n_mutations": "Missense mutations",
        "returned": "Returned by ActiveDriver", "n_T1_hotspots_P2_20": "Tier 1 single-position hotspots"}))]
def s39():
    a = clean(rt("P2_21_controls"), vals={"control": {"positive": "Positive control", "null": "Null, random non-sequon N"}},
              cols={"control": "Control", "rep": "Repeat", "fdr_below": "Proteins with FDR below 0.05"})
    b = clean(rt("P2_21_concordance"), cols={"ActiveDriver_enriched": "ActiveDriver call", "P2_20_hotspot": "Single-position hotspot in the protein",
                                             "Freq": "Proteins"})
    return [("Controls", a), ("Agreement with single-position hotspots", b)]

PY_DIST = {"sklearn": "scikit-learn", "yaml": "PyYAML", "docx": "python-docx", "Bio": "biopython", "PIL": "Pillow"}
def s40():
    prog = pd.DataFrame([
        ("UniProt", "Release 2026_03", "Proteome, features, variation file"), ("N-GlycositeAtlas", "HumanAll table", "Tier 1 evidence"),
        ("GDC masked MAF files", "GRCh38, GENCODE v36 transcripts", "TCGA and CPTAC somatic mutations"),
        ("AlphaFold Protein Structure Database", "v6, UP000005640", "Models, pLDDT, exposure"), ("FoldX", "5 (Windows)", DDG),
        ("AlphaMissense", "Amino-acid substitution scores", "Pathogenicity"), ("ESM1b", "LLR scores, all human isoforms", "Damage score"),
        ("NetNGlyc", "1.0 (web server)", "Sequon potential"), ("Re-Glyco", "Engine 0.2.6 (GlycoShape)", "Glycan placement"),
        ("SURFY", "In silico human surfaceome", "Surface annotation"), ("TCGA Clinical Data Resource", "Liu et al. 2018", "Progression-free interval"),
        ("Glycoprotein-Notebook", "Zhang et al. 2025 tables", "CPTAC protein and glycopeptide levels"),
        ("Driver genes", "Bailey et al. 2018, Table S1", "Driver annotation")], columns=["Program or resource", "Version or release", "Use"])
    sc = ROOT / "scripts"
    rp = {"ggplot2", "patchwork", "ggrepel", "scales", "ragg", "systemfonts", "survival", "dndscv"}
    for f in sc.glob("P2_*.R"):
        t = f.read_text(encoding="utf-8", errors="ignore")
        rp |= set(re.findall(r"(?:library|require|requireNamespace)\(\s*[\"']?([A-Za-z][A-Za-z0-9.]+)", t))
        rp |= set(re.findall(r"\b([A-Za-z][A-Za-z0-9.]{1,})::", t))
    rp -= {"base", "stats", "utils", "methods", "grDevices", "graphics", "tools", "parallel", "grid", "splines"}
    rrows = []
    if RSCRIPT.exists():
        code = "cat(R.version.string, '\\n'); for (p in c(%s)) cat(p, tryCatch(as.character(packageVersion(p)), error=function(e) 'not installed'), '\\n')" % ",".join("'%s'" % p for p in sorted(rp))
        try:
            o = subprocess.run([str(RSCRIPT), "-e", code], capture_output=True, text=True, timeout=120).stdout.splitlines()
            rrows.append(("R", o[0].replace("R version ", "").strip()))
            for l in o[1:]:
                s = l.split()
                if len(s) >= 2 and s[1] != "not": rrows.append((s[0], s[1]))
        except Exception as e:
            lg("Rscript call failed:", e); rrows.append(("R", "[confirm]"))
    else:
        rrows.append(("R", "4.6.0 [confirm: Rscript not found]"))
    pm = set()
    for f in sc.glob("P2_*.py"):
        for m in re.findall(r"^\s*(?:import|from)\s+([A-Za-z_][A-Za-z0-9_]*)", f.read_text(encoding="utf-8", errors="ignore"), re.M):
            pm.add(m)
    std = set(getattr(sys, "stdlib_module_names", ())) | {"__future__"}
    prow = [("Python", sys.version.split()[0])]
    import importlib.metadata as md
    for m in sorted(pm - std):
        if m.startswith("P2_"): continue
        try: prow.append((PY_DIST.get(m, m), md.version(PY_DIST.get(m, m))))
        except Exception: pass
    return [("Programs and data releases", prog), ("R and R packages", pd.DataFrame(rrows, columns=["Package", "Version"])),
            ("Python and Python packages", pd.DataFrame(prow, columns=["Package", "Version"]))]

BUILD = [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, s25,
         s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36, s37, s38, s39, s40]
assert len(BUILD) == 40

# ---------------------------------------------------------------- 4. workbook + scan
WT = re.compile(r"\bP2_\d+|\bQ[1-4]\b|\b[AB][1-4]\b|\bC1\b|\bT[123]\b|nohyper|sensB|uniform_restricted|norescue|nocaps|nodrivers|"
                r"\bunion\b|_FAILED|\bC-(?:N|Pl|P1|Cal)\b|\bscript|\bchat\b|\bitem\b|amendment|ctrl_|c4_|\bx_(?:lo|both|any|L|C)\b|"
                r"\bM[123]_|\bS[1-8]_|y_adj|y_unadj|\bgain_|\bloss_|\bbg_|calib_|[a-z]_[a-z]|\b3\d{7}\b", re.I)
scan = []
def scan_txt(where, cell, txt):
    if isinstance(txt, str):
        for m in WT.finditer(txt): scan.append((where, cell, txt[:120], m.group(0)))

wb = Workbook(); ws0 = wb.active; ws0.title = "README"
BOLD = Font(name="Times New Roman", bold=True); NORM = Font(name="Times New Roman"); TITLEF = Font(name="Times New Roman", bold=True, size=12)
nrows = {}
for (sid, title, desc), fn in zip(STAB, BUILD):
    ws = wb.create_sheet("Table " + sid)
    ws["A1"] = "Supplementary Table %s. %s" % (sid, title); ws["A1"].font = TITLEF
    ws["A2"] = desc; ws["A2"].font = NORM
    r = 4; total = 0; widths = {}
    try:
        blocks = fn()
    except Exception:
        tb = traceback.format_exc(); ERRORS.append((sid, tb)); lg("ERROR in", sid, "\n" + tb)
        ws["A4"] = "Not built (see log)"; nrows[sid] = 0; continue
    for sub, df in blocks:
        if sub:
            ws.cell(r, 1, sub).font = BOLD; scan_txt("Table " + sid, "subtitle", sub); r += 1
        hdr_row = r
        for j, c in enumerate(df.columns, 1):
            ws.cell(r, j, c).font = BOLD; ws.cell(r, j).alignment = Alignment(wrap_text=True, vertical="top")
            scan_txt("Table " + sid, "header", c); widths[j] = max(widths.get(j, 8), min(len(str(c)), 40))
        r += 1
        for row in df.itertuples(index=False):
            for j, v in enumerate(row, 1):
                if isinstance(v, (float, np.floating)) and (math.isnan(v) or math.isinf(v)): v = None
                elif isinstance(v, np.generic): v = v.item()
                cell = ws.cell(r, j, v); cell.font = NORM
                if isinstance(v, str):
                    scan_txt("Table " + sid, get_column_letter(j) + str(r), v); widths[j] = max(widths.get(j, 8), min(len(v), 45))
            r += 1
        total += len(df); r += 1
        if sub is None: ws.freeze_panes = ws.cell(hdr_row + 1, 1)
    for j, w in widths.items(): ws.column_dimensions[get_column_letter(j)].width = w + 2
    nrows[sid] = total
    scan_txt("Table " + sid, "title", title)
    lg("Table", sid, "rows", total)

TITLE = "Somatic loss of N-glycosylation sequons occurs at the rate expected from mutational processes across 33 tumour types"
readme = [("Supplementary Tables", ""), ("Paper", TITLE), ("Author", "Anishkiran Balasundar"), ("", ""),
          ("Notes", "Each sheet holds one table; the title is in row 1 and the column names in row 4. Set-up names are defined in Table 1 of the paper."),
          ("", "Labels (Equivalent; Changed, higher; Changed, lower; Inconclusive; Not tested; Not clear) follow the pre-specified rules in the Supplementary Methods (SM20). Labels with a note in brackets carry the result of a failed control."),
          ("", "Cancer types are TCGA codes (Supplementary Table S1) or CPTAC codes (Supplementary Table S25). CI, confidence interval; FDR, Benjamini-Hochberg false discovery rate; TOST, two one-sided tests; pLDDT, predicted local distance difference test; RSA, relative solvent accessibility; VAF, variant allele fraction."),
          ("", "Empty cells are values that do not apply or were not computed."), ("", "")]
r = 1
for a, b in readme:
    ws0.cell(r, 1, a).font = BOLD if a else NORM; ws0.cell(r, 2, b).font = NORM; r += 1
for j, h in enumerate(("Table", "Title", "Content", "Rows"), 1): ws0.cell(r, j, h).font = BOLD
r += 1
for sid, t, d in STAB:
    ws0.cell(r, 1, "Table " + sid).font = NORM; ws0.cell(r, 2, t).font = NORM; ws0.cell(r, 3, d).font = NORM; ws0.cell(r, 4, nrows.get(sid, 0)).font = NORM; r += 1
ws0.column_dimensions["A"].width = 12; ws0.column_dimensions["B"].width = 60; ws0.column_dimensions["C"].width = 70; ws0["A1"].font = TITLEF
XL = OUTT / "P2_Supplementary_Tables.xlsx"; wb.save(XL); lg("workbook written:", XL)

# ---------------------------------------------------------------- 5. Table 2 numbers
def f3(x): return ("%.3f" % x).replace("-", MINUS)
def f2(x): return ("%.2f" % x).replace("-", MINUS)
def pick(df, **kw):
    for k, v in kw.items(): df = df[df[k] == v]
    assert len(df) == 1, "expected one row for " + str(kw) + ", found " + str(len(df))
    return df.iloc[0]
def plain_lab(x): return tok(LABEL.get(x, x))
T2 = {}
try:
    for key, name in (("q2", "P2_09_Q2_results"), ("ctx", "P2_14_Q2ctx_results"), ("germ", "P2_15_germline_results"),
                      ("cys", "P2_17_constrained_results"), ("gain", "P2_19_gain_results")):
        q = pick(rt(name), config="main", scope="pooled", tier="T1", group="both")
        lab = plain_lab(q.final_label) + (" (germline biology control failed)" if key in ("q2", "ctx") else "")
        T2[key] = ("log2 rate ratio %s (%s to %s)" % (f3(q.log2rr), f3(q.lo90), f3(q.hi90)), lab)
    d = rt("P2_23_results"); q = pick(d, setup="main", mapping="strict", run=1, tier="T1", group="both")
    T2["dnds"] = ("log2 rate ratio %s (%s to %s)" % (f3(q.log2rr), f3(q.lo90), f3(q.hi90)), plain_lab(q.final_label))
    q = pick(rt("P2_16_differences"), contrast="T1_loss_ST minus T1_loss_N", category="CpG_CtoT")
    k = 100 if abs(q.difference) <= 1 else 1
    T2["sig"] = (("%+.1f (95%% CI %.1f to %.1f)" % (k * q.difference, k * q.lo95, k * q.hi95)).replace("-", MINUS), "Descriptive")
    q = pick(rt("P2_13_fan_positions"), tier="T1")
    T2["fan"] = ("%s (no CI)" % f3(q.log2_ratio_of_ratios), "Descriptive")
    a = rt("P2_18_A_results"); a = a[a.run == 1]
    for key, ms, unit in (("foldx", "FoldX", " kcal/mol"), ("am", "AM", ""), ("esm", "ESM", "")):
        q = pick(a, unit=1, measure=ms)
        T2[key] = ("%s%s (%s to %s)" % (f3(q.median_d), unit, f3(q.ci90_lo), f3(q.ci90_hi)), plain_lab(q.label))
    q = rt("P2_22_primary").iloc[0]
    T2["cptac"] = ("%s log2 (%s to %s)" % (f3(q.D_log2), f3(q.lo90), f3(q.hi90)), plain_lab(q.final_label))
    ps = rt("P2_12_primary_secondary")
    for key, term in (("pfi1", "x_lo"), ("pfi2", "x_any")):
        q = pick(ps, term=term)
        lo, hi = math.exp(q.b - 1.6448536 * q.se), math.exp(q.b + 1.6448536 * q.se)
        T2[key] = ("HR %s (%s to %s)" % (f2(q.HR), f2(lo), f2(hi)), plain_lab(q.label))
    q = pick(rt("P2_20_summary"), setup="main")
    T2["hot"] = ("OR %s (95%% CI %s to %s)" % (f2(q.OR_hotspot), f2(q.OR_lo), f2(q.OR_hi)), "Descriptive (calibration failed)")
    ad = rt("P2_21_activedriver_results"); nc = int((ad.label == "Enriched glycosite region").sum()); nt = len(ad)
    cn = rt("P2_21_controls"); nul = cn[cn.control == "null"]["share"].mean()
    T2["ad"] = ("%d of %s proteins (%.1f%%); null sets %.1f%%" % (nc, format(nt, ","), 100.0 * nc / nt, 100.0 * nul), "Descriptive")
except Exception:
    tb = traceback.format_exc(); ERRORS.append(("Table 2", tb)); lg("ERROR in Table 2\n" + tb)

t2rows = []
for l in B["Table_2"]:
    if l.startswith("T2ROW "):
        p = [x.strip() for x in l[6:].split("|")]
        key, an, unit, margin, ctrl, status, fig = p
        est, lab = T2.get(key, ("[not computed]", "[not computed]"))
        t2rows.append([an, unit, est, margin, lab, ctrl, status, fig])
T2H = [x.strip() for x in field("Table_2", "HEADER").split("|")]
pd.DataFrame(t2rows, columns=T2H).to_csv(OUTT / "P2_Table_2.tsv", sep="\t", index=False, encoding="utf-8")
for r_ in t2rows: lg("Table 2 |", r_[0], "|", r_[2], "|", r_[4])

# ---------------------------------------------------------------- 6. Word files
def setfont(doc, size):
    st = doc.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(size)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
def para(doc, bold_part, rest, size=None):
    p = doc.add_paragraph(); r1 = p.add_run(bold_part); r1.bold = True
    if rest: p.add_run(rest)
    if size:
        for r_ in p.runs: r_.font.size = Pt(size)
    return p
def table(doc, header, rows, size):
    t = doc.add_table(rows=1, cols=len(header)); t.style = "Table Grid"
    for j, h in enumerate(header):
        c = t.rows[0].cells[j]; c.text = ""; rr = c.paragraphs[0].add_run(h); rr.bold = True; rr.font.size = Pt(size)
    trPr = t.rows[0]._tr.get_or_add_trPr(); el = OxmlElement("w:tblHeader"); el.set(qn("w:val"), "true"); trPr.append(el)
    for row in rows:
        cells = t.add_row().cells
        for j, v in enumerate(row):
            cells[j].text = ""; rr = cells[j].paragraphs[0].add_run(str(v)); rr.font.size = Pt(size)
    return t

doc = docx.Document(); setfont(doc, 10)
sec = doc.sections[0]; sec.page_width, sec.page_height = Mm(210), Mm(297)
for s in ("left_margin", "right_margin", "top_margin", "bottom_margin"): setattr(sec, s, Mm(20))
para(doc, "Table 1. ", field("Table_1", "TITLE") + ".")
parts, cur = [], None
for l in B["Table_1"]:
    if l.startswith("HEADER "): cur = [[x.strip() for x in l[7:].split("|")], []]; parts.append(cur)
    elif l.startswith("ROW "): cur[1].append([x.strip() for x in l[4:].split("|")])
for i, (h, rows) in enumerate(parts):
    para(doc, "Part %s. " % "AB"[i], "Data sources" if i == 0 else "Set-ups", 10)
    table(doc, h, rows, 8.5); doc.add_paragraph()
for l in B["Table_1"]:
    if l.startswith("FOOT "): para(doc, "", l[5:], 8.5)
s2_ = doc.add_section(); s2_.orientation = WD_ORIENT.LANDSCAPE; s2_.page_width, s2_.page_height = Mm(297), Mm(210)
para(doc, "Table 2. ", field("Table_2", "TITLE") + ".")
table(doc, T2H, t2rows, 8)
doc.add_paragraph()
for l in B["Table_2"]:
    if l.startswith("FOOT "): para(doc, "", l[5:], 8.5)
D1 = OUTW / "P2_Tables_1_2.docx"; doc.save(D1); lg("tables written:", D1)

doc = docx.Document(); setfont(doc, 11)
doc.add_heading("Figure legends", level=1)
for f in FIGS[:8]:
    para(doc, "Figure %s. %s. " % (f.split("_")[1], field(f, "TITLE")), field(f, "LEGEND"))
doc.add_heading("Supplementary figure legends", level=1)
for f in FIGS[8:]:
    para(doc, "Supplementary Figure %s. %s. " % (f.split("_")[1], field(f, "TITLE")), field(f, "LEGEND"))
doc.add_heading("Alt text", level=1)
for f in FIGS:
    n = f.split("_")[1]
    para(doc, ("Figure %s. " if not n.startswith("S") else "Supplementary Figure %s. ") % n, field(f, "ALT"))
doc.add_heading("Table titles", level=1)
para(doc, "Table 1. ", field("Table_1", "TITLE") + "."); para(doc, "Table 2. ", field("Table_2", "TITLE") + ".")
doc.add_heading("Supplementary tables (one workbook, one sheet per table)", level=1)
for sid, t, d in STAB: para(doc, "Supplementary Table %s. " % sid, "%s. %s." % (t, d))
D2 = OUTW / "P2_Legends_AltText.docx"; doc.save(D2); lg("legends written:", D2)

# ---------------------------------------------------------------- 7. scan of the Word text and close
for f in FIGS:
    for k in ("TITLE", "LEGEND", "ALT"): scan_txt("Legends", f + " " + k, field(f, k).replace("Re-Glyco", "ReGlyco"))
for r_ in t2rows:
    for v in r_: scan_txt("Table 2", r_[0], v)
pd.DataFrame(scan, columns=["where", "cell", "text", "match"]).to_csv(OUTT / "P2_26_scan.tsv", sep="\t", index=False, encoding="utf-8")
lg("working-term scan:", len(scan), "hits (P2_26_scan.tsv);", len(set(s[3] for s in scan)), "distinct terms:",
   ", ".join(sorted(set(s[3] for s in scan)))[:600])
lg("sheets built:", sum(1 for s in nrows if nrows[s] > 0), "of 40; errors:", len(ERRORS))
lg("SqnPAN_26 done")
sys.exit(1 if ERRORS else 0)
