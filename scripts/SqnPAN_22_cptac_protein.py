# SqnPAN_22_cptac_protein.py - B2 protein-level check in CPTAC (Amendment 4, amended 2026-11-07, option A; added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-11-09; D1 A). Post hoc: no label of SqnPAN_09 to SqnPAN_21 changes.
# Question: in CPTAC tumours, is the gene-level protein amount lower in tumours that carry a T1 sequon-loss missense mutation in that protein than in tumours that carry
#   a control missense mutation in the same protein, given the RNA level? (Misfolding and degradation; an empirical protein-level check beside FoldX.)
# Data: Glycoprotein-Notebook cohorts (Zenodo 14019975), tumour columns (".T") of the gene-level proteome (log2 reference intensity) and RNA (RSEM UQ log2) tables:
#   BRCA, CCRCC (CCRCC_P zip; the CCRCC_G protein and RNA files are byte-identical), CRC, GBM, HNSCC, LSCC, LUAD (LUAD_P_PDC000149 zip), OV, PDAC, UCEC (UCEC_P_PDC000126 zip; the
#   UCEC_G and UCEC_P_PDC000441 protein and RNA files are byte-identical). Mutations: the 2,450 CPTAC masked somatic MAFs (GDC, open) with the GDC file map
#   GDC_CPTAC_masked_maf_case_map_2026-11-09.tsv: case submitter_id = Notebook column name without ".T" (1,017 of 1,030 tumours have a MAF, checked 2026-11-09).
#   Only files whose tumour aliquot (#tumor.aliquot line) is a "Primary Tumor" sample are used; the used files of a case are combined (union), with a flag for mutations
#   seen in every used file of the case. Hypermutated = more than 1,000 missense rows in a used file (flag).
# Mapping and classes: Missense_Mutation rows mapped to UniProt 2026_03 with map_missense (SqnPAN_06 rule 5A and Amendment 1), classify_change and the SqnPAN_06 tiers, unchanged.
#   Pair = one tumour x one protein. LOSS = at least one missense that removes a T1 sequon (N changed, X changed to P, or S/T at +2 changed to non-S/T).
#   CONTROL = at least one missense at a non-sequon N or S/T (res_role nonseqN or nonseqST; S<->T not counted; as SqnPAN_12 decision D3 A) and no LOSS mutation.
#   Both = BOTH (own group, not in the tests). Other missense are ignored. TRUNC (positive control) = Nonsense_Mutation, Frame_Shift_Del or Frame_Shift_Ins, by Hugo_Symbol.
#   A LOSS, CONTROL or BOTH pair whose tumour also has a TRUNC mutation in the same gene is left out (counted).
#   Protein set: proteins with at least one usable LOSS pair; CONTROL (and BOTH) pairs count only inside this set.
# Outcome: for cancer type c and gene g, reference tumours = tumours with a used MAF and no protein-changing mutation in g (Missense, Nonsense, Nonstop, Frame_Shift,
#   In_Frame, Splice_Site, Translation_Start_Site; by Hugo_Symbol and by the mapped UniProt gene). OLS protein = a + b RNA on the reference tumours with both values
#   (at least 20); outcome y = protein - (a + b RNA) in log2 units (below 0 = less protein than the RNA predicts). Unadjusted outcome: protein minus the median of the
#   reference tumours with a protein value (at least 20). Gene key: UniProt primary gene name if it is a row of the proteome table, else the MAF Hugo_Symbol.
#   A pair is usable if the carrier has protein and RNA values and the gene has a fit.
# Count step (section 1, counts only): usable LOSS, CONTROL, BOTH and TRUNC pairs (pooled and by cancer type), and the LOSS pairs whose protein also has a usable CONTROL
#   pair ("within-usable"). Route: within-usable LOSS >= 30 -> primary = within-protein; else usable LOSS >= 30 -> primary = pooled; else no test, case-level table only.
#   The counts and the route are written below as a dated line before section 2 can run (section 2 needs P2_22_counts_confirmed.txt).
# Primary (D1 A, within-protein): for every protein with >= 1 LOSS and >= 1 CONTROL pair, d_p = mean y(LOSS) - mean y(CONTROL), weight n_L n_C / (n_L + n_C); estimate
#   D = weighted mean (= the protein fixed-effects contrast). p: two-sided, 10,000 permutations of the LOSS/CONTROL labels within protein. CI: 2,000 bootstrap resamples of
#   proteins (percentile 95% and 90%). Pooled (beside; primary only if the route says so): D = mean y(LOSS) - mean y(CONTROL); p from 10,000 permutations within cancer
#   type; CI from 2,000 cluster bootstrap resamples of tumours.
# Labels (SqnPAN_08 thresholds unchanged, read as a protein ratio): Lower / Higher = p < 0.05 and |D| >= 0.32; else Equivalent = 90% CI inside -0.14 to 0.14 (ratio about
#   0.91 to 1.10); else Not_clear. One test, no BH family.
# Controls (section 2, judged before any LOSS outcome is read). Pseudo sets are drawn from CONTROL pairs only, sized like the real LOSS pairs per protein (within) or per
#   cancer type (pooled), at most n_C - 1 per protein or type:
#   Null (carrier-label permutation null): 500 pseudo sets, full label rule (2,000 permutations, 500 bootstrap resamples). PASSED if the share labelled Lower or Higher is
#     <= 0.07. If FAILED the primary reads "Inconclusive (null control failed)". The share labelled Equivalent is reported (capability of Equivalent).
#   Planted: the same kind of pseudo sets with -0.32, -0.5 and -1.0 added to the pseudo-LOSS outcomes, 200 each. PASSED if the share labelled Lower is >= 0.80 at -0.5.
#     If FAILED, Not_clear reads "Not_clear, test not sensitive".
#   Positive control (real data, all proteins): TRUNC against CONTROL pairs, unadjusted outcome, within-protein (unit = gene) if >= 30 within-usable TRUNC pairs, else
#     pooled. PASSED if labelled Lower. Not blind (nonsense-mediated decay lowers the protein). If FAILED every label gets the suffix _positive_control_FAILED.
#     The RNA-adjusted TRUNC contrast is descriptive.
# Beside (descriptive, no label changed): the other estimator; unadjusted outcome; nohyper (hypermutated tumours out); mutations seen in every used file only; LOSS by
#   mechanism (N, S/T, X to P) against all CONTROL pairs; BOTH against CONTROL; each cancer type (pooled, >= 5 LOSS pairs).
# Secondary glycopeptide readout (kept as secondary, no test): for LOSS pairs whose lost T1 site is in the cohort's cleaned glycoform table (ENSP of
#   GENCODE v42 whose sequence equals the UniProt sequence, so the positions agree): carrier value minus the median of the reference tumours (>= 5), averaged over the
#   site's glycoforms, beside the value expected by design, log2(1 - VAF) (heterozygous, equal expression assumed). The mutant glycopeptide cannot be identified (the
#   Notebook searched reference sequences), so a missing or lower glycopeptide is never read as loss of glycosylation.
# Not blind: the pooled Q2, FoldX and B1 results (Equivalent) were read before; no CPTAC protein value was read before.
# Seeds (block 33100000 to 33199999 reserved for B2): primary permutations 33100000, primary bootstrap 33110000, null 33120000 (+ r; bootstrap + 500 + r), planted
#   33130000, 33140000, 33150000 (d = -0.32, -0.5, -1.0), positive control 33160000 (+ 0 and + 1; bootstrap + 500), other estimator 33170000, beside 33180000
#   (+ i; bootstrap + 500 + i; cancer types + 100 + k and + 600 + k). Overlap check in the code.
# COUNTS (2026-11-09, section 1, counts only; recorded here before section 2 was run): 1,161 used MAF files, 1,017 tumours; usable LOSS pairs 141 in 133
#   proteins (mechanism N 67, S/T 69, X to P 5; by type BRCA 5, CCRCC 7, CRC 17, GBM 5, HNSCC 17, LSCC 22, LUAD 20, OV 6, PDAC 4, UCEC 38); usable CONTROL pairs in
#   the protein set 245; within-usable LOSS pairs 86 in 81 proteins -> route within (D1 A primary as approved); BOTH usable in set 1; TRUNC usable (unadjusted)
#   14,034, 8,882 in genes with a usable CONTROL pair -> positive control within-protein, against 8,648 CONTROL pairs.
# Usage (from the repository folder): python scripts/SqnPAN_22_cptac_protein.py 1     (counts only; compare with the COUNTS block above)
#   after the counts are confirmed (see README):  python scripts/SqnPAN_22_cptac_protein.py 2     (runs section 2 controls and section 3 real result)
# Logs: results\logs\P2_22_section<k>_<name>.txt. Outputs: data\processed P2_22_pairs.tsv.gz, P2_22_case_genes.tsv.gz, P2_22_controls_verdict.txt, gate files
#   P2_22_section1_ok to section3_ok; results\tables P2_22_counts.tsv, P2_22_controls.tsv, P2_22_primary.tsv, P2_22_beside.tsv, P2_22_loss_case_level.tsv,
#   P2_22_glycopeptide.tsv, P2_22_summary.tsv.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, re, csv, gzip, io, glob, math, time, zipfile, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_06_helpers as P
import SqnPAN_08_selection_controls as S8
import SqnPAN_14_context_helpers as C

BLOCK = 33100000
SEEDS = dict(perm=33100000, boot=33110000, null=33120000, pl1=33130000, pl2=33140000, pl3=33150000, pos=33160000, other=33170000, desc=33180000)
P1RAW = H.RAW   # LUAD protein zip and GENCODE v42 translations (shared with the companion LUAD study)
MAP = os.path.join(H.RAW, "GDC_CPTAC_masked_maf_case_map_2026-11-09.tsv"); MAP_MD5 = "cf1cd301ad6ed2746fb91c531479801e"
MANC = os.path.join(H.RAW, "CPTAC_masked_maf_manifest.txt"); MANC_MD5 = "b10990250fd4d89d77977ec4a1dba003"
MAFC = os.path.join(H.RAW, "maf_cptac")
G42 = os.path.join(P1RAW, "gencode.v42.pc_translations.fa.gz"); G42_MD5 = "2b83b9c1751f60a7cc3ad499973a80af"
COHORTS = [("BRCA", os.path.join(H.RAW, "BRC_P_PDC000121.zip"), "9266f7f1dbaab7272cf417577a9d9c00"),
           ("CCRCC", os.path.join(H.RAW, "CCRCC_P_PDC000128.zip"), "0e6efc6795a575c877629cdfa4da6c97"),
           ("CRC", os.path.join(H.RAW, "CRC_P_PDC000117.zip"), "42bbe35ac99c3c03836b077954b7609b"),
           ("GBM", os.path.join(H.RAW, "GBM_P_PDC000205.zip"), "9a91371d09d91859961e8e6cf6884969"),
           ("HNSCC", os.path.join(H.RAW, "HNSCC_P_PDC000222.zip"), "5ffc664b61c35dd4974cc6ce3f70bf37"),
           ("LSCC", os.path.join(H.RAW, "LSCC_P_PDC000232.zip"), "13a137f91ac9c7e58749910da1a6fd25"),
           ("LUAD", os.path.join(P1RAW, "LUAD_P_PDC000149.zip"), "8d29ffce58d8dacd8307df31c931b61c"),
           ("OV", os.path.join(H.RAW, "OVC_G_PDC000250.zip"), "3ca47a71c659bb24b1bec5a2ac71fa82"),
           ("PDAC", os.path.join(H.RAW, "PDAC_PandG_271_272.zip"), "a4e25b62509f40f5d3334070b64e4c7d"),
           ("UCEC", os.path.join(H.RAW, "UCEC_P_PDC000126.zip"), "39be02a0ab90544cfcdbd7088ac93f7a")]
PAIRS = os.path.join(H.PROC_OUT, "P2_22_pairs.tsv.gz"); CASEG = os.path.join(H.PROC_OUT, "P2_22_case_genes.tsv.gz")
CONFIRM = os.path.join(H.PROC_OUT, "P2_22_counts_confirmed.txt"); VERDICT = os.path.join(H.PROC_OUT, "P2_22_controls_verdict.txt")
OK1, OK2, OK3 = (os.path.join(H.PROC_OUT, f"P2_22_section{k}_ok.txt") for k in (1, 2, 3))
T = lambda n: os.path.join(H.TABLES, f"P2_22_{n}.tsv")
NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
          "Translation_Start_Site"}
TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins"}
MIN_REF = 20; MIN_ROUTE = 30; HYPER = 1000; MIN_GREF = 5; MIN_CT = 5
N_PERM = 10000; N_BOOT = 2000; N_NULL = 500; N_PL = 200; CP = 2000; CB = 500; NULL_MAX = 0.07; PLANT_MIN = 0.80
PLANT = [("pl1", -0.32), ("pl2", -0.5), ("pl3", -1.0)]
CHG = S8.CHG; EQ = S8.EQ; fmt = S8.fmt
NA = ("", "NA", "None", "nan", "NaN")
PCOLS = ["cancer", "case", "accession", "gene", "group", "mutations", "lost_T1_sites", "mechanism", "vaf", "in_all_files", "n_used_files", "hypermutated",
         "excluded_trunc_same_gene", "in_loss_set", "y_adj", "y_unadj", "pct_adj_in_reference", "why_not_usable"]
SCOLS = ["n_units", "n_loss", "n_control", "D_log2", "ratio", "lo95", "hi95", "lo90", "hi90", "p"]
NAN = float("nan")

# ---------------- small helpers ----------------
def need_ok(path, what):
    if not os.path.exists(path) or not open(path).read().startswith("OK"): raise RuntimeError(f"{what} has not finished ({os.path.basename(path)} missing)")

def write_ok(path):
    with open(path, "w") as o: o.write("OK " + time.strftime("%Y-%m-%d %H:%M:%S") + "\n")

def check_md5(log, pairs):
    bad = 0
    for path, rec in pairs:
        got = H.md5_file(path); log(f"md5 {os.path.basename(path)}: {got} ({'same as recorded' if got == rec else 'DIFFERENT, recorded ' + rec})"); bad += got != rec
    if bad: raise RuntimeError("an input md5 differs from the recorded value: stop and decide")

def tsv(path, header, rows):
    o = gzip.open(path, "wt", newline="", encoding="utf-8") if path.endswith(".gz") else open(path, "w", newline="", encoding="utf-8")
    with o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(header)
        for r in rows: w.writerow([fmt(x) for x in r])

def zmember(zf, key):
    m = [n for n in zf.namelist() if key in n and n.endswith("_GPNB_cleaned.tsv")]
    if len(m) != 1: raise RuntimeError(f"{zf.filename}: {len(m)} cleaned members contain '{key}' (expected 1)")
    return m[0]

def fnum(s): return float(s) if s not in NA else NAN

def read_matrix(zf, name):
    """{row id: float array over the tumour columns}, tumour ids (column name without .T), repeated row ids (first kept)."""
    out = {}; dup = 0
    with zf.open(name) as fb:
        f = io.TextIOWrapper(fb, encoding="utf-8", newline="")
        hdr = f.readline().rstrip("\r\n").split("\t"); tcol = [i for i, h in enumerate(hdr) if h.endswith(".T")]
        for line in f:
            v = line.rstrip("\r\n").split("\t")
            if v[0] in out: dup += 1; continue
            out[v[0]] = np.array([fnum(v[i]) if i < len(v) else NAN for i in tcol])
    return out, [hdr[i][:-2] for i in tcol], dup

def rna_vec(M, g):
    v = M["rna"].get(g)
    if v is None: return None
    a = np.full(len(M["pid"]), NAN); m = M["take"] >= 0; a[m] = v[M["take"][m]]; return a

def load_cohorts(log):
    mats = {}; case_ct = {}
    for ct, zp, _ in COHORTS:
        with zipfile.ZipFile(zp) as zf:
            prot, pid, d1 = read_matrix(zf, zmember(zf, "proteomics_gene_abundance"))
            rna, rid, d2 = read_matrix(zf, zmember(zf, "RNAseq_gene"))
        if len(set(pid)) != len(pid) or len(set(rid)) != len(rid): raise RuntimeError(f"{ct}: repeated tumour column")
        rix = {x: i for i, x in enumerate(rid)}; take = np.array([rix.get(x, -1) for x in pid])
        mats[ct] = dict(prot=prot, pid=pid, rna=rna, take=take, pidx={x: i for i, x in enumerate(pid)})
        for x in pid:
            if x in case_ct: raise RuntimeError(f"tumour {x} in two cohorts ({case_ct[x]}, {ct})")
            case_ct[x] = ct
        log(f"{ct}: proteome {len(prot)} genes x {len(pid)} tumours (repeated gene rows {d1}); RNA {len(rna)} genes x {len(rid)} tumours (repeated {d2}); "
            f"proteome tumours with RNA {int((take >= 0).sum())}")
    return mats, case_ct

def read_map():
    fmap = {}
    with open(MAP, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            al = {}
            for k in range(20):
                a = r.get(f"cases.0.samples.{k}.portions.0.analytes.0.aliquots.0.aliquot_id")
                if a: al[a] = r.get(f"cases.0.samples.{k}.sample_type") or ""
            fmap[r["file_name"]] = (r["cases.0.submitter_id"], r["cases.0.project.project_id"], al)
    return fmap

# ---------------- statistics ----------------
def subset_sums(v, m, n, rng):
    """n sums of m values drawn without replacement from v."""
    k = len(v)
    if m == 0: return np.zeros(n)
    if m >= k: return np.full(n, float(v.sum()))
    out = np.empty(n); step = max(1, 4000000 // k); i = 0
    while i < n:
        b = min(step, n - i); idx = np.argpartition(rng.random((b, k)), m - 1, axis=1)[:, :m]; out[i:i + b] = v[idx].sum(1); i += b
    return out

def stat_dict(D, p, bs, n_units, nl, nc):
    bs = np.where(np.isfinite(bs), bs, np.nan); ok = np.isfinite(bs).any()
    q = lambda a: float(np.nanquantile(bs, a)) if ok else NAN
    return dict(D=float(D), p=float(p), lo95=q(0.025), hi95=q(0.975), lo90=q(0.05), hi90=q(0.95), n_units=n_units, n_loss=nl, n_control=nc)

def within_stats(y, isL, unit, rp, rb, n_perm, n_boot):
    grp = collections.defaultdict(list)
    for i, u in enumerate(unit.tolist()): grp[u].append(i)
    d = []; w = []; nl = nc = 0; perm = np.zeros(n_perm)
    for u in sorted(grp):
        ix = np.array(grp[u]); l = isL[ix]; a = int(l.sum()); b = len(ix) - a
        if a == 0 or b == 0: continue
        v = y[ix]; wp = a * b / (a + b); d.append(v[l].mean() - v[~l].mean()); w.append(wp); nl += a; nc += b
        if n_perm: S = subset_sums(v, a, n_perm, rp); perm += wp * (S / a - (v.sum() - S) / b)
    if not d: return None
    d = np.array(d); w = np.array(w); W = w.sum(); D = (w * d).sum() / W
    p = (1 + np.sum(np.abs(perm / W) >= abs(D) - 1e-12)) / (1 + n_perm)
    bi = rb.integers(0, len(d), size=(n_boot, len(d))); bs = (w[bi] * d[bi]).sum(1) / w[bi].sum(1)
    return stat_dict(D, p, bs, len(d), nl, nc)

def pooled_stats(y, isL, strata, clus, rp, rb, n_perm, n_boot):
    a = int(isL.sum()); b = len(y) - a
    if a == 0 or b == 0: return None
    D = y[isL].mean() - y[~isL].mean(); Tt = y.sum(); S = np.zeros(n_perm)
    for s in sorted(set(strata.tolist())):
        ix = np.where(strata == s)[0]; m = int(isL[ix].sum())
        if m: S += subset_sums(y[ix], m, n_perm, rp)
    perm = S / a - (Tt - S) / b; p = (1 + np.sum(np.abs(perm) >= abs(D) - 1e-12)) / (1 + n_perm)
    cl, inv = np.unique(clus, return_inverse=True); k = len(cl); fl = isL.astype(float)
    sL = np.bincount(inv, weights=y * fl, minlength=k); nL = np.bincount(inv, weights=fl, minlength=k)
    sC = np.bincount(inv, weights=y * (1 - fl), minlength=k); nC = np.bincount(inv, weights=1 - fl, minlength=k)
    bs = np.empty(n_boot); step = max(1, 2000000 // k); i = 0
    while i < n_boot:
        bb = min(step, n_boot - i); c = rb.multinomial(k, np.full(k, 1.0 / k), size=bb).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"): bs[i:i + bb] = (c @ sL) / (c @ nL) - (c @ sC) / (c @ nC)
        i += bb
    return stat_dict(D, p, bs, k, a, b)

def estimate(kind, y, isL, unit, strata, clus, rp, rb, n_perm, n_boot):
    if kind == "within": return within_stats(y, isL, unit, rp, rb, n_perm, n_boot)
    return pooled_stats(y, isL, strata, clus, rp, rb, n_perm, n_boot)

def label(s):
    if s is None: return "Not tested (too few data)"
    if s["p"] < 0.05 and abs(s["D"]) >= CHG: return "Lower" if s["D"] < 0 else "Higher"
    if s["lo90"] > -EQ and s["hi90"] < EQ: return "Equivalent"
    return "Not_clear"

def srow(s):
    if s is None: return [""] * len(SCOLS)
    return [s["n_units"], s["n_loss"], s["n_control"], s["D"], 2.0 ** s["D"], s["lo95"], s["hi95"], s["lo90"], s["hi90"], s["p"]]

def sline(s):
    if s is None: return "not tested (too few data)"
    return (f"D {s['D']:.4f} log2 (ratio {2 ** s['D']:.3f}), 95% CI {s['lo95']:.4f} to {s['hi95']:.4f}, 90% CI {s['lo90']:.4f} to {s['hi90']:.4f}, p {s['p']:.4g}, "
            f"units {s['n_units']}, LOSS-side pairs {s['n_loss']}, CONTROL pairs {s['n_control']}")

def rngs(base, i): return np.random.default_rng(base + i), np.random.default_rng(base + 500 + i)

# ---------------- pairs table ----------------
def load_R():
    with gzip.open(PAIRS, "rt", newline="", encoding="utf-8") as f: rows = list(csv.DictReader(f, delimiter="\t"))
    R = {k: np.array([r[k] for r in rows], dtype=object) for k in ("cancer", "case", "accession", "gene", "group", "mechanism", "lost_T1_sites", "mutations")}
    for k in ("y_adj", "y_unadj", "vaf"): R[k] = np.array([fnum(r[k]) for r in rows])
    for k in ("in_loss_set", "in_all_files", "hypermutated", "excluded_trunc_same_gene"): R[k] = np.array([r[k] == "1" for r in rows])
    R["rows"] = rows; return R

def main_sel(R, yk="y_adj", grp="loss"):
    ok = ~np.isnan(R[yk]) & R["in_loss_set"] & ~R["excluded_trunc_same_gene"]
    return ok & (R["group"] == grp), ok & (R["group"] == "control")

def route_of(R):
    Lm, Cm = main_sel(R); cacc = set(R["accession"][Cm].tolist())
    nL = int(Lm.sum()); nW = sum(1 for a in R["accession"][Lm].tolist() if a in cacc)
    route = "within" if nW >= MIN_ROUTE else ("pooled" if nL >= MIN_ROUTE else "none")
    return route, dict(n_loss=nL, n_control=int(Cm.sum()), n_within=nW, n_within_proteins=len(set(R["accession"][Lm].tolist()) & cacc),
                       n_loss_proteins=len(set(R["accession"][Lm].tolist())))

# ======================= Section 1: counts =======================
def section1(log):
    check_md5(log, [(MAP, MAP_MD5), (MANC, MANC_MD5), (H.UNIPROT_MAIN, S8.MD5[H.UNIPROT_MAIN]), (S8.SITES, S8.MD5[S8.SITES])] + [(z, m) for _, z, m in COHORTS])
    C.check_seed_blocks(log, BLOCK, SEEDS.values(), n_rep=1000)
    U = H.load_uniprot_main(); revmap = P.build_revmap(U)
    tier_of = {(r["accession"], int(r["N_position"])): r["tier"] for r in csv.DictReader(open(S8.SITES, newline=""), delimiter="\t")}
    mats, case_ct = load_cohorts(log); fmap = read_map()
    man = {r["filename"] for r in csv.DictReader(open(MANC, newline=""), delimiter="\t")}
    paths = {os.path.basename(p): p for p in glob.glob(os.path.join(MAFC, "*", "*.maf.gz"))}
    log(f"MAF files on disk {len(paths)}; map rows {len(fmap)}; manifest rows {len(man)}; disk names equal map names {set(paths) == set(fmap)}; map names equal manifest {set(fmap) == man}")
    if set(paths) != set(fmap) or set(fmap) != man: raise RuntimeError("MAF files, map and manifest differ")
    need = ["Hugo_Symbol", "Variant_Classification", "SWISSPROT", "Protein_position", "Amino_acids", "all_effects", "t_depth", "t_alt_count", "Tumor_Sample_UUID"]
    FC = collections.Counter(); VC = collections.Counter(); MC = collections.Counter(); G = collections.Counter()
    used = collections.defaultdict(list); mis_max = collections.Counter(); nonsyn = collections.defaultdict(set)
    muts = collections.defaultdict(dict); trunc = collections.defaultdict(lambda: collections.defaultdict(set)); t0 = time.time()
    for fi, fn in enumerate(sorted(paths), 1):
        case, proj, al = fmap[fn]
        if case not in case_ct: FC["file not used: case not in a Notebook proteome table"] += 1; continue
        hdr = None; tal = None; nmis = 0; skip = False
        with gzip.open(paths[fn], "rt", encoding="utf-8", newline="") as f:
            for line in f:
                if line[0] == "#":
                    if line.startswith("#tumor.aliquot"):
                        q = line.split(); tal = q[1].strip() if len(q) > 1 else None
                    continue
                if hdr is None:
                    hdr = line.rstrip("\r\n").split("\t"); ix = {k: hdr.index(k) for k in need}
                    st = (al.get(tal) or "not in the map") if tal else "no #tumor.aliquot line"
                    if st != "Primary Tumor": FC["file not used: tumour aliquot sample type " + st] += 1; skip = True; break
                    FC["file used"] += 1; used[case].append(fn); continue
                v = line.rstrip("\r\n").split("\t"); vc = v[ix["Variant_Classification"]]; VC[vc] += 1
                if v[ix["Tumor_Sample_UUID"]] != tal: G["row tumour UUID differs from the #tumor.aliquot line"] += 1
                hugo = v[ix["Hugo_Symbol"]]
                if vc in NONSYN: nonsyn[case].add(hugo)
                if vc in TRUNC: trunc[case][hugo].add(fn)
                if vc != "Missense_Mutation": continue
                nmis += 1
                key, how = C.map_missense(v[ix["SWISSPROT"]], v[ix["Protein_position"]], v[ix["Amino_acids"]], hugo, v[ix["all_effects"]], U, revmap); MC[how] += 1
                if key is None: continue
                nonsyn[case].add(U[key[0]]["gene"])
                e = muts[case].setdefault(key, {"files": set(), "hugo": hugo, "depth": -1, "vaf": NAN}); e["files"].add(fn)
                try: dp = int(v[ix["t_depth"]]); al_ = int(v[ix["t_alt_count"]])
                except ValueError: dp, al_ = -1, 0
                if dp > e["depth"]: e["depth"] = dp; e["vaf"] = al_ / dp if dp > 0 else NAN
        if hdr is None and not skip: FC["file not used: no header line"] += 1
        if hdr is not None and not skip: mis_max[case] = max(mis_max[case], nmis)
        if fi % 500 == 0 or fi == len(paths): log(f"  {fi} of {len(paths)} MAF files, {round(time.time() - t0)} s")
    log("--- MAF files ---")
    for k, v in sorted(FC.items()): log(f"  {k}: {v}")
    log(f"  cases with at least one used file: {len(used)}; used files per case: {dict(sorted(collections.Counter(len(x) for x in used.values()).items()))}")
    log("--- variant classes in used files ---"); [log(f"  {k}: {v}") for k, v in VC.most_common()]
    log("--- missense mapping (map_missense) ---"); [log(f"  {k}: {v}") for k, v in sorted(MC.items())]
    log("--- guards (all must be 0) ---", dict(G) if G else "all 0")
    if G: raise RuntimeError("a guard is not 0")
    cache = {}
    def fit(ct, g):
        if (ct, g) in cache: return cache[(ct, g)]
        M = mats[ct]; p = M["prot"][g]; out = {"p": p, "why": ""}
        refm = np.array([(x in used) and (g not in nonsyn.get(x, ())) for x in M["pid"]]); pm = ~np.isnan(p)
        if (refm & pm).sum() >= MIN_REF: out["med"] = float(np.median(p[refm & pm]))
        r = rna_vec(M, g)
        if r is None: out["why"] = "gene not in the RNA table"
        else:
            ok = refm & pm & ~np.isnan(r)
            if ok.sum() < MIN_REF: out["why"] = "fewer than 20 reference tumours with protein and RNA"
            elif np.var(r[ok]) <= 0: out["why"] = "RNA constant in the reference tumours"
            else:
                b, a = np.polyfit(r[ok], p[ok], 1); res = p - (a + b * r); out.update(r=r, res=res, refres=np.sort(res[ok]))
        cache[(ct, g)] = out; return out
    def outcome(ct, g, j):
        if not g: return NAN, NAN, NAN, "no gene key in the proteome table"
        F = fit(ct, g); pv = F["p"][j]
        if np.isnan(pv): return NAN, NAN, NAN, "carrier has no protein value"
        yu = float(pv - F["med"]) if "med" in F else NAN
        if "res" not in F: return NAN, yu, NAN, F["why"]
        if np.isnan(F["r"][j]): return NAN, yu, NAN, "carrier has no RNA value"
        ya = float(F["res"][j]); return ya, yu, float(np.searchsorted(F["refres"], ya) / len(F["refres"])), ""
    rows = []; X = collections.Counter()
    for case in sorted(used):
        ct = case_ct[case]; M = mats[ct]; nf = len(used[case]); hyp = int(mis_max[case] > HYPER); j = M["pidx"][case]
        byacc = collections.defaultdict(list)
        for key, e in muts[case].items(): byacc[key[0]].append((key, e))
        for acc in sorted(byacc):
            L = []; Cl = []
            for (a, pos, ref, alt), e in sorted(byacc[acc], key=lambda t: t[0]):
                lost, gained = P.classify_change(U[a]["seq"], pos, ref, alt)
                if (lost, gained) != P.classify_literal(U[a]["seq"], pos, ref, alt): X["guard: classifier disagreement"] += 1
                lt1 = [q for q in lost if tier_of.get((a, q)) == "T1"]
                if ref == "N": role = ("seqN_" + tier_of[(a, pos)]) if (a, pos) in tier_of else "nonseqN"
                elif ref in "ST": role = ("seqST_" + tier_of[(a, pos - 2)]) if (a, pos - 2) in tier_of else "nonseqST"
                else: role = "other"
                if lt1:
                    mech = "N" if (ref == "N" and pos in lt1) else ("XP" if (alt == "P" and pos - 1 in lt1) else "ST")
                    L.append((f"{ref}{pos}{alt}", lt1, mech, e))
                elif role in ("nonseqN", "nonseqST") and not (ref in "ST" and alt in "ST"): Cl.append((f"{ref}{pos}{alt}", e))
                elif lost: X["missense removing a T2 or T3 sequon only (ignored)"] += 1
            if not L and not Cl: continue
            grp = "both" if (L and Cl) else ("loss" if L else "control")
            hugos = [t[-1]["hugo"] for t in L] + [t[-1]["hugo"] for t in Cl]
            g = U[acc]["gene"] if U[acc]["gene"] in M["prot"] else next((h for h in hugos if h in M["prot"]), "")
            excl = int(any(h in trunc[case] for h in hugos + [U[acc]["gene"]]))
            defin = [t[-1] for t in L] if L else [t[-1] for t in Cl]
            vaf = max((e["vaf"] for e in defin if not math.isnan(e["vaf"])), default=NAN)
            inall = int(any(len(e["files"]) == nf for e in defin))
            ya, yu, pc, why = outcome(ct, g, j)
            rows.append([ct, case, acc, g, grp, ";".join([t[0] for t in L] + [t[0] for t in Cl]), ";".join(str(q) for t in L for q in t[1]),
                         ";".join(sorted({t[2] for t in L})), vaf, inall, nf, hyp, excl, 0, ya, yu, pc, why])
        for h in sorted(trunc[case]):
            g = h if h in M["prot"] else ""; ya, yu, pc, why = outcome(ct, g, j)
            rows.append([ct, case, "", g, "trunc", "", "", "", NAN, int(len(trunc[case][h]) == nf), nf, hyp, 0, 0, ya, yu, pc, why])
    if X["guard: classifier disagreement"]: raise RuntimeError("classifier disagreement")
    lossset = {r[2] for r in rows if r[4] == "loss" and not r[12] and not math.isnan(r[14])}
    for r in rows:
        if r[4] in ("loss", "control", "both") and r[2] in lossset: r[13] = 1
    tsv(PAIRS, PCOLS, rows)
    tsv(CASEG, ["case", "cancer", "n_used_files", "missense_rows_max", "hypermutated", "nonsyn_genes"],
        [[c, case_ct[c], len(used[c]), mis_max[c], int(mis_max[c] > HYPER), ";".join(sorted(nonsyn[c]))] for c in sorted(used)])
    for k, v in sorted(X.items()): log(f"  {k}: {v}")
    R = load_R(); route, rc = route_of(R)
    out = []; log("--- tumours ---")
    for ct, _, _ in COHORTS:
        tum = mats[ct]["pid"]; wm = [x for x in tum if x in used]; nh = sum(1 for x in wm if mis_max[x] > HYPER)
        log(f"  {ct}: proteome tumours {len(tum)}, with a used MAF {len(wm)}, hypermutated {nh}"); out += [("tumours_proteome", ct, len(tum)), ("tumours_with_used_MAF", ct, len(wm)), ("hypermutated", ct, nh)]
    log("--- pairs (counts only) ---")
    gr = R["group"]; ua = ~np.isnan(R["y_adj"]); uu = ~np.isnan(R["y_unadj"]); ex = R["excluded_trunc_same_gene"]; ins = R["in_loss_set"]
    for g in ("loss", "control", "both", "trunc"):
        m = gr == g; why = collections.Counter(r["why_not_usable"] for r, k in zip(R["rows"], m) if k and r["why_not_usable"])
        line = (f"  {g}: pairs {int(m.sum())}, left out (TRUNC in the same gene) {int((m & ex).sum())}, usable adjusted {int((m & ua & ~ex).sum())}, "
                f"usable unadjusted {int((m & uu & ~ex).sum())}" + (f", usable adjusted inside the protein set {int((m & ua & ~ex & ins).sum())}" if g != "trunc" else ""))
        log(line); log(f"     not usable, reasons: {dict(why.most_common())}")
        out += [(f"{g}_pairs", "all", int(m.sum())), (f"{g}_excluded_trunc", "all", int((m & ex).sum())), (f"{g}_usable_adj", "all", int((m & ua & ~ex).sum())),
                (f"{g}_usable_unadj", "all", int((m & uu & ~ex).sum()))]
        for ct, _, _ in COHORTS:
            mc = m & (R["cancer"] == ct) & ua & ~ex & (ins if g != "trunc" else True); out.append((f"{g}_usable_adj" + ("_in_set" if g != "trunc" else ""), ct, int(mc.sum())))
    Lm, Cm = main_sel(R)
    log("  usable LOSS pairs by cancer type: " + ", ".join(f"{ct} {int((Lm & (R['cancer'] == ct)).sum())}" for ct, _, _ in COHORTS))
    log("  usable CONTROL pairs in the protein set by cancer type: " + ", ".join(f"{ct} {int((Cm & (R['cancer'] == ct)).sum())}" for ct, _, _ in COHORTS))
    log("  usable LOSS pairs by mechanism: " + str(dict(collections.Counter(R["mechanism"][Lm].tolist()))))
    Tm = (gr == "trunc") & uu; Cu = (gr == "control") & uu & ~ex; tg = set(R["gene"][Cu].tolist()); nTW = sum(1 for g in R["gene"][Tm].tolist() if g in tg)
    log(f"  positive control: usable (unadjusted) TRUNC pairs {int(Tm.sum())}, of these in genes with a usable CONTROL pair {nTW} -> "
        f"{'within-protein' if nTW >= MIN_ROUTE else 'pooled'}; CONTROL pairs (all proteins, unadjusted) {int(Cu.sum())}")
    log(f"ROUTE: usable LOSS pairs {rc['n_loss']} in {rc['n_loss_proteins']} proteins; usable CONTROL pairs in the protein set {rc['n_control']}; within-usable LOSS pairs "
        f"{rc['n_within']} in {rc['n_within_proteins']} proteins -> primary route: {route}")
    out += [(k, "all", v) for k, v in rc.items()] + [("trunc_within_usable", "all", nTW), ("route", "all", route)]
    tsv(T("counts"), ["item", "cancer", "value"], out)
    for p in (PAIRS, CASEG, T("counts")): log(f"written {p} | md5 {H.md5_file(p)}")
    write_ok(OK1); log("Section 1 done. Check the counts, then write the confirmation file (see README) before section 2 runs.")

# ======================= Section 2: controls =======================
def read_confirm():
    if not os.path.exists(CONFIRM): raise RuntimeError("P2_22_counts_confirmed.txt is missing: confirm the section 1 counts first (see README)")
    t = open(CONFIRM).read(); m = re.search(r"route=(\w+)", t)
    if not t.startswith("CONFIRMED") or not m: raise RuntimeError("P2_22_counts_confirmed.txt is not a confirmation")
    return m.group(1)

def pseudo_sets(route, R, Lm, Cm, rng):
    ci = np.where(Cm)[0]; key = R["accession"] if route == "within" else R["cancer"]
    need = collections.Counter(key[Lm].tolist()); byk = collections.defaultdict(list)
    for n, i in enumerate(ci): byk[key[i]].append(n)
    pl = np.zeros(len(ci), bool)
    for k2, m in sorted(need.items()):
        ix = byk.get(k2, [])
        if len(ix) < 2: continue
        for t in rng.choice(len(ix), size=min(m, len(ix) - 1), replace=False): pl[ix[t]] = True
    return ci, pl

def control_rep(route, R, Lm, Cm, base, r, d):
    rp, rb = rngs(base, r); ci, pl = pseudo_sets(route, R, Lm, Cm, rp)
    y = R["y_adj"][ci].copy(); y[pl] += d
    s = estimate(route, y, pl, R["accession"][ci], R["cancer"][ci], R["case"][ci], rp, rb, CP, CB); return s, label(s)

def read_verdict():
    return {l.split("\t")[0]: l.split("\t")[1].strip() for l in open(VERDICT) if "\t" in l}

def section2(log):
    need_ok(OK1, "section 1"); conf = read_confirm(); R = load_R(); route, rc = route_of(R)
    log(f"route recomputed from the pairs table: {route} ({rc}); confirmed route: {conf}")
    if route != conf: raise RuntimeError("the confirmed route differs from the pairs table")
    C.check_seed_blocks(log, BLOCK, SEEDS.values(), n_rep=1000)
    Lm, Cm = main_sel(R); V = {}; rows = []
    if route == "none":
        V["null"] = V["planted"] = "not run (no test)"
    else:
        for name, base, reps, d in [("null", SEEDS["null"], N_NULL, 0.0)] + [(f"planted {d}", SEEDS[k], N_PL, d) for k, d in PLANT]:
            t0 = time.time(); labs = collections.Counter(); Ds = []
            for r in range(1, reps + 1):
                s, lb = control_rep(route, R, Lm, Cm, base, r, d); labs[lb] += 1
                if s is not None: Ds.append(s["D"])
                if r == 10: log(f"  {name}: 10 repeats {time.time() - t0:.1f} s, estimate for {reps} about {(time.time() - t0) * reps / 10 / 60:.1f} min")
            sh = {k: labs[k] / reps for k in ("Lower", "Higher", "Equivalent", "Not_clear", "Not tested (too few data)")}
            log(f"  {name}: repeats {reps}, shares {sh}, median D {np.median(Ds) if Ds else NAN:.4f}, {time.time() - t0:.0f} s")
            rows.append([name, d, reps] + [sh[k] for k in ("Lower", "Higher", "Equivalent", "Not_clear", "Not tested (too few data)")] + [float(np.median(Ds)) if Ds else NAN])
            if name == "null": V["null"] = "PASSED" if sh["Lower"] + sh["Higher"] <= NULL_MAX else "FAILED"
            if d == -0.5: V["planted"] = "PASSED" if sh["Lower"] >= PLANT_MIN else "FAILED"
    ex = R["excluded_trunc_same_gene"]; gr = R["group"]
    for yk, idx in (("y_unadj", 0), ("y_adj", 1)):
        ok = ~np.isnan(R[yk]); Tm = ok & (gr == "trunc"); Cu = ok & (gr == "control") & ~ex
        tg = set(R["gene"][Cu].tolist()); nTW = sum(1 for g in R["gene"][Tm].tolist() if g in tg); kind = "within" if nTW >= MIN_ROUTE else "pooled"
        sel = Tm | Cu; rp, rb = rngs(SEEDS["pos"], idx)
        s = estimate(kind, R[yk][sel], Tm[sel], R["gene"][sel], R["cancer"][sel], R["case"][sel], rp, rb, N_PERM, N_BOOT); lb = label(s)
        log(f"  positive control TRUNC against CONTROL, {yk} ({'judged' if idx == 0 else 'descriptive'}), {kind}: {sline(s)} -> {lb}")
        rows.append([f"positive control {yk} {kind}", "", ""] + [""] * 5 + [lb]); rows[-1][3:8] = srow(s)[3:8] if s else [""] * 5
        if idx == 0: V["positive"] = "PASSED" if lb == "Lower" else "FAILED"
    tsv(T("controls"), ["control", "planted_d", "repeats", "share_Lower_or_D", "share_Higher_or_ratio", "share_Equivalent_or_lo95", "share_Not_clear_or_hi95",
                        "share_not_tested_or_lo90", "median_D_or_label"], rows)
    with open(VERDICT, "w") as o:
        o.write(f"P2_22 controls verdict {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        for k in ("null", "planted", "positive"): o.write(f"{k}\t{V[k]}\n")
    log("VERDICT:", V); log(f"written {VERDICT} and {T('controls')}"); write_ok(OK2)

# ======================= Section 3: real result =======================
def glyco_secondary(log, R):
    U = H.load_uniprot_main(); check_md5(log, [(G42, G42_MD5)])
    cg = {r["case"]: set(r["nonsyn_genes"].split(";")) for r in csv.DictReader(gzip.open(CASEG, "rt", newline="", encoding="utf-8"), delimiter="\t")}
    Li = [i for i in range(len(R["group"])) if R["group"][i] == "loss" and not R["excluded_trunc_same_gene"][i]]
    seq2acc = collections.defaultdict(set)
    for i in Li: seq2acc[U[R["accession"][i]]["seq"]].add(R["accession"][i])
    ensp_of = collections.defaultdict(set)
    with gzip.open(G42, "rt", encoding="utf-8") as f:
        hd = None; parts = []
        for line in itertools.chain(f, [">"]):
            if line.startswith(">"):
                if hd:
                    for a in seq2acc.get("".join(parts), ()): ensp_of[a].add(hd)
                hd = line[1:].split("|")[0].split(".")[0].strip(); parts = []
            else: parts.append(line.strip())
    log(f"glycopeptide readout: LOSS pairs {len(Li)}; proteins {len({R['accession'][i] for i in Li})}, with an identical GENCODE v42 ENSP {len(ensp_of)}")
    out = []
    for ct, zp, _ in COHORTS:
        W = collections.defaultdict(list)
        for i in Li:
            if R["cancer"][i] != ct: continue
            for st in [int(x) for x in R["lost_T1_sites"][i].split(";") if x]:
                for e in ensp_of.get(R["accession"][i], ()): W[(e, st)].append((i, st))
        got = collections.defaultdict(list); tid = []
        if W:
            with zipfile.ZipFile(zp) as zf, zf.open(zmember(zf, "nglycoform-glycosite")) as fb:
                f = io.TextIOWrapper(fb, encoding="utf-8", newline=""); hdr = f.readline().rstrip("\r\n").split("\t")
                tcol = [j for j, x in enumerate(hdr) if x.endswith(".T")]; tid = [hdr[j][:-2] for j in tcol]
                for line in f:
                    q = line.split("\t", 1)[0].split("_")
                    if len(q) < 4: continue
                    try: kk = (q[0], int(q[-2]))
                    except ValueError: continue
                    if kk in W:
                        v = line.rstrip("\r\n").split("\t")
                        for t in W[kk]: got[t].append(np.array([fnum(v[j]) if j < len(v) else NAN for j in tcol]))
        tix = {x: j for j, x in enumerate(tid)}
        for i in Li:
            if R["cancer"][i] != ct: continue
            for st in [int(x) for x in R["lost_T1_sites"][i].split(";") if x]:
                case = R["case"][i]; acc = R["accession"][i]; g1 = R["gene"][i]; g2 = U[acc]["gene"]
                refj = [tix[x] for x in tid if x in cg and g1 not in cg[x] and g2 not in cg[x]]
                vals = got.get((i, st), []); diffs = []
                for v in vals:
                    if case not in tix or np.isnan(v[tix[case]]): continue
                    rv = v[refj]; rv = rv[~np.isnan(rv)]
                    if len(rv) >= MIN_GREF: diffs.append(v[tix[case]] - float(np.median(rv)))
                vaf = R["vaf"][i]; exp = math.log2(1 - vaf) if (not math.isnan(vaf) and 0 <= vaf < 1) else NAN; obs = float(np.mean(diffs)) if diffs else NAN
                out.append([ct, case, acc, g2, st, R["mechanism"][i], R["mutations"][i], vaf, exp, int(acc in ensp_of), len(vals), len(diffs), obs, obs - exp,
                            int(R["in_all_files"][i]), int(R["hypermutated"][i])])
    tsv(T("glycopeptide"), ["cancer", "case", "accession", "gene", "lost_T1_site", "mechanism", "mutations", "vaf", "expected_log2_by_design", "protein_has_identical_ENSP",
                            "n_glycoform_rows_at_site", "n_glycoforms_with_carrier_value", "observed_log2_vs_reference_median", "observed_minus_expected",
                            "in_all_files", "hypermutated"], out)
    ob = np.array([r[12] for r in out]); dif = np.array([r[13] for r in out]); ex = np.array([r[8] for r in out])
    n_in = sum(1 for r in out if r[10] > 0); n_val = int(np.sum(~np.isnan(ob))) if len(out) else 0
    summ = [("glyco_loss_sites", len(out)), ("glyco_sites_in_table", n_in), ("glyco_sites_with_carrier_value", n_val),
            ("glyco_median_observed", float(np.nanmedian(ob)) if n_val else NAN), ("glyco_median_expected", float(np.nanmedian(ex[~np.isnan(ob)])) if n_val else NAN),
            ("glyco_median_observed_minus_expected", float(np.nanmedian(dif)) if np.isfinite(dif).any() else NAN)]
    for k, v in summ: log(f"  {k}: {fmt(v)}")
    return summ

def section3(log):
    need_ok(OK2, "section 2"); V = read_verdict(); R = load_R(); route, rc = route_of(R); Lm, Cm = main_sel(R); summ = [("route", route)] + list(rc.items())
    log("controls verdict:", V)
    tsv(T("loss_case_level"), PCOLS, [[r[k] for k in PCOLS] for r, g in zip(R["rows"], R["group"].tolist()) if g == "loss"])
    if route == "none":
        log("route none: no test; the case-level table of LOSS pairs is the result (descriptive)"); summ.append(("primary_label", "no test (case-level description only)"))
    else:
        sel = Lm | Cm
        s = estimate(route, R["y_adj"][sel], Lm[sel], R["accession"][sel], R["cancer"][sel], R["case"][sel],
                     np.random.default_rng(SEEDS["perm"]), np.random.default_rng(SEEDS["boot"]), N_PERM, N_BOOT)
        raw = label(s)
        if V.get("null") == "FAILED": final = "Inconclusive (null control failed)"
        else: final = "Not_clear, test not sensitive" if (raw == "Not_clear" and V.get("planted") == "FAILED") else raw
        if V.get("positive") == "FAILED": final += "_positive_control_FAILED"
        log(f"PRIMARY ({route}, RNA-adjusted protein, LOSS against CONTROL): {sline(s)} -> raw label {raw}; FINAL LABEL {final}")
        tsv(T("primary"), ["estimator", "outcome"] + SCOLS + ["raw_label", "final_label", "null", "planted", "positive"],
            [[route, "y_adj"] + srow(s) + [raw, final, V.get("null"), V.get("planted"), V.get("positive")]])
        summ += [("primary_D", s["D"] if s else NAN), ("primary_raw_label", raw), ("primary_final_label", final)]
        other = "pooled" if route == "within" else "within"; Lu, Cu = main_sel(R, "y_unadj"); Bm, _ = main_sel(R, grp="both")
        hy = R["hypermutated"]; ia = R["in_all_files"]; mech = R["mechanism"]
        B = [("other estimator", other, "y_adj", Lm, Cm, SEEDS["other"], 0), ("unadjusted outcome", route, "y_unadj", Lu, Cu, SEEDS["desc"], 1),
             ("nohyper", route, "y_adj", Lm & ~hy, Cm & ~hy, SEEDS["desc"], 2), ("mutations in all used files only", route, "y_adj", Lm & ia, Cm & ia, SEEDS["desc"], 3),
             ("mechanism N", route, "y_adj", Lm & (mech == "N"), Cm, SEEDS["desc"], 4), ("mechanism ST", route, "y_adj", Lm & (mech == "ST"), Cm, SEEDS["desc"], 5),
             ("mechanism XP", route, "y_adj", Lm & (mech == "XP"), Cm, SEEDS["desc"], 6), ("BOTH against CONTROL", route, "y_adj", Bm, Cm, SEEDS["desc"], 7)]
        for k, (ct, _, _) in enumerate(COHORTS):
            cm = R["cancer"] == ct
            if int((Lm & cm).sum()) >= MIN_CT: B.append((f"cancer type {ct}", "pooled", "y_adj", Lm & cm, Cm & cm, SEEDS["desc"], 100 + k))
        brows = []
        for name, kind, yk, a, b, base, i in B:
            sel2 = a | b; rp, rb = rngs(base, i); npm, nbt = (N_PERM, N_BOOT) if i < 100 else (2000, 1000)
            s2 = estimate(kind, R[yk][sel2], a[sel2], R["accession"][sel2], R["cancer"][sel2], R["case"][sel2], rp, rb, npm, nbt)
            log(f"  beside, {name} ({kind}, {yk}): {sline(s2)} [descriptive pattern {label(s2)}]"); brows.append([name, kind, yk] + srow(s2) + [label(s2)])
        tsv(T("beside"), ["analysis", "estimator", "outcome"] + SCOLS + ["descriptive_pattern_not_a_result_label"], brows)
    summ += glyco_secondary(log, R)
    tsv(T("summary"), ["item", "value"], summ)
    for n in ("loss_case_level", "primary", "beside", "glycopeptide", "summary"):
        if os.path.exists(T(n)): log(f"written {T(n)} | md5 {H.md5_file(T(n))}")
    write_ok(OK3)

SECTIONS = {"1": ("counts", section1), "2": ("controls", section2), "3": ("real_result", section3)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else ""
    if which == "1": H.run_section("P2_22", "1", "counts", section1)
    elif which in ("2", "3"): C.run_chain("P2_22", {k: SECTIONS[k] for k in ("2", "3")}, which)
    else: raise SystemExit("usage: SqnPAN_22_cptac_protein.py 1 (counts), then 2 (controls and real result) after the counts are confirmed")
# ---- END OF SCRIPT ----
