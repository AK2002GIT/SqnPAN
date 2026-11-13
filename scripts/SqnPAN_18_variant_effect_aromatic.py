# SqnPAN_18_variant_effect_aromatic.py - B1 AlphaMissense and ESM1b scores, C1 aromatic sequon split (Amendment 4; added after the primary analysis, no label changes)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (fixed 2026-11-01). Seed blocks: B1 33000000 to 33099999, C1 33400000 to 33499999.
# Scores: AlphaMissense pathogenicity as is (file line: CC BY-NC-SA 4.0; Zenodo says CC BY 4.0); ESM1b damage = -LLR (Brandes et al. 2023, Nature Genetics 55:1512; the Hugging Face Space
#   ntranoslab/esm_variants lists licence cc-by-nc-4.0 in its README header, no further licence text). The wild-type residue must equal UniProt 2026_03 (mismatches counted, left out).
# Section 1: structure check (zip list, names, coverage of the mapped proteins and of all loss and control keys, long proteins) and the margin per score = 0.10 x (median damage of
#   humsavar LP/P minus median of LB/B), computed BEFORE any TCGA pair is scored. Section 2 controls: humsavar LP/P vs LB/B AUROC >= 0.85 (AlphaMissense) and >= 0.80 (ESM1b), a failing
#   score is not interpreted; draw 1 vs draw 2 controls (median of the per-loss difference of the mean control damage, 95% CI inside +-margin); >= 95% of the primary pairs scored.
# Section 3 (a): the six SqnPAN_11 units unchanged (same pairs; bootstrap over proteins 2,000; runs 1 and 2 identical; run 3 other seed): d = damage(loss) - mean damage(its controls); 95% and
#   90% CI; labels Higher (95% CI above +margin) / Lower (95% CI below -margin) / Equivalent (90% CI inside +-margin) / Inconclusive, as FoldX with its 0.5 kcal/mol margin; shares above
#   AlphaMissense 0.564 and ESM1b LLR < -7.5 (Brandes 2023, text: "an LLR threshold of -7.5"); FoldX (taken from the SqnPAN_11 tables), AlphaMissense and ESM1b side by side; Spearman
#   correlation of the three d over the primary pairs (descriptive).
# Section 4 (b): proteome-wide, no TCGA. Per sequon N of tier T1 (T2 and T3 separate units): mean damage over all 19 substitutions minus the mean over the background N of the same
#   protein (N outside every N-X-S/T sequon and not at a sequon middle position; usable protein, pLDDT >= 70 at the position, as the Q2 background). Primary unit = all 19 substitutions;
#   beside it the single-nucleotide-reachable substitutions only (an N codon reaches D, H, I, K, S, T, Y whatever the codon: AAT or AAC), and a strict version matched also on pLDDT bin
#   (70 to 90, >= 90) and relative exposure bin (< 0.25, >= 0.25). Median paired difference, bootstrap over proteins (2,000), labels with the same rule and margin as (a) (the
#   pre-specification names no label rule for (b); this choice was made while writing the script and is logged), 50 within-protein label shuffles as null. Confounders: conservation and burial; both scores already
#   encode conservation, so a difference means predicted tolerance, not selection.
# Section 5 C1: on the existing SqnPAN_11 pairs (unit 1; unit 5 = strict matching as sensitivity for the primary split). Primary split: F/Y/W at N-2 vs the rest. Secondary: F only at N-2,
#   F/Y/W at N-3, F/Y/W at N-2 or N-3, F/Y/W at N-2 with T at +2, F at N-1 (Price 2011). FoldX, AlphaMissense and ESM1b d per group with 95% and 90% CI, difference of medians (group minus
#   rest) with bootstrap CI over proteins; no equivalence labels; a stratum needs >= 20 proteins and >= 30 pairs or it reads "Not tested". Limit for the text: FoldX does not model the glycan,
#   so C1 cannot test the aromatic-GlcNAc stabilisation of Culyba 2011 (Science 331:571); Y and W are an extension made here (the papers tested Phe).
# Seeds: bootstrap SeedSequence([base, unit or tier, variant, score]) with bases a12 33000000 (runs 1 and 2), a3 33020000 (run 3), b12 33040000, b3 33060000, null 33080000, C1 33400000 and 33420000.
# Usage (from the repository folder):  python scripts/SqnPAN_18_variant_effect_aromatic.py all [workers=6]     (or one section: 1 to 5)
# Expected run time: section 1 about 10 to 15 min (a guess: the AlphaMissense pass is a Python loop over 216 million lines; an awk pass took 1.9 min), sections 2 to 5 about 1 to 3 min each (guess).
# Logs: results\logs\P2_18_section<k>_<name>.txt.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, re, csv, gzip, math, time, zipfile, collections
import multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_14_context_helpers as C

ARGS = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
WORKERS = int(ARGS.get("workers", "6"))
AAS = "ACDEFGHIKLMNPQRSTVWY"; AAI = {a: i for i, a in enumerate(AAS)}
THREE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M",
         "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
PV = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
BLOCK_B = 33000000; SEEDS_B = dict(a12=33000000, a3=33020000, b12=33040000, b3=33060000, null=33080000)
BLOCK_C = 33400000; SEEDS_C = dict(c12=33400000, c3=33420000)
N_BOOT = 2000; N_SHUF = 50; MARGIN_FRAC = 0.10; FOLDX_MARGIN = 0.5
AM_CUT = 0.564; ESM_CUT = 7.5; FOLDX_CUT = 1.0
AUC_MIN = {"AM": 0.85, "ESM": 0.80}; MIN_PROT = 20; MIN_PAIRS_C1 = 30
SID = {"FoldX": 0, "AM": 1, "ESM": 2}; MEASURES = ["FoldX", "AM", "ESM"]; SCORES = ["AM", "ESM"]
SNV_N = [AAI[a] for a in "DHIKSTY"]
NOTINT = "Not interpreted (control failed)"; NOTTEST = "Not tested (too few data)"
ESM_ZIP = os.path.join(H.RAW, "esm1b", "ALL_hum_isoforms_ESM1b_LLR.zip"); ESM_MD5 = "12c4a5d2e8d2466196748da038cb0617"
AM_FILE = os.path.join(H.RAW, "alphamissense", "AlphaMissense_aa_substitutions.tsv.gz")
HUMSAVAR = os.path.join(H.RAW, "Uniprot_variation", "humsavar.txt")
PLAN = os.path.join(H.PROC, "P2_11_foldx_plan.tsv"); DDG = os.path.join(H.TABLES, "P2_11_foldx_ddg.tsv"); SITES = os.path.join(H.PROC, "P2_06_sites.tsv")
RESID = os.path.join(H.PROC, "P2_07_residues.tsv.gz")
ENG14 = os.path.join(H.PROC_OUT, "P2_14_engine.npz"); POS14 = os.path.join(H.PROC_OUT, "P2_14_positions.npz"); TMAP14 = os.path.join(H.TABLES, "P2_14_transcript_map.tsv")
CACHE = os.path.join(H.PROC_OUT, "P2_18_scores.npz"); MARGINS = os.path.join(H.PROC_OUT, "P2_18_margins.txt"); VERDICT = os.path.join(H.PROC_OUT, "P2_18_controls_verdict.txt")
def OKF(k): return os.path.join(H.PROC_OUT, f"P2_18_section{k}_ok.txt")
CORE_SETS = ("loss", "ctrl_d1", "ctrl_d2", "ctrl_strict")
UNITS = [(1, "primary: N + S/T losses, draw 1, all models", "ctrl_d1", ("N", "ST"), True), (2, "N losses, draw 1", "ctrl_d1", ("N",), True),
         (3, "S/T losses, draw 1", "ctrl_d1", ("ST",), True), (4, "N + S/T losses, draw 1, without fragment models", "ctrl_d1", ("N", "ST"), False),
         (5, "stricter matching (same mutant residue)", "ctrl_strict", ("N", "ST"), True), (6, "N + S/T losses, draw 2 (seed check of the control draw)", "ctrl_d2", ("N", "ST"), True)]
def P211_RUN(k): return os.path.join(H.TABLES, f"P2_11_comparison_run{k}.tsv")
def T_A(k): return os.path.join(H.TABLES, f"P2_18_A_run{k}.tsv")
def T_B(k): return os.path.join(H.TABLES, f"P2_18_B_run{k}.tsv")
def T_C(k): return os.path.join(H.TABLES, f"P2_18_C1_run{k}.tsv")

def need_ok(k, what):
    p = OKF(k)
    if not os.path.exists(p) or not open(p).read().startswith("OK"): raise RuntimeError(f"{what} has not finished ({os.path.basename(p)} missing)")
def write_ok(k, extra=""):
    with open(OKF(k), "w") as o: o.write("OK\n" + extra)
def write_tsv(path, hdr, rows):
    with open(path, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(hdr); w.writerows(rows)
def fnum(x, nd=4):
    try: v = float(x)
    except Exception: return str(x)
    return "nan" if math.isnan(v) else f"{v:.{nd}f}"
def offsets(U, accs):
    offs = np.zeros(len(accs) + 1, np.int64)
    for i, a in enumerate(accs): offs[i + 1] = offs[i] + len(U[a]["seq"])
    return offs

# ---------------- statistics ----------------
def rank_avg(x):
    u, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    cum = np.cumsum(cnt); return (cum - (cnt - 1) / 2.0)[inv.ravel()]
def auroc(pos, neg):
    r = rank_avg(np.concatenate([pos, neg])); n1, n0 = len(pos), len(neg)
    return float((r[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))
def spearman(a, b): return float(np.corrcoef(rank_avg(a), rank_avg(b))[0, 1])
def _wmed(ds, w):
    cw = np.cumsum(w); tot = int(cw[-1]) if len(cw) else 0
    if tot == 0: return np.nan
    lo = (tot + 1) // 2; hi = tot // 2 + 1
    return 0.5 * (ds[np.searchsorted(cw, lo, side="left")] + ds[np.searchsorted(cw, hi, side="left")])
def boot_sets(sets, m, seed, B=N_BOOT):
    """sets: list of (values, protein index 0..m-1). Bootstrap over proteins (the same resample for every set); returns B x len(sets) medians (median of the resampled data)."""
    prep = []
    for d, pl in sets:
        o = np.argsort(d, kind="stable"); prep.append((d[o], pl[o]))
    rng = np.random.default_rng(np.random.SeedSequence(seed)); out = np.full((B, len(sets)), np.nan)
    for b in range(B):
        cnt = np.bincount(rng.integers(0, m, m), minlength=m)
        for k, (ds, ps) in enumerate(prep):
            if len(ds): out[b, k] = _wmed(ds, cnt[ps])
    return out
def ci(col):
    c = col[~np.isnan(col)]
    if len(c) < 10: return (np.nan, np.nan), (np.nan, np.nan)
    return tuple(np.percentile(c, [2.5, 97.5])), tuple(np.percentile(c, [5, 95]))
def label_rule(c95, c90, margin):
    if any(math.isnan(x) for x in (c95[0], c95[1], c90[0], c90[1])): return "Inconclusive"
    if c95[0] > margin: return "Higher"
    if c95[1] < -margin: return "Lower"
    if c90[0] > -margin and c90[1] < margin: return "Equivalent"
    return "Inconclusive"
def local_index(accs_list):
    codes = {a: i for i, a in enumerate(sorted(set(accs_list)))}
    return np.array([codes[a] for a in accs_list], np.int64), len(codes)

# ---------------- humsavar ----------------
def read_humsavar(U, aidx, T):
    out = []; started = False
    with open(HUMSAVAR, encoding="utf-8") as f:
        for line in f:
            if line.startswith("_____"): started = True; continue
            if not started: continue
            p = line.split()
            if len(p) < 5: continue
            m = PV.match(p[3])
            if not m: continue
            T["humsavar missense rows"] += 1
            if p[4] not in ("LP/P", "LB/B"): T["humsavar rows with another category (US, not used)"] += 1; continue
            ai = aidx.get(p[1])
            if ai is None: T["humsavar rows on accessions not in the reviewed set (left out)"] += 1; continue
            ref = THREE.get(m.group(1)); alt = THREE.get(m.group(3)); pos = int(m.group(2)); sq = U[p[1]]["seq"]
            if ref is None or alt is None or ref == alt or pos > len(sq) or sq[pos - 1] != ref: T["humsavar rows whose reference residue differs from UniProt 2026_03 (left out)"] += 1; continue
            out.append((ai, pos, AAI[alt], 1 if p[4] == "LP/P" else 0)); T["humsavar LP/P and LB/B rows kept"] += 1
    return out

# ---------------- ESM1b extraction (pool workers) ----------------
_Z = {}
def _esm_init(path): _Z["z"] = zipfile.ZipFile(path)
def _esm_job(job):
    acc, seq, positions = job
    try: raw = _Z["z"].read(f"content/ALL_hum_isoforms_ESM1b_LLR/{acc}_LLR.csv")
    except KeyError: return acc, "no file", 0, False, {}
    lines = raw.decode("utf-8").split("\n")
    hdr = [t.strip().strip('"') for t in lines[0].split(",")][1:]; wt = []; okh = True
    for i, t in enumerate(hdr):
        pr = t.split()
        if len(pr) != 2 or pr[1] != str(i + 1): okh = False; wt.append("?")
        else: wt.append(pr[0])
    wts = "".join(wt); L = len(hdr)
    if not okh: return acc, "unexpected header", L, False, {}
    cells = {p: [float("nan")] * 20 for p in positions if p <= L and p <= len(seq) and wts[p - 1] == seq[p - 1]}; nrow = 0
    for ln in lines[1:]:
        if not ln.strip(): continue
        f = ln.split(","); lab = f[0].strip().strip('"')
        if lab not in AAI: return acc, "unexpected row label " + lab, L, False, {}
        a_ = AAI[lab]; nrow += 1
        for p, vec in cells.items():
            if wts[p - 1] != lab: vec[a_] = float(f[p])
    return acc, ("ok" if nrow == 20 else f"{nrow} rows"), L, wts == seq, cells

def read_am(log, U, accs, aidx, offs, keys):
    rowof = dict(zip(keys.tolist(), range(len(keys)))); V = np.full((len(keys), 20), np.nan, np.float32)
    aidx_b = {a.encode(): i for a, i in aidx.items()}; offl = offs.tolist(); seqb = [U[a]["seq"].encode() for a in accs]; aai = {ord(a): i for a, i in AAI.items()}
    T = collections.Counter(); seen = set(); last = None; n = 0; t0 = time.time()
    with gzip.open(AM_FILE, "rb") as f:
        for line in f:
            n += 1
            if line[0:1] == b"#" or line.startswith(b"uniprot_id"): continue
            parts = line.split(b"\t", 3)
            if len(parts) < 4: continue
            a, pv, am, _c = parts
            if a != last: last = a; seen.add(a)
            ai = aidx_b.get(a)
            if ai is None: continue
            pos = int(pv[1:-1]); row = rowof.get(offl[ai] + pos - 1)
            if row is None: continue
            sq = seqb[ai]
            if pos > len(sq) or sq[pos - 1] != pv[0]: T["requested cells whose reference residue differs from UniProt 2026_03 (left out)"] += 1; continue
            V[row, aai[pv[-1]]] = float(am); T["cells read"] += 1
            if n % 20000000 == 0: log(f"  AlphaMissense lines read {n} | {time.time() - t0:.0f} s")
    log(f"AlphaMissense file: {n} lines | {time.time() - t0:.0f} s | accessions in the file {len(seen)}, of them in the reviewed set {sum(1 for a in seen if a in aidx_b)}")
    for k in sorted(T): log(f"  {k}: {T[k]}")
    return V

# ---------------- score store ----------------
class Store:
    def __init__(self):
        z = C.load_npz(CACHE); self.keys = z["keys"]
        self.dam = {"AM": z["am"].astype(np.float32), "ESM": (-z["esm"]).astype(np.float32)}
        self.hum_gk = z["hum_gk"]; self.hum_alt = z["hum_alt"]; self.hum_cls = z["hum_cls"]
    def rows(self, gk):
        gk = np.asarray(gk, np.int64); i = np.minimum(np.searchsorted(self.keys, gk), len(self.keys) - 1); return np.where(self.keys[i] == gk, i, -1)
    def val(self, score, gk, alt):
        r = self.rows(gk); v = np.full(len(r), np.nan, np.float32); ok = r >= 0
        v[ok] = self.dam[score][r[ok], np.asarray(alt)[ok]]; return v
    def hum(self, score):
        v = self.val(score, self.hum_gk, self.hum_alt); return v, self.hum_cls
def read_margins():
    out = {}
    for line in open(MARGINS):
        p = line.rstrip("\n").split("\t")
        if p[0] in SCORES: out[p[0]] = float(p[1])
    return out
def read_verdict():
    out = {}
    for line in open(VERDICT):
        p = line.split()
        if len(p) >= 2 and p[0] in SCORES: out[p[0]] = p[1] == "PASSED"
    return out

# ---------------- pairs of SqnPAN_11 ----------------
def load_pairs(log, U, aidx, offs, st):
    with open(PLAN, newline="", encoding="utf-8") as f: plan = list(csv.DictReader(f, delimiter="\t"))
    dd = {}
    with open(DDG, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["set"] in CORE_SETS and r["ddg_mean"] != "": dd[r["key"]] = float(r["ddg_mean"])
    # note: P2_11_foldx_ddg.tsv holds ddG rounded to 4 decimals; the FoldX rows of section 3 come from the SqnPAN_11 comparison tables, the FoldX d used here (correlation, C1) from the rounded values
    info = {}; losses = {}; ctrl = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in plan:
        if r["set"] not in CORE_SETS: continue
        info[r["key"]] = (r["accession"], int(r["position"]), r["wt"], r["mt"], r["role"], r["fragment_model"])
        if r["set"] == "loss": losses[r["key"]] = r
        else: ctrl[r["set"]][r["loss_key"]].append(r["key"])
    keys = sorted(info); gk = []; alt = []; bad = 0
    for k in keys:
        acc, pos, wt, mt, _r, _f = info[k]; ai = aidx.get(acc)
        if ai is None or pos > len(U[acc]["seq"]) or U[acc]["seq"][pos - 1] != wt: gk.append(-1); alt.append(0); bad += 1
        else: gk.append(int(offs[ai]) + pos - 1); alt.append(AAI[mt])
    gk = np.array(gk, np.int64); alt = np.array(alt, np.int64); VAL = {"FoldX": dict(dd)}
    for s in SCORES:
        v = st.val(s, np.where(gk >= 0, gk, 0), alt); v[gk < 0] = np.nan; VAL[s] = {k: float(x) for k, x in zip(keys, v) if not np.isnan(x)}
    log(f"plan substitutions {len(keys)} | reference residue differs from UniProt 2026_03: {bad} | with FoldX ddG {len(dd)} | with AlphaMissense {len(VAL['AM'])} | with ESM1b {len(VAL['ESM'])}")
    PAIRS = {}
    for u, desc, cs, kinds, allm in UNITS:
        pr = []
        for k, L in losses.items():
            if L["role"] not in kinds or (not allm and L["fragment_model"] == "yes") or k not in dd: continue
            c_ = [c for c in ctrl[cs].get(k, []) if c in dd]
            if c_: pr.append((k, c_))
        PAIRS[u] = pr
    with open(P211_RUN(1), newline="", encoding="utf-8") as f: ref = {int(r["unit"]): r for r in csv.DictReader(f, delimiter="\t")}
    for u, desc, cs, kinds, allm in UNITS:
        npair = len(PAIRS[u]); nprot = len({info[k][0] for k, _c in PAIRS[u]})
        ok = npair == int(ref[u]["n_pairs"]) and nprot == int(ref[u]["n_proteins"])
        log(f"  unit {u}: {npair} pairs in {nprot} proteins (SqnPAN_11 table: {ref[u]['n_pairs']} in {ref[u]['n_proteins']}) {'equal' if ok else 'DIFFERENT'}")
        if not ok: raise RuntimeError("pairs differ from the SqnPAN_11 comparison table: stop")
    return info, VAL, PAIRS, ref
def pair_arrays(pairs, V, info):
    d = []; acc = []; lv = []; cu = {}; keys = []
    for k, cs in pairs:
        x = V.get(k)
        if x is None: continue
        cv = [V[c] for c in cs if c in V]
        if not cv: continue
        d.append(x - sum(cv) / len(cv)); acc.append(info[k][0]); lv.append(x); keys.append(k)
        for c in cs:
            if c in V: cu[c] = V[c]
    return np.array(d, float), acc, np.array(lv, float), cu, keys

# ======================= Section 1: structure check, extraction, margins =======================
def section1(log):
    log("Pre-specification: header of this script (fixed 2026-11-01). Added after the primary analysis.")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__} | workers {WORKERS}")
    U = H.load_uniprot_main(); accs = sorted(U); aidx = {a: i for i, a in enumerate(accs)}; offs = offsets(U, accs)
    got = H.md5_file(ESM_ZIP); log(f"md5 ESM1b zip: {got} ({'same as recorded' if got == ESM_MD5 else 'DIFFERENT, recorded ' + ESM_MD5})")
    if got != ESM_MD5: raise RuntimeError("ESM1b zip md5 differs: stop")
    log(f"md5 AlphaMissense file: {H.md5_file(AM_FILE)} (v3; the file header says CC BY-NC-SA 4.0)")
    z = zipfile.ZipFile(ESM_ZIP); names = z.namelist(); pat = re.compile(r"^content/ALL_hum_isoforms_ESM1b_LLR/(.+)_LLR\.csv$"); ids = {}; other = []
    for nm in names:
        m = pat.match(nm)
        if m: ids[m.group(1)] = nm
        else: other.append(nm)
    canon = {i for i in ids if "-" not in i}
    log(f"zip: {len(names)} entries | per-protein files {len(ids)} (canonical IDs {len(canon)}, isoform IDs {len(ids) - len(canon)}) | other entries {len(other)}: {other[:3]}")
    mapped = set()
    with open(TMAP14, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"): mapped.add(r["accession"])
    log(f"UniProt 2026_03 reviewed proteins {len(accs)} | with a canonical ESM1b file {sum(1 for a in accs if a in canon)} | mapped proteins of SqnPAN_14 {len(mapped)}, of them with a file {sum(1 for a in mapped if a in canon)}")
    long_ = [a for a in accs if len(U[a]["seq"]) > 1022]; log(f"proteins longer than 1,022 residues {len(long_)}, of them with a canonical file {sum(1 for a in long_ if a in canon)}")
    T = collections.Counter(); hum = read_humsavar(U, aidx, T)
    for k in sorted(T): log(f"  {k}: {T[k]}")
    Eg = C.load_npz(ENG14); Pz = C.load_npz(POS14); isN = Eg["P_res"] == 0
    gkN = offs[Eg["P_prot"][isN].astype(np.int64)] + Pz["P_posn"][isN].astype(np.int64) - 1
    with open(PLAN, newline="", encoding="utf-8") as f: plan = [r for r in csv.DictReader(f, delimiter="\t") if r["set"] in CORE_SETS]
    gkP = []; badw = 0
    for r in plan:
        ai = aidx.get(r["accession"]); pos = int(r["position"])
        if ai is None or pos > len(U[r["accession"]]["seq"]) or U[r["accession"]]["seq"][pos - 1] != r["wt"]: badw += 1; continue
        gkP.append(int(offs[ai]) + pos - 1)
    gkH = [int(offs[a]) + p - 1 for a, p, _alt, _c in hum]
    keys = np.unique(np.concatenate([gkN, np.array(gkP, np.int64), np.array(gkH, np.int64)]))
    log(f"requested positions: all N of the Q2 position list {len(gkN)} | plan positions {len(set(gkP))} (reference residue differs: {badw}) | humsavar {len(set(gkH))} | together {len(keys)}")
    ai_arr = np.searchsorted(offs, keys, side="right") - 1; pos_arr = keys - offs[ai_arr] + 1
    cut = np.flatnonzero(np.diff(ai_arr)) + 1; starts = np.concatenate([[0], cut]); ends = np.concatenate([cut, [len(keys)]])
    jobs = []; nofile = 0; nofile_pos = 0
    for s_, e_ in zip(starts.tolist(), ends.tolist()):
        a = accs[int(ai_arr[s_])]
        if a in canon: jobs.append((a, U[a]["seq"], pos_arr[s_:e_].tolist()))
        else: nofile += 1; nofile_pos += e_ - s_
    log(f"proteins with requested positions {len(starts)} | without a canonical ESM1b file {nofile} ({nofile_pos} positions) | jobs {len(jobs)}")
    rowof = dict(zip(keys.tolist(), range(len(keys)))); ES = np.full((len(keys), 20), np.nan, np.float32)
    st_ = collections.Counter(); t0 = time.time(); done = 0; long_ok = 0; long_n = 0
    with mp.get_context("spawn").Pool(WORKERS, initializer=_esm_init, initargs=(ESM_ZIP,)) as pool:
        for acc, status, L, seq_eq, cells in pool.imap_unordered(_esm_job, jobs, chunksize=10):
            done += 1; st_["status: " + status] += 1
            if status == "ok":
                n = len(U[acc]["seq"]); st_["file length equals UniProt 2026_03"] += int(L == n); st_["file sequence equals UniProt 2026_03"] += int(seq_eq)
                if n > 1022: long_n += 1; long_ok += int(L == n)
                ai = aidx[acc]
                for p, vec in cells.items(): ES[rowof[int(offs[ai]) + p - 1]] = vec
            if done % 4000 == 0: log(f"  ESM1b files read {done} of {len(jobs)} | {time.time() - t0:.0f} s")
    log(f"ESM1b extraction {time.time() - t0:.0f} s")
    for k in sorted(st_): log(f"  {k}: {st_[k]}")
    log(f"  long proteins read {long_n}, with the right length {long_ok}")
    AMv = read_am(log, U, accs, aidx, offs, keys)
    hum_gk = np.array(gkH, np.int64); hum_alt = np.array([h[2] for h in hum], np.int64); hum_cls = np.array([h[3] for h in hum], np.int8)
    np.savez_compressed(CACHE, keys=keys, am=AMv, esm=ES, hum_gk=hum_gk, hum_alt=hum_alt, hum_cls=hum_cls)
    log(f"saved {CACHE} | md5 {H.md5_file(CACHE)}")
    log(f"requested cells with an AlphaMissense value (at least one substitution): {int(np.isfinite(AMv).any(axis=1).sum())} of {len(keys)} | ESM1b: {int(np.isfinite(ES).any(axis=1).sum())}")
    st = Store(); lines = []
    for s in SCORES:
        v, c = st.hum(s); ok = np.isfinite(v); a = v[ok & (c == 1)]; b = v[ok & (c == 0)]
        mg = MARGIN_FRAC * (float(np.median(a)) - float(np.median(b)))
        log(f"margin {s}: humsavar LP/P scored {len(a)} (median damage {np.median(a):.4f}), LB/B scored {len(b)} (median {np.median(b):.4f}) -> margin 0.10 x difference = {mg:.4f}")
        lines.append(f"{s}\t{mg:.6f}\t{np.median(a):.6f}\t{np.median(b):.6f}\t{len(a)}\t{len(b)}")
    with open(MARGINS, "w") as o: o.write("\n".join(lines) + "\n")
    info, VAL, PAIRS, ref = load_pairs(log, U, aidx, offs, st)
    for s in SCORES: log(f"coverage of all loss and control substitutions (core sets of SqnPAN_11): {s} {len(VAL[s])} of {len(info)} scored")
    log("wild-type check, structure check and margins done; no TCGA pair was scored in this section")
    write_ok(1, f"{os.path.basename(CACHE)}\t{H.md5_file(CACHE)}\n")

# ======================= Section 2: controls =======================
def section2(log):
    need_ok(1, "section 1"); C.check_seed_blocks(log, BLOCK_B, SEEDS_B.values())
    U = H.load_uniprot_main(); accs = sorted(U); aidx = {a: i for i, a in enumerate(accs)}; offs = offsets(U, accs)
    st = Store(); mg = read_margins(); info, VAL, PAIRS, ref = load_pairs(log, U, aidx, offs, st); res = {}; lines = []
    for s in SCORES:
        v, c = st.hum(s); ok = np.isfinite(v); a = v[ok & (c == 1)]; b = v[ok & (c == 0)]; auc = auroc(a, b)
        d1, a1, _l1, _c1, k1 = pair_arrays(PAIRS[1], VAL[s], info); cov = len(d1) / len(PAIRS[1])
        m6 = dict(PAIRS[6]); e = []; ea = []
        for k, c1 in PAIRS[1]:
            c2 = m6.get(k)
            if not c2 or k not in VAL[s]: continue
            v1 = [VAL[s][c] for c in c1 if c in VAL[s]]; v2 = [VAL[s][c] for c in c2 if c in VAL[s]]
            if v1 and v2: e.append(sum(v1) / len(v1) - sum(v2) / len(v2)); ea.append(info[k][0])
        pl, m = local_index(ea); bt = boot_sets([(np.array(e), pl)], m, [SEEDS_B["a12"], 99, SID[s]]); c95, c90 = ci(bt[:, 0]); med = float(np.median(e))
        ok_auc = auc >= AUC_MIN[s]; ok_draw = (not math.isnan(c95[0])) and c95[0] > -mg[s] and c95[1] < mg[s]; ok_cov = cov >= 0.95
        res[s] = ok_auc and ok_draw and ok_cov
        log(f"CONTROL {s}: humsavar AUROC {auc:.4f} (LP/P {len(a)}, LB/B {len(b)}; need >= {AUC_MIN[s]}) {'PASSED' if ok_auc else 'FAILED'}")
        log(f"CONTROL {s}: draw 1 vs draw 2, per-loss difference of the mean control damage: median {med:.4f}, 95% CI {c95[0]:.4f} to {c95[1]:.4f}, {len(e)} losses; margin {mg[s]:.4f} {'PASSED' if ok_draw else 'FAILED'}")
        log(f"CONTROL {s}: primary pairs scored {len(d1)} of {len(PAIRS[1])} = {cov:.3f} (need >= 0.95) {'PASSED' if ok_cov else 'FAILED'}")
        lines.append(f"{s} {'PASSED' if res[s] else 'FAILED'} auroc {auc:.4f} draw_median {med:.4f} draw_ci {c95[0]:.4f} {c95[1]:.4f} coverage {cov:.3f} margin {mg[s]:.4f}")
    lines.append(f"written {time.strftime('%Y-%m-%d %H:%M:%S')}")
    with open(VERDICT, "w") as o: o.write("\n".join(lines) + "\n")
    for x in lines: log("VERDICT " + x)
    write_ok(2)

# ======================= Section 3: (a) the SqnPAN_11 pairs =======================
def section3(log):
    need_ok(2, "section 2"); C.check_seed_blocks(log, BLOCK_B, SEEDS_B.values())
    U = H.load_uniprot_main(); accs = sorted(U); aidx = {a: i for i, a in enumerate(accs)}; offs = offsets(U, accs)
    st = Store(); mg = read_margins(); vd = read_verdict(); info, VAL, PAIRS, ref = load_pairs(log, U, aidx, offs, st)
    hdr = ["run", "unit", "description", "measure", "n_pairs_unit", "n_scored", "n_proteins", "median_d", "ci95_lo", "ci95_hi", "ci90_lo", "ci90_hi", "margin", "label",
           "median_loss", "median_controls", "share_loss_above_cut", "share_controls_above_cut", "cut", "n_unique_controls"]
    cuts = {"FoldX": FOLDX_CUT, "AM": AM_CUT, "ESM": ESM_CUT}; tabs = {}
    for run in (1, 2, 3):
        base = SEEDS_B["a12"] if run < 3 else SEEDS_B["a3"]; rows = []
        with open(P211_RUN(run), newline="", encoding="utf-8") as f: p211 = {int(r["unit"]): r for r in csv.DictReader(f, delimiter="\t")}
        for u, desc, cs, kinds, allm in UNITS:
            r = p211[u]
            rows.append([run, u, desc, "FoldX", r["n_pairs"], r["n_pairs"], r["n_proteins"], r["median_d"], r["ci95_lo"], r["ci95_hi"], r["ci90_lo"], r["ci90_hi"], FOLDX_MARGIN, r["label"],
                         r["median_ddg_losses"], r["median_ddg_controls"], r["losses_gt1"], r["controls_gt1"], FOLDX_CUT, r["n_unique_controls"]])
            for s in SCORES:
                d, ac, lv, cu, _k = pair_arrays(PAIRS[u], VAL[s], info); pl, m = local_index(ac); cvv = np.array(sorted(cu.values()))
                if m < MIN_PROT: rows.append([run, u, desc, s, len(PAIRS[u]), len(d), m] + [""] * 5 + [mg[s], NOTTEST] + [""] * 6); continue
                bt = boot_sets([(d, pl)], m, [base, u, SID[s]]); c95, c90 = ci(bt[:, 0]); lab = label_rule(c95, c90, mg[s]) if vd[s] else NOTINT
                rows.append([run, u, desc, s, len(PAIRS[u]), len(d), m, fnum(np.median(d)), fnum(c95[0]), fnum(c95[1]), fnum(c90[0]), fnum(c90[1]), fnum(mg[s], 4), lab,
                             fnum(np.median(lv)), fnum(np.median(cvv)), fnum(np.mean(lv > cuts[s])), fnum(np.mean(cvv > cuts[s])), cuts[s], len(cvv)])
        tabs[run] = rows; write_tsv(T_A(run), hdr, rows); log(f"run {run} (bootstrap base {base}): {T_A(run)} md5 {H.md5_file(T_A(run))} | content md5 {C.body_md5(T_A(run))}")
    same = C.body_md5(T_A(1)) == C.body_md5(T_A(2)); log(f"run 1 and run 2 identical apart from the run column: {'yes' if same else 'NO'}")
    res_rows = [a + [b[13], b[8], b[9]] for a, b in zip(tabs[1], tabs[3])]
    write_tsv(os.path.join(H.TABLES, "P2_18_A_results.tsv"), hdr + ["label_run3", "ci95_lo_run3", "ci95_hi_run3"], res_rows)
    log("=== (a) TCGA pairs: d = damage(loss) - mean damage(controls); FoldX in kcal/mol (> 0 = losses more destabilising), AlphaMissense pathogenicity, ESM1b -LLR; labels from run 1 ===")
    for a, b in zip(tabs[1], tabs[3]):
        if a[13] == NOTTEST: log(f"  unit {a[1]} {a[3]}: {NOTTEST}"); continue
        log(f"  unit {a[1]} {a[3]}: {a[5]} pairs scored in {a[6]} proteins; median d {a[7]} (95% CI {a[8]} to {a[9]}; 90% CI {a[10]} to {a[11]}; margin {a[12]}) -> {a[13]} | run 3: {b[13]}"
            f" | median loss {a[14]}, controls {a[15]}; share above cut {a[18]}: losses {a[16]}, controls {a[17]}")
    nd = sum(1 for a, b in zip(tabs[1], tabs[3]) if a[13] != b[13]); log(f"units whose label differs in run 3: {nd}")
    d_all = {}
    for s in MEASURES:
        d, ac, lv, cu, ks = pair_arrays(PAIRS[1], VAL[s], info); d_all[s] = dict(zip(ks, d))
    cor = []
    for x, y in (("FoldX", "AM"), ("FoldX", "ESM"), ("AM", "ESM")):
        ks = [k for k in d_all[x] if k in d_all[y]]; r = spearman(np.array([d_all[x][k] for k in ks]), np.array([d_all[y][k] for k in ks]))
        lx = [VAL[x][k] for k in ks]; ly = [VAL[y][k] for k in ks]; r2 = spearman(np.array(lx), np.array(ly))
        cor.append([x, y, len(ks), fnum(r, 3), fnum(r2, 3)]); log(f"  Spearman {x} vs {y} over {len(ks)} primary pairs: d {r:.3f}; loss value {r2:.3f} (descriptive)")
    write_tsv(os.path.join(H.TABLES, "P2_18_A_correlation.tsv"), ["measure_x", "measure_y", "n_pairs", "spearman_d", "spearman_loss_value"], cor)
    write_ok(3)

# ======================= Section 4: (b) proteome-wide =======================
def section4(log):
    need_ok(2, "section 2"); C.check_seed_blocks(log, BLOCK_B, SEEDS_B.values())
    U = H.load_uniprot_main(); accs = sorted(U); aidx = {a: i for i, a in enumerate(accs)}; offs = offsets(U, accs)
    st = Store(); mg = read_margins(); vd = read_verdict()
    Eg = C.load_npz(ENG14); Pz = C.load_npz(POS14); isN = Eg["P_res"] == 0
    prot = Eg["P_prot"][isN].astype(np.int64); posn = Pz["P_posn"][isN].astype(np.int64); role = Eg["P_role"][isN].astype(np.int64); usable = Eg["P_usable"][isN]
    gk = offs[prot] + posn - 1; idx_of = dict(zip(gk.tolist(), range(len(gk)))); pl_ = np.full(len(gk), np.nan); rs_ = np.full(len(gk), np.nan); offl = offs.tolist()
    with gzip.open(RESID, "rt", newline="", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t"); hd = next(rd); ia, ip, ix, il, ir = (hd.index(k) for k in ("accession", "position", "aa", "plddt", "rsa"))
        for row in rd:
            if row[ix] != "N": continue
            i = aidx.get(row[ia])
            if i is None: continue
            j = idx_of.get(offl[i] + int(row[ip]) - 1)
            if j is not None and row[il] != "" and row[ir] != "": pl_[j] = float(row[il]); rs_[j] = float(row[ir])
    base_ok = usable & np.isfinite(pl_) & (pl_ >= 70.0) & np.isfinite(rs_)
    log(f"N positions of the Q2 list {len(gk)} | usable protein with pLDDT and exposure and pLDDT >= 70: {int(base_ok.sum())}")
    for t, nm in ((0, "background"), (1, "T1"), (2, "T2"), (3, "T3")): log(f"  {nm} N positions: all {int((role == t).sum())}, kept {int(((role == t) & base_ok).sum())}")
    pbin = (pl_ >= 90).astype(np.int64); rbin = (rs_ >= 0.25).astype(np.int64)
    dam = {}
    for s in SCORES:
        r = st.rows(gk); D = np.full((len(gk), 20), np.nan, np.float32); okr = r >= 0; D[okr] = st.dam[s][r[okr]]
        cnt = np.isfinite(D).sum(axis=1); sm = np.where(np.isfinite(D), D, 0).sum(axis=1)
        m19 = np.where(cnt == 19, sm / np.maximum(cnt, 1), np.nan); D7 = D[:, SNV_N]; m7 = np.where(np.isfinite(D7).all(axis=1), D7.mean(axis=1), np.nan)
        dam[s] = {"all19": m19, "snv7": m7}
        log(f"{s}: N positions kept with all 19 substitutions scored {int((base_ok & np.isfinite(m19)).sum())} of {int(base_ok.sum())} | with the 7 reachable ones {int((base_ok & np.isfinite(m7)).sum())}")
    hdr = ["run", "score", "tier", "variant", "n_sequons", "n_proteins", "median_d", "ci95_lo", "ci95_hi", "ci90_lo", "ci90_hi", "margin", "label", "null_min", "null_median", "null_max",
           "null_share_abs_ge_observed", "n_background_positions"]
    tabs = {1: [], 2: [], 3: []}; nullrows = []
    for s in SCORES:
        for t in (1, 2, 3):
            for vi, var in enumerate(("all19", "snv7", "strict19")):
                x = dam[s]["snv7" if var == "snv7" else "all19"]
                pool = np.flatnonzero(base_ok & ((role == t) | (role == 0)) & np.isfinite(x))
                strat = prot[pool] if var != "strict19" else (prot[pool] * 4 + pbin[pool] * 2 + rbin[pool])
                uq, sti = np.unique(strat, return_inverse=True); sti = sti.ravel(); lab = role[pool] == t
                n1 = np.bincount(sti, weights=lab.astype(float), minlength=len(uq)).round().astype(np.int64); n0 = np.bincount(sti, minlength=len(uq)) - n1
                elig = (n1 > 0) & (n0 > 0); keep = elig[sti]; pool, sti, lab = pool[keep], sti[keep], lab[keep]
                uq2, sti = np.unique(sti, return_inverse=True); sti = sti.ravel(); n1 = np.bincount(sti, weights=lab.astype(float)).round().astype(np.int64); n0 = np.bincount(sti) - n1
                xv = x[pool].astype(float)
                def dvals(L):
                    s0 = np.bincount(sti, weights=np.where(L, 0.0, xv), minlength=len(n1)); mean0 = s0 / np.maximum(np.bincount(sti, weights=(~L).astype(float), minlength=len(n1)), 1)
                    return (xv - mean0[sti])[L]
                dlab = dvals(lab); pl, m = local_index([str(prot[p]) for p in pool[lab]])
                rng = np.random.default_rng(np.random.SeedSequence([SEEDS_B["null"], t, vi, SID[s]])); nul = np.empty(N_SHUF)
                for q in range(N_SHUF):
                    u_ = rng.random(len(sti)); order = np.lexsort((u_, sti)); sts = sti[order]; rank = np.arange(len(sts)) - np.searchsorted(sts, sts, side="left")
                    L2 = np.empty(len(sts), bool); L2[order] = rank < n1[sts]; nul[q] = float(np.median(dvals(L2)))
                obs = float(np.median(dlab)); shr = float(np.mean(np.abs(nul) >= abs(obs)))
                nullrows.append([s, t, var, fnum(obs), " ".join(fnum(v) for v in nul)])
                for run in (1, 2, 3):
                    base = SEEDS_B["b12"] if run < 3 else SEEDS_B["b3"]
                    if m < MIN_PROT: tabs[run].append([run, s, f"T{t}", var, len(dlab), m] + [""] * 6 + [NOTTEST, "", "", "", "", int((n0).sum())]); continue
                    bt = boot_sets([(dlab, pl)], m, [base, t, vi, SID[s]]); c95, c90 = ci(bt[:, 0]); lb = label_rule(c95, c90, mg[s]) if vd[s] else NOTINT
                    tabs[run].append([run, s, f"T{t}", var, len(dlab), m, fnum(obs), fnum(c95[0]), fnum(c95[1]), fnum(c90[0]), fnum(c90[1]), fnum(mg[s]), lb,
                                      fnum(nul.min()), fnum(np.median(nul)), fnum(nul.max()), fnum(shr, 3), int(n0.sum())])
    for run in (1, 2, 3): write_tsv(T_B(run), hdr, tabs[run]); log(f"run {run}: {T_B(run)} md5 {H.md5_file(T_B(run))} | content md5 {C.body_md5(T_B(run))}")
    log(f"run 1 and run 2 identical apart from the run column: {'yes' if C.body_md5(T_B(1)) == C.body_md5(T_B(2)) else 'NO'}")
    write_tsv(os.path.join(H.TABLES, "P2_18_B_null.tsv"), ["score", "tier", "variant", "observed_median_d", "null_medians_50"], nullrows)
    write_tsv(os.path.join(H.TABLES, "P2_18_B_results.tsv"), hdr + ["label_run3", "ci95_lo_run3", "ci95_hi_run3"], [a + [b[12], b[7], b[8]] for a, b in zip(tabs[1], tabs[3])])
    log("=== (b) proteome-wide: d = mean damage over the substitutions at a sequon N minus the mean at the background N of the same protein (positive = sequon N predicted less tolerant); labels from run 1 ===")
    for a, b in zip(tabs[1], tabs[3]):
        if a[12] == NOTTEST: log(f"  {a[1]} {a[2]} {a[3]}: {a[4]} sequons in {a[5]} proteins: {NOTTEST}"); continue
        log(f"  {a[1]} {a[2]} {a[3]}: {a[4]} sequons in {a[5]} proteins; median d {a[6]} (95% CI {a[7]} to {a[8]}; 90% CI {a[9]} to {a[10]}; margin {a[11]}) -> {a[12]} | run 3: {b[12]} | null (50 shuffles) {a[13]} to {a[15]}, share |null| >= |observed| {a[16]}")
    log("Confounders: conservation and burial; both scores already encode conservation, so a difference means predicted tolerance, not selection.")
    write_ok(4)

# ======================= Section 5: C1 aromatic sequon split =======================
def section5(log):
    need_ok(2, "section 2"); C.check_seed_blocks(log, BLOCK_C, SEEDS_C.values())
    U = H.load_uniprot_main(); accs = sorted(U); aidx = {a: i for i, a in enumerate(accs)}; offs = offsets(U, accs)
    st = Store(); mg = read_margins(); vd = read_verdict(); info, VAL, PAIRS, ref = load_pairs(log, U, aidx, offs, st)
    t1 = set()
    with open(SITES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["tier"] == "T1": t1.add((r["accession"], int(r["N_position"])))
    def ctx_of(k):
        acc, pos, wt, mt, rl, fm = info[k]; n = pos if rl == "N" else pos - 2
        if (acc, n) not in t1: return None
        s = U[acc]["seq"]
        def r_(off): j = n - 1 - off; return s[j] if j >= 0 else "-"
        p2 = s[n + 1] if n + 1 < len(s) else "-"
        return dict(m1=r_(1), m2=r_(2), m3=r_(3), p2=p2)
    FYW = "FYW"
    SPLITS = [("F/Y/W at N-2 (primary)", lambda c: c["m2"] in FYW), ("F only at N-2", lambda c: c["m2"] == "F"), ("F/Y/W at N-3", lambda c: c["m3"] in FYW),
              ("F/Y/W at N-2 or N-3", lambda c: c["m2"] in FYW or c["m3"] in FYW), ("F/Y/W at N-2 with T at +2", lambda c: c["m2"] in FYW and c["p2"] == "T"),
              ("F at N-1 (Price 2011)", lambda c: c["m1"] == "F")]
    hdr = ["run", "unit", "split", "measure", "n_group_pairs", "n_group_proteins", "median_group", "g_ci95_lo", "g_ci95_hi", "g_ci90_lo", "g_ci90_hi", "n_rest_pairs", "n_rest_proteins",
           "median_rest", "diff_group_minus_rest", "diff_ci95_lo", "diff_ci95_hi", "diff_ci90_lo", "diff_ci90_hi", "status"]
    tabs = {1: [], 2: [], 3: []}
    for u in (1, 5):
        cx = {}; unres = 0
        for k, _c in PAIRS[u]:
            c = ctx_of(k)
            if c is None: unres += 1
            else: cx[k] = c
        log(f"unit {u}: {len(PAIRS[u])} pairs | sequon context resolved {len(cx)} | no T1 sequon found for the pair {unres} (left out)")
        splits = SPLITS if u == 1 else SPLITS[:1]
        for si, (sname, fn) in enumerate(splits):
            for s in MEASURES:
                d, ac, lv, cu, ks = pair_arrays(PAIRS[u], VAL[s], info); inn = np.array([k in cx for k in ks], bool)
                g = np.array([bool(fn(cx[k])) if k in cx else False for k in ks], bool); rest = inn & ~g
                pl, m = local_index([a for a, i_ in zip(ac, inn) if i_]); amap = {a: i for i, a in enumerate(sorted({a for a, i_ in zip(ac, inn) if i_}))}
                ac_a = np.array(ac, object); dg = d[g]; dr = d[rest]
                pg = np.array([amap[a] for a in ac_a[g]], np.int64); pr = np.array([amap[a] for a in ac_a[rest]], np.int64)
                ng = len(set(ac_a[g].tolist())); nr = len(set(ac_a[rest].tolist()))
                for run in (1, 2, 3):
                    base = SEEDS_C["c12"] if run < 3 else SEEDS_C["c3"]
                    if ng < MIN_PROT or len(dg) < MIN_PAIRS_C1 or nr < MIN_PROT or len(dr) < MIN_PAIRS_C1:
                        tabs[run].append([run, u, sname, s, len(dg), ng, "", "", "", "", "", len(dr), nr] + [""] * 6 + [NOTTEST]); continue
                    bt = boot_sets([(dg, pg), (dr, pr)], m, [base, u, si, SID[s]]); diff = bt[:, 0] - bt[:, 1]
                    c95g, c90g = ci(bt[:, 0]); c95d, c90d = ci(diff)
                    tabs[run].append([run, u, sname, s, len(dg), ng, fnum(np.median(dg)), fnum(c95g[0]), fnum(c95g[1]), fnum(c90g[0]), fnum(c90g[1]), len(dr), nr, fnum(np.median(dr)),
                                      fnum(np.median(dg) - np.median(dr)), fnum(c95d[0]), fnum(c95d[1]), fnum(c90d[0]), fnum(c90d[1]), "Tested" if (s == "FoldX" or vd[s]) else NOTINT])
    for run in (1, 2, 3): write_tsv(T_C(run), hdr, tabs[run]); log(f"run {run}: {T_C(run)} md5 {H.md5_file(T_C(run))} | content md5 {C.body_md5(T_C(run))}")
    log(f"run 1 and run 2 identical apart from the run column: {'yes' if C.body_md5(T_C(1)) == C.body_md5(T_C(2)) else 'NO'}")
    write_tsv(os.path.join(H.TABLES, "P2_18_C1_results.tsv"), hdr[1:] + ["diff_ci95_lo_run3", "diff_ci95_hi_run3"], [a[1:] + [b[15], b[16]] for a, b in zip(tabs[1], tabs[3])])
    log("=== C1 aromatic sequon split (d = damage(loss) - mean damage(controls); group = pairs whose sequon has the aromatic context; no equivalence labels) ===")
    for a, b in zip(tabs[1], tabs[3]):
        if a[19] == NOTTEST: log(f"  unit {a[1]} | {a[2]} | {a[3]}: group {a[4]} pairs in {a[5]} proteins: {NOTTEST}"); continue
        log(f"  unit {a[1]} | {a[2]} | {a[3]}: group {a[4]} pairs ({a[5]} proteins) median d {a[6]} (95% CI {a[7]} to {a[8]}); rest {a[11]} pairs median {a[13]}; difference {a[14]} (95% CI {a[15]} to {a[16]}; 90% CI {a[17]} to {a[18]}; run 3: {b[15]} to {b[16]}) [{a[19]}]")
    log("Limit for the text: FoldX does not model the glycan, so C1 cannot test the aromatic-GlcNAc stabilisation of Culyba 2011; Y and W are an extension of the papers, which tested Phe.")
    write_ok(5)

SECTIONS = {"1": ("structure_scores_margins", section1), "2": ("controls", section2), "3": ("A_tcga_pairs", section3), "4": ("B_proteome_wide", section4), "5": ("C1_aromatic", section5)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    C.run_chain("P2_18", SECTIONS, "1" if which == "all" else which)
# ---- END OF SCRIPT ----
