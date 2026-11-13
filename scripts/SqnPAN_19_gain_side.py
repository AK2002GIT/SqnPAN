# SqnPAN_19_gain_side.py - B3 gain side and shortlist (Amendment 4; added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (fixed 2026-11-07). Post hoc: no primary label changes.
# (a) Question: do tumours carry more or fewer substitutions that CREATE an N-glycosylation sequon than expected from trinucleotide context, compared with the same
#     substitution elsewhere in the same protein that creates no sequon?
#   Entries (one row per protein position x gain type, proteins with a SqnPAN_14 transcript; Engine slot N = gain type "new N", slot ST = gain type "new S/T"):
#     new N  : residue at p is not N, target N; creating (gain rule 1) if p+1 exists and is not P and p+2 is S or T.
#     new S/T: residue at p is not S or T, target S or T; creating (gain rule 2) if p-2 is N and p-1 is not P.
#     Left out of both types (counted): positions 1 to 3 of every existing sequon (N-X-S/T, X not P; so no loss-and-gain) and the P of every N-P-S/T (gain rule 3).
#     Opportunities = single-nucleotide changes of the GENCODE v36 codon (SqnPAN_14 transcript map) that make the target; flanks from the transcript; SqnPAN_14 channels.
#     Entries without a usable change are not stored (counted).
#   Plausibility (analysis plan Amendment 2, SqnPAN_06 helpers) of a creating entry is taken at its new sequon N (p for new N, p-2 for new S/T).
#   Roles (Engine tiers): T1 = creating and plausible = DECISION UNIT (main, T1, both, pooled); T2 = creating, plausible or unknown (plausible rows duplicated with their events);
#     T3 = creating and implausible (comparison only, not judged: is a shift specific to sites that can be glycosylated?); 0 = background (not creating, any region).
#   Exposure = expected rate from R_main / R_nohyper of P2_14_engine.npz (A1, unchanged); Engine of SqnPAN_14_context_helpers.py unchanged; strata protein x gain type
#     (strict: also pLDDT >= 70 and RSA >= 0.25 bins at p from P2_07_residues, proteins with status use, as SqnPAN_08). Set-ups main, strict, nohyper, union, norescue,
#     uniform_restricted (entry counts); pooled units only (no per-cancer-type units).
#   Events: P2_06_mutations rows with alt N at a new-N entry or alt S/T (ref not S/T) at a new-S/T entry; one count per case x entry; chosen, norescue, union and
#     hypermutated flags as SqnPAN_08. GUARD (must be 0): SqnPAN_06 classify_change gains exactly the entry's N (creating) or changes nothing (background); SqnPAN_06 class and gain state agree.
#   Controls: calibration as A1 (background entries split at the stratum median of the pooled expected rate; judged: 90% CI inside +-0.32 log2, else every label gets the suffix
#     _calibration_FAILED); C2 natural null and C3 planted change as SqnPAN_08 / SqnPAN_14 (capability, Equivalent allowed). Labels, BH families, final-label rule 1A unchanged.
#     No suffix _biology_control_FAILED (the A2 rule belongs to the loss test; the text names the A2 result beside the gain result).
#   Reading: Enriched = more sequon-creating substitutions than expected (positive selection or a bias the model misses); Depleted = fewer (selection against new sequons);
#     Equivalent = a small shift at most. Not blind: Q1 gain counts (SqnPAN_09), FoldX gain ddG (SqnPAN_11) and the A3 gain signature shares were seen; no gain rate ratio before.
# (b) Shortlist (section 3, descriptive, hypotheses only): SqnPAN_06 rows of class gain, chosen file, gain_state_best plausible, one gained N; SURFY label surface
#     (Bausch-Fluck et al. 2018, data\raw\table_S3_surfaceome.xlsx); AlphaFold status use with one fragment; pLDDT >= 70 and RSA >= 0.25 at the new N; no P at N+3;
#     FoldX ddG (SqnPAN_11 gain table, set gain_plausible) <= 3 kcal/mol (rows without ddG left out); mutant length <= 4,000. FASTA of the full mutant sequences for
#     NetNGlyc 1.0 (DTU web server; Gupta and Brunak 2002), each file <= 190,000 aa and <= 2,000 sequences.
#   Section 6 (after the NetNGlyc runs, see README): positive = potential > 0.5 at the new N; ranked by potential, one substitution per protein, ties: more cases, higher RSA, then accession
#     and position; top 10 -> FoldX BuildModel (1 run, PDB kept) on the SqnPAN_10 repaired model; Re-Glyco sheet (glycan G54258NG = GS00086, Process Build, seeds 33290001 to
#     33290010; readout clash status and completeness, descriptive). Predicted occupancy and a clash-free glycan are not evidence of glycosylation.
# Seeds (block 33200000 to 33299999): C2/C3 runs 1 and 2 33200000, run 3 33220000; real runs 1 and 2 33240000, run 3 33260000; calibration 33280000; Re-Glyco 33290001 to 33290010.
# Usage (from the repository folder):  python scripts/SqnPAN_19_gain_side.py all [workers=4]   (sections 1 to 5; "4" restarts at the controls)
#                          python scripts/SqnPAN_19_gain_side.py 6               (after the NetNGlyc result pages are saved in data\raw\netnglyc)
# Memory: about 2 GB per worker (about 10 million entries); use workers=2 if Python stops with a memory error. Logs: results\logs\P2_19_section<k>_<name>.txt.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, re, csv, gzip, time, glob, html, shutil, array, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_06_helpers as P
import SqnPAN_08_selection_controls as S8
import SqnPAN_14_context_helpers as C

ARGS = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
WORKERS = int(ARGS.get("workers", "4"))
BLOCK = 33200000
SEEDS = dict(cc1=33200000, cc3=33220000, real1=33240000, real3=33260000, calib=33280000, reglyco=33290000)
CTRL_BASES = {1: SEEDS["cc1"], 2: SEEDS["cc1"], 3: SEEDS["cc3"]}; REAL_BASES = {1: SEEDS["real1"], 2: SEEDS["real1"], 3: SEEDS["real3"]}
CAL_LIM = 0.32; PLDDT_CUT = 70.0; RSA_CUT = 0.25; DDG_MAX = 3.0; NG_POT = 0.5; TOP_N = 10; COVER_MIN = 0.999
FA_MAX_AA = 190000; FA_MAX_SEQ = 2000; NG_MAX_LEN = 4000; GLYCAN = "G54258NG"; GLYCAN_GS = "GS00086"

ENG14 = os.path.join(H.PROC_OUT, "P2_14_engine.npz"); ENG14_MD5 = "11da2bc41ad5271d27d099b1dace5f3e"
TMAP14 = os.path.join(H.TABLES, "P2_14_transcript_map.tsv"); TMAP14_MD5 = "4f5e35bf97046ca2b119725fd057a4e1"
GAINS11 = os.path.join(H.TABLES, "P2_11_foldx_gains_ddg.tsv"); GAINS11_MD5 = "327d591cb05e3d590edb55cb65d78f65"
SURFY = os.path.join(H.RAW, "table_S3_surfaceome.xlsx"); SURFY_MD5 = "bf95ed25ffc933fbe0c0383ab19b8652"
MD5_1 = [(H.UNIPROT_MAIN, S8.MD5[H.UNIPROT_MAIN]), (S8.SITES, S8.MD5[S8.SITES]), (S8.RESID, S8.MD5[S8.RESID]), (H.USE_TSV, S8.MD5[H.USE_TSV]),
         (ENG14, ENG14_MD5), (TMAP14, TMAP14_MD5), (C.GENCODE, C.GENCODE_MD5)]
MD5_2 = [(ENG14, ENG14_MD5), (S8.MUTS, S8.MD5[S8.MUTS]), (S8.CASES, S8.MD5[S8.CASES])]
MD5_3 = [(GAINS11, GAINS11_MD5), (SURFY, SURFY_MD5), (S8.MUTS, S8.MD5[S8.MUTS]), (S8.RESID, S8.MD5[S8.RESID]), (H.USE_TSV, S8.MD5[H.USE_TSV])]

ENT_NPZ = os.path.join(H.PROC_OUT, "P2_19_entries.npz"); ENGINE_NPZ = os.path.join(H.PROC_OUT, "P2_19_engine.npz")
OK1 = os.path.join(H.PROC_OUT, "P2_19_section1_ok.txt"); OK2 = os.path.join(H.PROC_OUT, "P2_19_section2_ok.txt"); OK3 = os.path.join(H.PROC_OUT, "P2_19_section3_ok.txt")
PARTS = os.path.join(H.PROC_OUT, "P2_19_parts"); CAPAB = os.path.join(H.TABLES, "P2_19_unit_capability.tsv")
VERDICT = os.path.join(H.PROC_OUT, "P2_19_controls_verdict.txt"); CALIB = os.path.join(H.TABLES, "P2_19_calibration.tsv")
CALIB_VERDICT = os.path.join(H.PROC_OUT, "P2_19_calibration_verdict.txt"); RESULTS = os.path.join(H.TABLES, "P2_19_gain_results.tsv")
DECISION = os.path.join(H.PROC_OUT, "P2_19_decision.txt")
SHORT = os.path.join(H.TABLES, "P2_19_shortlist.tsv"); NG_DIR = os.path.join(H.RAW, "netnglyc"); NG_TABLE = os.path.join(H.TABLES, "P2_19_netnglyc_scores.tsv")
TOP = os.path.join(H.TABLES, "P2_19_reglyco_input.tsv"); MUT_DIR = os.path.join(H.PAPER2, "results", "reglyco_gain")
def CTRL_TABLE(k): return os.path.join(H.TABLES, f"P2_19_C2C3_run{k}.tsv")
def REAL_TABLE(k): return os.path.join(H.TABLES, f"P2_19_gain_run{k}.tsv")
def FASTA(k): return os.path.join(H.PROC_OUT, f"P2_19_netnglyc_input_{k}.fasta")

CONFIGS = [dict(name="main", kind="main", mask="main", mode="ctx", spec="main", cancer=False),
           dict(name="strict", kind="strict", mask="main", mode="ctx", spec="main", cancer=False),
           dict(name="nohyper", kind="main", mask="nohyper", mode="ctx", spec="nohyper", cancer=False),
           dict(name="union", kind="main", mask="union", mode="ctx", spec="main", cancer=False),
           dict(name="norescue", kind="main", mask="norescue", mode="ctx", spec="main", cancer=False),
           dict(name="uniform_restricted", kind="main", mask="main", mode="uniform", spec=None, cancer=False)]
CFG_NAMES = [c["name"] for c in CONFIGS]
ROLE_OF = {"plausible": 1, "unknown": 2, "implausible": 3}; STATE_NAME = {1: "plausible", 2: "unknown", 3: "implausible"}
GTYPE = {"N": "new N", "ST": "new S/T", "both": "new N + new S/T"}
TMEAN = {"T1": "creating, plausible", "T2": "creating, plausible or unknown", "T3": "creating, implausible (comparison)"}
EMPTY_TOPO = {"topo": [], "tm": [], "im": [], "sig": []}

def _chg_to(codon, targets):
    """Single-nucleotide changes of codon that make a residue in targets (codon residue not in targets, not a stop)."""
    aa = C.CODON[codon]; out = []
    if aa == "*" or aa in targets: return out
    for k in range(3):
        for alt in "ACGT":
            if alt == codon[k]: continue
            aa2 = C.CODON[codon[:k] + alt + codon[k + 1:]]
            if aa2 in targets: out.append((k, alt, aa2))
    return out
CHG_T = {0: {c: _chg_to(c, "N") for c in C.CODON}, 1: {c: _chg_to(c, "ST") for c in C.CODON}}
assert max(len(v) for d in CHG_T.values() for v in d.values()) <= 9

def need_ok(path, what):
    if not os.path.exists(path) or not open(path).read().startswith("OK"): raise RuntimeError(f"{what} has not finished ({os.path.basename(path)} missing)")

def check_md5(log, pairs):
    for path, rec in pairs:
        got = H.md5_file(path); same = got == rec
        log(f"md5 {os.path.basename(path)}: {got} ({'same as recorded' if same else 'DIFFERENT, recorded ' + rec})")
        if not same: raise RuntimeError("an input md5 differs from the recorded value: stop and decide")

def read_tmap(U, G):
    M = {}
    with open(TMAP14, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            vl = G[r["transcript"]]; vi = [x["ver"] for x in vl].index(r["version"]); M[r["accession"]] = (r["transcript"], vi)
    bad = sum(1 for a, (tid, vi) in M.items() if C.translate(G[tid][vi]["seq"][G[tid][vi]["cs"]:G[tid][vi]["ce"]]) != U[a]["seq"])
    return M, bad

def offsets(U, accs):
    offs = np.zeros(len(accs) + 1, np.int64)
    for i, a in enumerate(accs): offs[i + 1] = offs[i] + len(U[a]["seq"])
    return offs

def composition(aa_codes, mask):
    cnt = np.bincount(aa_codes[mask].astype(np.int64), minlength=26); tot = cnt.sum()
    if not tot: return "(none)"
    top = np.argsort(-cnt)[:8]
    return " ".join(f"{chr(65 + i)} {cnt[i] / tot:.3f}" for i in top if cnt[i])

# ======================= Section 1: entries, roles, opportunities, strict bins =======================
def section1(log):
    log("Pre-specification: header of this script (fixed 2026-11-07). Post hoc, added after the primary analysis.")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__}")
    check_md5(log, MD5_1); t0 = time.time()
    U = H.load_uniprot_main(); accs = sorted(U); D14 = C.load_npz(ENG14)
    if accs != [str(x) for x in D14["accs"]]: raise RuntimeError("accession order differs from SqnPAN_14")
    G, _, _ = C.load_gencode(); M, bad = read_tmap(U, G)
    log(f"SqnPAN_14 transcript map: {len(M)} proteins | translation differs from UniProt: {bad} (must be 0) | {time.time() - t0:.0f} s")
    if bad: raise RuntimeError("transcript map does not reproduce")
    topo = P.load_topology(H.load_uniprot_features())
    hasT1 = {r["accession"] for r in csv.DictReader(open(S8.SITES, newline="", encoding="utf-8"), delimiter="\t") if r["tier"] == "T1"}
    use = {r["accession"] for r in csv.DictReader(open(H.USE_TSV, newline="", encoding="utf-8"), delimiter="\t") if r["status"] == "use"}
    log(f"proteins with a T1 site: {len(hasT1)} | AlphaFold status use: {len(use)}")
    offs = offsets(U, accs); aidx = {a: i for i, a in enumerate(accs)}
    PB = np.full(int(offs[-1]), -1, np.int8); RB = np.full(int(offs[-1]), -1, np.int8); rc = collections.Counter(); t1 = time.time()
    with gzip.open(S8.RESID, "rt", newline="", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t"); hdr = next(rd); ix = {k: hdr.index(k) for k in ("accession", "position", "aa", "plddt", "rsa")}
        for row in rd:
            a = row[ix["accession"]]
            if a not in use or a not in M: continue
            i = aidx.get(a)
            if i is None: continue
            pos = int(row[ix["position"]]); s = U[a]["seq"]
            if not (1 <= pos <= len(s)) or s[pos - 1] != row[ix["aa"]]: rc["residue rows whose letter differs from UniProt (not used)"] += 1; continue
            pl = row[ix["plddt"]]; rs = row[ix["rsa"]]
            if pl == "" or rs == "": rc["residue rows with a blank pLDDT or RSA (not used)"] += 1; continue
            g = offs[i] + pos - 1; PB[g] = 1 if float(pl) >= PLDDT_CUT else 0; RB[g] = 1 if float(rs) >= RSA_CUT else 0; rc["residue rows used"] += 1
    for k in sorted(rc): log(f"  {k}: {rc[k]}")
    log(f"residue table read in {time.time() - t1:.0f} s")
    E_prot = array.array("i"); E_posn = array.array("i"); E_res = array.array("b"); E_cls = array.array("b"); E_aa = array.array("b"); OPP = array.array("h")
    c = collections.Counter(); st = collections.Counter(); t1 = time.time()
    for a_i, acc in enumerate(accs):
        if acc not in M: continue
        s = U[acc]["seq"]; n = len(s); d = G[M[acc][0]][M[acc][1]]; Tx = d["seq"]; cs = d["cs"]; tp = topo.get(acc, EMPTY_TOPO); t1p = acc in hasT1
        excl = set()
        for m in re.finditer(r"(?=N[^P][ST])", s): j = m.start() + 1; excl.update((j, j + 1, j + 2))
        nps = {m.start() + 2 for m in re.finditer(r"(?=NP[ST])", s)}
        for p in range(1, n + 1):
            if p in excl: c["positions left out: inside an existing sequon"] += 1; continue
            if p in nps: c["positions left out: P of an N-P-S/T (gain rule 3)"] += 1; continue
            i = cs + 3 * (p - 1); cod = Tx[i:i + 3]
            if C.CODON.get(cod) != s[p - 1]: c["codon does not code the UniProt residue (must be 0)"] += 1; continue
            for gt in (0, 1):
                chg = CHG_T[gt][cod]
                if not chg: continue
                row = [-1] * 9; j = 0
                for kk, alt, _a in chg:
                    idx = i + kk
                    if idx < 1 or idx + 1 >= len(Tx): c["changes without both flanking bases (left out)"] += 1; continue
                    ch = C.chan_tx(Tx[idx - 1:idx + 2], alt)
                    if ch < 0: c["changes with a non-ACGT flank (left out)"] += 1; continue
                    row[j] = ch; j += 1
                if not j: c["entries without a usable change (not stored)"] += 1; continue
                if gt == 0: creating = p + 2 <= n and s[p] != "P" and s[p + 1] in "ST"; npos = p
                else: creating = p >= 3 and s[p - 3] == "N" and s[p - 2] != "P"; npos = p - 2
                if creating:
                    state = P.plausibility(tp, npos, t1p); cl = ROLE_OF[state]; st[(gt, state)] += 1
                else: cl = 0
                E_prot.append(a_i); E_posn.append(p); E_res.append(gt); E_cls.append(cl); E_aa.append(ord(s[p - 1]) - 65); OPP.extend(row)
    log(f"entries built in {time.time() - t1:.0f} s")
    for k in sorted(c): log(f"  {k}: {c[k]}")
    if c["codon does not code the UniProt residue (must be 0)"]: raise RuntimeError("codon check failed")
    prot = np.frombuffer(E_prot, np.int32).copy(); posn = np.frombuffer(E_posn, np.int32).copy(); res = np.frombuffer(E_res, np.int8).copy()
    cls = np.frombuffer(E_cls, np.int8).copy(); aa = np.frombuffer(E_aa, np.int8).copy(); O = np.frombuffer(OPP, np.int16).copy().reshape(-1, 9)
    gidx = offs[prot.astype(np.int64)] + posn.astype(np.int64) - 1; pb = PB[gidx]; rb = RB[gidx]
    use_arr = np.array([a in use for a in accs], bool); inuse = use_arr[prot]; usable = inuse & (pb >= 0) & (rb >= 0)
    share = float(usable.sum()) / max(int(inuse.sum()), 1)
    log(f"strict bins: entries in proteins with status use {int(inuse.sum())} | with pLDDT and RSA {int(usable.sum())} | share {share:.5f} (gate >= {COVER_MIN})")
    if share < COVER_MIN: raise RuntimeError("residue table does not cover the entries of usable proteins: stop and decide")
    log(f"--- entries stored: {len(prot)} in {len(np.unique(prot))} proteins ---")
    for gt, gn in ((0, "new N"), (1, "new S/T")):
        for cl, nm in ((0, "background (not creating)"), (1, "creating, plausible"), (2, "creating, unknown"), (3, "creating, implausible")):
            mk = (res == gt) & (cls == cl); log(f"  {gn} | {nm}: entries {int(mk.sum())} | proteins {len(np.unique(prot[mk]))}")
    log("--- wild-type residue composition at p (top 8 shares; descriptive check that the classes are alike) ---")
    for gt, gn in ((0, "new N"), (1, "new S/T")):
        for cl, nm in ((0, "background"), (1, "creating plausible"), (3, "creating implausible")):
            log(f"  {gn} | {nm}: " + composition(aa, (res == gt) & (cls == cl)))
    np.savez_compressed(ENT_NPZ, accs=np.array(accs), E_prot=prot, E_posn=posn, E_res=res, E_cls=cls, E_aa=aa, OPP=O, P_usable=usable,
                        P_pb=np.where(usable, pb, -1).astype(np.int8), P_rb=np.where(usable, rb, -1).astype(np.int8))
    with open(OK1, "w") as o: o.write(f"OK\n{os.path.basename(ENT_NPZ)}\t{H.md5_file(ENT_NPZ)}\n")
    log(f"saved {ENT_NPZ} | md5 {H.md5_file(ENT_NPZ)} | {time.time() - t0:.0f} s")

# ======================= Section 2: events, engine inputs, calibration =======================
def section2(log):
    need_ok(OK1, "section 1"); check_md5(log, MD5_2); C.check_seed_blocks(log, BLOCK, SEEDS.values())
    U = H.load_uniprot_main(); accs = sorted(U); Z = C.load_npz(ENT_NPZ); D14 = C.load_npz(ENG14)
    if accs != [str(x) for x in Z["accs"]]: raise RuntimeError("accession order differs from section 1")
    prot = Z["E_prot"].astype(np.int64); posn = Z["E_posn"].astype(np.int64); gt = Z["E_res"].astype(np.int64); cls = Z["E_cls"]; n0 = len(prot)
    offs = offsets(U, accs); emap = np.full(2 * int(offs[-1]), -1, np.int64); emap[(offs[prot] + posn - 1) * 2 + gt] = np.arange(n0)
    cases = {}
    with open(S8.CASES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"): cases[r["case_id"]] = (r["project"], r["hypermutated"])
    projects = [str(x) for x in D14["projects"]]
    if sorted({p for p, _ in cases.values()}) != projects: raise RuntimeError("projects of the case table differ from SqnPAN_14")
    pj = {p: i for i, p in enumerate(projects)}; aidx = {a: i for i, a in enumerate(accs)}
    ev = {}; g = collections.Counter(); g2 = collections.Counter(); t0 = time.time()
    with gzip.open(S8.MUTS, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ref = r["ref"]; alt = r["alt"]
            if alt == "N": t = 0
            elif alt in ("S", "T") and ref not in ("S", "T"): t = 1
            else: continue
            a = r["accession"]; i = aidx.get(a)
            if i is None: continue
            pos = int(r["position"]); s = U[a]["seq"]
            if s[pos - 1] != ref: g["reference residue differs from UniProt (must be 0)"] += 1; continue
            k = int(emap[(offs[i] + pos - 1) * 2 + t])
            if k < 0: g2["rows to N or S/T not at a stored entry (left-out position, protein without transcript or no usable change)"] += 1; continue
            lost, gained = P.classify_change(s, pos, ref, alt)
            want_g = ([pos] if t == 0 else [pos - 2]) if cls[k] > 0 else []
            if lost or gained != want_g: g["classify_change disagrees with the entry role (must be 0)"] += 1
            if r["class"] != ("gain" if cls[k] > 0 else "none"): g["P2_06 class differs from the entry role (must be 0)"] += 1
            if cls[k] > 0 and r["gain_state_best"] != STATE_NAME[int(cls[k])]: g["P2_06 gain state differs from the entry state (must be 0)"] += 1
            proj, hyp = cases[r["case_id"]]
            if r["project"] != proj or r["hypermutated"] != hyp: g["project or hypermutated flag differs from the case table (must be 0)"] += 1
            chosen = r["in_chosen_file"] == "1"; norescue = chosen and r["only_rescued"] == "0"
            key = (r["case_id"], k); e = ev.get(key)
            if e is None: ev[key] = [chosen, norescue]
            else: e[0] = e[0] or chosen; e[1] = e[1] or norescue; g2["second row of the same case x entry (counted once)"] += 1
    log(f"mutation table read in {time.time() - t0:.0f} s")
    log("--- mutation table cross-checks (must all be 0) ---")
    for k in ("reference residue differs from UniProt (must be 0)", "classify_change disagrees with the entry role (must be 0)", "P2_06 class differs from the entry role (must be 0)",
              "SqnPAN_06 gain state differs from the entry state (must be 0)", "project or hypermutated flag differs from the case table (must be 0)"): log(f"  {k}: {g[k]}")
    for k in sorted(g2): log(f"  {k}: {g2[k]}")
    if any(g.values()): raise RuntimeError("mutation table cross-check failed: stop and decide")
    keys = sorted(ev); ev_pos = np.array([k for _, k in keys], np.int64)
    ev_cancer = np.array([pj[cases[cid][0]] for cid, _ in keys], np.int8); ev_hyper = np.array([cases[cid][1] == "1" for cid, _ in keys], bool)
    ev_main = np.array([ev[k][0] for k in keys], bool); ev_norescue = np.array([ev[k][1] for k in keys], bool); ev_union = np.ones(len(keys), bool)
    log(f"counted case x entry events: union {len(keys)}, main {int(ev_main.sum())}, nohyper {int((ev_main & ~ev_hyper).sum())}, norescue {int(ev_norescue.sum())}")
    dup = np.flatnonzero(cls == 1); newidx = np.full(n0, -1, np.int64); newidx[dup] = n0 + np.arange(len(dup))
    rows = np.concatenate([np.arange(n0), dup]); role = np.concatenate([cls.astype(np.int8), np.full(len(dup), 2, np.int8)])
    dm = newidx[ev_pos] >= 0; ev_pos2 = np.concatenate([ev_pos, newidx[ev_pos[dm]]])
    def cat(x): return np.concatenate([x, x[dm]])
    E = dict(accs=Z["accs"], projects=D14["projects"], P_prot=Z["E_prot"][rows], P_res=Z["E_res"][rows], P_role=role, P_usable=Z["P_usable"][rows],
             P_pb=Z["P_pb"][rows], P_rb=Z["P_rb"][rows], P_ctx=np.ones(len(rows), bool), OPP=Z["OPP"][rows], ev_pos=ev_pos2.astype(np.int32),
             ev_cancer=cat(ev_cancer), ev_hyper=cat(ev_hyper), ev_main=cat(ev_main), ev_norescue=cat(ev_norescue), ev_union=cat(ev_union),
             R_main=D14["R_main"], R_nohyper=D14["R_nohyper"])
    log(f"engine rows {len(rows)} (entries {n0} + {len(dup)} plausible rows copied into T2) | events {len(ev_pos2)} (+{int(dm.sum())} copies)")
    log("--- coverage (main-model events) ---")
    evm = E["ev_pos"][E["ev_main"]].astype(np.int64)
    for t, nm in ((0, "background"), (1, "T1 creating plausible"), (2, "T2 creating plausible or unknown"), (3, "T3 creating implausible")):
        for rc_, rn in ((0, "new N"), (1, "new S/T")):
            mk = (E["P_role"] == t) & (E["P_res"] == rc_); log(f"  {nm} {rn}: entries {int(mk.sum())} | events {int(mk[evm].sum())}")
    np.savez_compressed(ENGINE_NPZ, **E); log(f"engine inputs {ENGINE_NPZ} | md5 {H.md5_file(ENGINE_NPZ)}")
    # ---- calibration control (copy of the A1 rule) ----
    Rm = E["R_main"]; Oa = E["OPP"].astype(np.int64); wp = np.where(Oa >= 0, Rm[np.maximum(Oa, 0), 0], 0.0).sum(axis=1)
    strat = E["P_prot"].astype(np.int64) * 2 + E["P_res"].astype(np.int64); idx = np.flatnonzero(E["P_role"] == 0)
    order = np.lexsort((wp[idx], strat[idx])); si = strat[idx][order]; wi = wp[idx][order]
    uq, start, cn = np.unique(si, return_index=True, return_counts=True); med = 0.5 * (wi[start + (cn - 1) // 2] + wi[start + cn // 2])
    high = wp[idx] > med[np.searchsorted(uq, strat[idx])]
    role_cal = np.full(len(rows), -1, np.int8); role_cal[idx] = np.where(high, 1, 0)
    log(f"calibration split: background entries {len(idx)} | high {int(high.sum())} | low {int((~high).sum())}")
    cfgs = [dict(name="calib_ctx", kind="main", mask="main", mode="ctx", spec="main", cancer=False), dict(name="calib_uniform", kind="main", mask="main", mode="uniform", spec=None, cancer=False)]
    eng = C.Engine(E, cfgs, restricted=True, tiers=(1,), role=role_cal)
    recs = eng.evaluate(eng.real_labels(), np.random.default_rng(np.random.SeedSequence(SEEDS["calib"])), None, C.LEVEL1, extras=True); del eng
    with open(CALIB, "w", newline="", encoding="utf-8") as o: o.write("\t".join(C.COLS_REAL) + "\n" + "".join(C.real_line(0, r) + "\n" for r in recs))
    log(f"calibration table {CALIB} | md5 {H.md5_file(CALIB)} (glyco = high-rate half, background = low-rate half)")
    for r in recs: log(f"  {r['cfg']} {GTYPE[r['group']]}: log2RR {r['est']:.3f} (90% CI {r['lo90']:.3f} to {r['hi90']:.3f}) events high {r['ev_g']} low {r['ev_b']} tested {r['tested']}")
    pr = [r for r in recs if r["cfg"] == "calib_ctx" and r["group"] == "both" and r["scope"] == "pooled"][0]
    ok = (not np.isnan(pr["lo90"])) and pr["lo90"] >= -CAL_LIM and pr["hi90"] <= CAL_LIM
    with open(CALIB_VERDICT, "w") as o: o.write(("PASSED" if ok else "FAILED") + f"\nmain, pooled, new N + new S/T, context exposure: log2RR {pr['est']:.4f}, 90% CI {pr['lo90']:.4f} to {pr['hi90']:.4f}, limit +-{CAL_LIM}\nwritten {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    log(f"CALIBRATION {'PASSED' if ok else 'FAILED (every B3 final label gets the suffix _calibration_FAILED)'}")
    with open(OK2, "w") as o: o.write(f"OK\n{os.path.basename(ENGINE_NPZ)}\t{H.md5_file(ENGINE_NPZ)}\n")

# ======================= Section 3: shortlist and NetNGlyc input (descriptive; no test result is read) =======================
def section3(log):
    need_ok(OK1, "section 1"); check_md5(log, MD5_3)
    U = H.load_uniprot_main()
    rows = H.read_xlsx(SURFY, "SurfaceomeMasterTable"); hi = next(i for i, r in enumerate(rows[:6]) if "UniProt accession" in r)
    hdr = rows[hi]; ia = hdr.index("UniProt accession"); il = hdr.index("Surfaceome Label"); isrc = hdr.index("Surfaceome Label Source")
    surf = {r[ia]: (r[isrc] if len(r) > isrc else "") for r in rows[hi + 1:] if len(r) > il and r[il] == "surface"}
    log(f"SURFY SurfaceomeMasterTable: header row {hi + 1} | label surface {len(surf)} accessions")
    use = {r["accession"]: r for r in csv.DictReader(open(H.USE_TSV, newline="", encoding="utf-8"), delimiter="\t")}
    t1pos = collections.defaultdict(list)
    for r in csv.DictReader(open(S8.SITES, newline="", encoding="utf-8"), delimiter="\t"):
        if r["tier"] == "T1": t1pos[r["accession"]].append(int(r["N_position"]))
    ddg = {}
    with open(GAINS11, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["set"] == "gain_plausible" and r["fragment"] == "1": ddg[(r["accession"], int(r["position"]), r["wt"], r["mt"])] = float(r["ddg_mean"])
    log(f"SqnPAN_11 plausible gain ddG rows on single models: {len(ddg)}")
    G = {}; c = collections.Counter()
    with gzip.open(S8.MUTS, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["class"] != "gain" or r["in_chosen_file"] != "1" or r["gain_state_best"] != "plausible": continue
            c["plausible gain rows of the chosen files"] += 1
            gs = [x for x in r["gained_sequon_N"].split(";") if x]
            if len(gs) != 1: c["rows with more than one gained sequon (left out)"] += 1; continue
            k = (r["accession"], int(r["position"]), r["ref"], r["alt"]); e = G.setdefault(k, dict(N=int(gs[0]), cases=set(), proj=set()))
            e["cases"].add(r["case_id"]); e["proj"].add(r["project"])
    for k in sorted(c): log(f"  {k}: {c[k]}")
    need = {(k[0], e["N"]) for k, e in G.items()}; RV = {}
    with gzip.open(S8.RESID, "rt", newline="", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t"); h2 = next(rd); ix = {k: h2.index(k) for k in ("accession", "position", "plddt", "rsa")}
        for row in rd:
            kk = (row[ix["accession"]], int(row[ix["position"]]))
            if kk in need: RV[kk] = (row[ix["plddt"]], row[ix["rsa"]])
    steps = collections.OrderedDict((x, 0) for x in ("unique plausible gain substitutions", "on a SURFY surface protein", "AlphaFold status use, one fragment",
                                                     "pLDDT >= 70 and RSA >= 0.25 at the new N", "no proline at N+3", "FoldX ddG present and <= 3 kcal/mol", "mutant length <= 4,000"))
    out = []
    for k in sorted(G):
        acc, pos, ref, alt = k; e = G[k]; Np = e["N"]; s = U[acc]["seq"]; steps["unique plausible gain substitutions"] += 1
        if acc not in surf: continue
        steps["on a SURFY surface protein"] += 1
        u = use.get(acc)
        if u is None or u["status"] != "use" or u["n_fragments"] != "1": continue
        steps["AlphaFold status use, one fragment"] += 1
        pl, rs = RV.get((acc, Np), ("", ""))
        if pl == "" or rs == "" or float(pl) < PLDDT_CUT or float(rs) < RSA_CUT: continue
        steps["pLDDT >= 70 and RSA >= 0.25 at the new N"] += 1
        if Np + 3 <= len(s) and s[Np + 2] == "P": continue
        steps["no proline at N+3"] += 1
        dg = ddg.get(k)
        if dg is None or dg > DDG_MAX: continue
        steps["FoldX ddG present and <= 3 kcal/mol"] += 1
        if len(s) > NG_MAX_LEN: continue
        steps["mutant length <= 4,000"] += 1
        m = s[:pos - 1] + alt + s[pos:]
        if not (m[Np - 1] == "N" and m[Np] != "P" and m[Np + 1] in "ST"): raise RuntimeError(f"no sequon at the new N in the mutant: {k}")
        near = min(abs(Np - t) for t in t1pos[acc]) if t1pos.get(acc) else ""
        out.append(dict(accession=acc, gene=U[acc]["gene"], position=pos, ref=ref, alt=alt, mutation=f"{ref}{pos}{alt}", gain_type="new N" if Np == pos else "new S/T",
                        new_sequon_N=Np, sequon_mutant=m[Np - 1:Np + 3], n_cases=len(e["cases"]), projects=";".join(sorted(e["proj"])), plddt_at_N=pl, rsa_at_N=rs,
                        foldx_ddg=f"{dg:.3f}", nearest_T1_distance=near, surfy_label_source=surf[acc], length=len(s), seq=m))
    log("--- shortlist filters (unique substitutions remaining after each step) ---")
    for k, v in steps.items(): log(f"  {k}: {v}")
    for i, r in enumerate(out): r["id"] = f"g{i + 1:04d}"
    cols = ["id", "accession", "gene", "position", "ref", "alt", "mutation", "gain_type", "new_sequon_N", "sequon_mutant", "n_cases", "projects", "plddt_at_N", "rsa_at_N",
            "foldx_ddg", "nearest_T1_distance", "surfy_label_source", "length"]
    with open(SHORT, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(cols); w.writerows([[r[c_] for c_ in cols] for r in out])
    log(f"shortlist {SHORT}: {len(out)} substitutions in {len({r['accession'] for r in out})} proteins | md5 {H.md5_file(SHORT)}")
    for old in glob.glob(FASTA("*")): os.remove(old)
    chunks = [[]]; aa_ = 0
    for r in out:
        if chunks[-1] and (aa_ + r["length"] > FA_MAX_AA or len(chunks[-1]) >= FA_MAX_SEQ): chunks.append([]); aa_ = 0
        chunks[-1].append(r); aa_ += r["length"]
    for k, ch in enumerate(chunks, 1):
        if not ch: continue
        with open(FASTA(k), "w", newline="\n", encoding="ascii") as o:
            for r in ch: o.write(f">{r['id']}\n" + "\n".join(r["seq"][i:i + 60] for i in range(0, len(r["seq"]), 60)) + "\n")
        log(f"NetNGlyc input {FASTA(k)}: {len(ch)} sequences ({ch[0]['id']} to {ch[-1]['id']}), {sum(r['length'] for r in ch)} aa | md5 {H.md5_file(FASTA(k))}")
    with open(OK3, "w") as o: o.write(f"OK\n{os.path.basename(SHORT)}\t{H.md5_file(SHORT)}\n")

# ======================= Section 4: controls C2, C3 =======================
def section4(log):
    need_ok(OK2, "section 2"); C.check_seed_blocks(log, BLOCK, SEEDS.values())
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__} | workers {WORKERS} | repeats {S8.N_REP} | bootstrap {S8.B}")
    C.run_controls(log, ENGINE_NPZ, CONFIGS, True, CTRL_BASES, PARTS, CTRL_TABLE, WORKERS)
    cap, c2ok, same12 = C.capability(log, CTRL_TABLE, CFG_NAMES, CAPAB)
    pk = ("main", "T1", "both", "pooled")
    lines = ["PASSED" if same12 else "FAILED", f"C2/C3 run1 and run2 identical: {'yes' if same12 else 'no'}",
             *[f"C2 run{k} {cfg}: {'PASSED' if c2ok[(k, cfg)] else 'FAILED'}" for k in (1, 3) for cfg in CFG_NAMES],
             f"decision unit (main, T1, both, pooled): {cap[pk][0]}; Equivalent allowed: {'yes' if cap[pk][1] else 'no'}", f"written {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    with open(VERDICT, "w") as o: o.write("\n".join(lines) + "\n")
    for x in lines: log("VERDICT " + x)

# ======================= Section 5: real result =======================
READING = {"Enriched": "more sequon-creating substitutions than expected from context (positive selection for new sequons, or a bias the model misses)",
           "Depleted": "fewer sequon-creating substitutions than expected from context (selection against new sequons)",
           "Equivalent": "sequon-creating substitutions occur at the expected rate (a small shift at most; not 'no selection')"}
def section5(log):
    if not os.path.exists(VERDICT) or open(VERDICT).readline().strip() != "PASSED": raise RuntimeError("controls verdict is not PASSED: no real result is computed")
    calib = open(CALIB_VERDICT).readline().strip(); suffix = "" if calib == "PASSED" else "_calibration_FAILED"
    log(f"controls verdict PASSED | calibration {calib} | label suffix '{suffix}'"); C.check_seed_blocks(log, BLOCK, SEEDS.values())
    eng = C.Engine(C.load_npz(ENGINE_NPZ), CONFIGS, restricted=True)
    C.real_runs(log, eng, REAL_BASES, REAL_TABLE); cap = C.read_cap(CAPAB)
    runs, fin = C.assemble(log, REAL_TABLE, cap, CFG_NAMES, RESULTS, suffix, ("calibration", "gain_type", "tier_meaning"), lambda key: [calib, GTYPE[key[2]], TMEAN[key[1]]])
    pk = ("main", "T1", "both", "pooled")
    log("=== DECISION UNIT (main, T1 creating plausible, new N + new S/T, pooled; context-aware exposure); added after the primary analysis ===")
    log("  " + C.unit_line(runs, fin, cap, pk))
    log("=== POOLED UNITS, all set-ups (T1 plausible, T2 plausible or unknown, T3 implausible = comparison) ===")
    for cfg in CFG_NAMES:
        for t in ("T1", "T2", "T3"):
            for g in ("N", "ST", "both"):
                k = (cfg, t, g, "pooled")
                if k in runs[1]: log(f"  {cfg} {t} {GTYPE[g]}: " + C.unit_line(runs, fin, cap, k))
    lab = fin[1][pk][0]; base = lab[:-len(suffix)] if suffix and lab.endswith(suffix) else lab
    t3 = runs[1][("main", "T3", "both", "pooled")]
    lines = [f"decision unit final label: {lab}" + (f" [{fin[1][pk][1]}]" if fin[1][pk][1] else ""), "READING: " + READING.get(base, "no conclusion (" + base + ")"),
             f"comparison T3 (implausible, cannot be glycosylated): log2RR {C.ff(t3['log2rr'])} (90% CI {C.ff(t3['lo90'])} to {C.ff(t3['hi90'])}), final {fin[1][('main', 'T3', 'both', 'pooled')][0]}",
             "No primary label changes. Calibration: " + calib + ".", f"written {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    with open(DECISION, "w") as o: o.write("\n".join(lines) + "\n")
    for x in lines: log("B3 " + x)

# ======================= Section 6: NetNGlyc scores, top 10, FoldX mutant models, Re-Glyco sheet =======================
NG_RE = re.compile(r"\b(g\d{4})\s+(\d+)\s+([A-Z]{4})\s+([0-9]*\.[0-9]+)\s+\((\d)/9\)\s+(\S+)")
def section6(log):
    need_ok(OK3, "section 3"); C.check_seed_blocks(log, BLOCK, SEEDS.values())
    import SqnPAN_10_foldx_helpers as FX
    if time.strftime("%Y-%m-%d") > FX.FOLDX_EXPIRY: raise RuntimeError(f"FoldX licence expired on {FX.FOLDX_EXPIRY}")
    files = sorted(f for f in glob.glob(os.path.join(NG_DIR, "*")) if os.path.isfile(f))
    if not files: raise RuntimeError(f"no NetNGlyc result files in {NG_DIR}")
    short = list(csv.DictReader(open(SHORT, newline="", encoding="utf-8"), delimiter="\t")); byid = {r["id"]: r for r in short}
    hits = {}; seen_ids = set(); dupl = 0
    for fp in files:
        txt = open(fp, encoding="utf-8", errors="replace").read(); txt = html.unescape(re.sub(r"<[^>]+>", " ", txt)); n = 0
        for m in NG_RE.finditer(txt):
            k = (m.group(1), int(m.group(2))); n += 1; seen_ids.add(m.group(1))
            v = (float(m.group(4)), int(m.group(5)), m.group(6), m.group(3), os.path.basename(fp))
            if k in hits and hits[k][:2] != v[:2]: raise RuntimeError(f"two different NetNGlyc rows for {k}")
            if k in hits: dupl += 1
            hits[k] = v
        log(f"{os.path.basename(fp)}: md5 {H.md5_file(fp)} | sequon rows {n}")
    log(f"sequon rows read {len(hits)} (repeated identical rows {dupl}) | sequences with rows {len(seen_ids)} of {len(short)}")
    found = 0; bad_motif = 0; miss = []
    for r in short:
        h = hits.get((r["id"], int(r["new_sequon_N"])))
        if h is None: miss.append(r["id"]); continue
        found += 1
        if h[3][:3] != r["sequon_mutant"][:3]: bad_motif += 1
        r["netnglyc_potential"] = f"{h[0]:.4f}"; r["jury"] = f"{h[1]}/9"; r["netnglyc_result"] = h[2]; r["netnglyc_file"] = h[4]
    log(f"shortlist substitutions with a NetNGlyc row at the new N: {found} of {len(short)} | motif differs (must be 0): {bad_motif}")
    if miss: log("  without a row: " + " ".join(miss[:60]) + (" ..." if len(miss) > 60 else ""))
    if bad_motif or found < 0.95 * len(short): raise RuntimeError("NetNGlyc output does not cover the shortlist (or sequons differ): check the saved pages")
    cols = list(short[0].keys())
    for x in ("netnglyc_potential", "jury", "netnglyc_result", "netnglyc_file"):
        if x not in cols: cols.append(x)
    with open(NG_TABLE, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(cols); w.writerows([[r.get(c_, "") for c_ in cols] for r in short])
    log(f"scores table {NG_TABLE} | md5 {H.md5_file(NG_TABLE)}")
    pos_ = [r for r in short if r.get("netnglyc_potential") and float(r["netnglyc_potential"]) > NG_POT]
    log(f"positive (potential > {NG_POT}): {len(pos_)} of {found} substitutions in {len({r['accession'] for r in pos_})} proteins")
    pos_.sort(key=lambda r: (-float(r["netnglyc_potential"]), -int(r["n_cases"]), -float(r["rsa_at_N"]), r["accession"], int(r["position"])))
    top = []; acc_seen = set()
    for r in pos_:
        if r["accession"] in acc_seen: continue
        acc_seen.add(r["accession"]); top.append(r)
        if len(top) == TOP_N: break
    os.makedirs(MUT_DIR, exist_ok=True); outrows = []
    for rank, r in enumerate(top, 1):
        acc = r["accession"]; pos = int(r["position"]); mut = r["mutation"]; tag = re.sub(r"[^A-Za-z0-9_-]", "_", f"{r['gene']}_{mut}")
        rep = os.path.join(FX.REPAIRED_DIR, f"{acc}_F1_Repair.pdb")
        if not os.path.exists(rep): raise RuntimeError(f"repaired model missing: {rep}")
        wd = os.path.join(FX.FOLDX_WORK, "P2_19_work", tag)
        if os.path.isdir(wd): shutil.rmtree(wd)
        t0 = time.time(); fx = FX.build(rep, [(r["ref"], "A", pos, r["alt"])], wd, runs=1, tag=f"BuildModel_{tag}", keep_pdbs=True)
        pdbs = sorted(p for p in glob.glob(os.path.join(wd, fx["base"] + "_1*.pdb")) if re.fullmatch(re.escape(fx["base"]) + r"_1(_0)?\.pdb", os.path.basename(p)))   # one run: FoldX writes <base>_1.pdb (seen 2026-11-07)
        if fx["rc"] != 0 or len(pdbs) != 1: raise RuntimeError(f"FoldX BuildModel failed for {tag}: rc {fx['rc']}, mutant PDBs {len(pdbs)} (see {fx['stdout']})")
        with open(pdbs[0], encoding="utf-8", errors="replace") as f: resl = FX.pdb_ca_residues(f.read())
        got = {(ch, num): FX.THREE.get(nm) for ch, num, ic, nm in resl}.get(("A", pos))
        if got != r["alt"]: raise RuntimeError(f"mutant PDB of {tag} has {got} at A{pos}, expected {r['alt']}")
        avg = FX.average_ddg(fx["average"]); odir = os.path.join(MUT_DIR, tag); os.makedirs(odir, exist_ok=True)
        dst = os.path.join(odir, f"{tag}_mutant.pdb"); shutil.copy2(pdbs[0], dst)
        log(f"  {rank}. {r['gene']} {mut} ({acc}) new sequon N{r['new_sequon_N']} {r['sequon_mutant']} | NetNGlyc {r['netnglyc_potential']} {r['jury']} {r['netnglyc_result']} | cases {r['n_cases']}"
            f" | pLDDT {r['plddt_at_N']} RSA {r['rsa_at_N']} | ddG P2_11 {r['foldx_ddg']}, this run {avg.get(1, ('', ''))[1]} | {time.time() - t0:.0f} s | {dst} md5 {H.md5_file(dst)}")
        outrows.append([rank, r["id"], r["gene"], acc, mut, r["new_sequon_N"], r["sequon_mutant"], r["netnglyc_potential"], r["jury"], r["n_cases"], r["projects"], r["plddt_at_N"],
                        r["rsa_at_N"], r["foldx_ddg"], os.path.relpath(dst, H.PAPER2), GLYCAN, GLYCAN_GS, SEEDS["reglyco"] + rank, os.path.relpath(odir, H.PAPER2), "", ""])
    hdr = ["rank", "id", "gene", "accession", "mutation", "glycosite_N", "sequon", "netnglyc_potential", "jury", "n_cases", "projects", "plddt_at_N", "rsa_at_N", "foldx_ddg",
           "mutant_pdb", "glycan_glytoucan", "glycoshape_id", "reglyco_seed", "reglyco_folder", "clash_status", "notes"]
    with open(TOP, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(hdr); w.writerows(outrows)
    log(f"Re-Glyco sheet {TOP}: {len(outrows)} rows | md5 {H.md5_file(TOP)} (hypotheses only: predicted occupancy and a clash-free glycan are not evidence of glycosylation)")

SECTIONS = {"1": ("entries", section1), "2": ("events_engine_calibration", section2), "3": ("shortlist_fasta", section3), "4": ("controls_C2_C3", section4), "5": ("real_result", section5)}
S6 = {"6": ("netnglyc_top10_mutants", section6)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which == "6": C.run_chain("P2_19", S6, "6")
    else: C.run_chain("P2_19", SECTIONS, "1" if which == "all" else which)
# ---- END OF SCRIPT ----
