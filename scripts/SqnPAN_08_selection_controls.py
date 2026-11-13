# SqnPAN_08_selection_controls.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-10-24; decisions 1A, 2A, 3A and 4A;
#   the analysis plan sections 6 and 10, with Amendments 1 and 2). This script runs the Q2 controls C2 and C3 only. It NEVER computes the rate ratio on the
#   real (unshuffled) glycosite labels: the real Q2 result is first computed in SqnPAN_09, after the verdict of this script.
# Question (Q2, primary test): are the residues of N-glycosites (the N, and the S/T at +2) hit by somatic missense mutations at a different rate from the same residue
#   type elsewhere in the same protein?
# Counting: chosen-file rows of all cases, rescued rows included (main data). One count per case, protein and position. At S/T only substitutions to a non-S/T residue
#   count (glycosite and background alike); at N every missense substitution counts.
# Position classes for tier t (rebuilt here from the UniProt 2026_03 sequence and P2_06_sites.tsv; the SqnPAN_06 column res_role is a cross-check only and must agree):
#   glycosite N = position 1 of a tier-t sequon; glycosite S/T = position 3 of a tier-t sequon; background = N, S or T of the same protein at no position of any sequon
#   of any tier (decision 1A: N, S or T at the middle position X of a sequon is left out of both classes, because X changed to P is itself a loss). Glycosite residues of
#   the other tiers are left out of both classes of the tier-t test.
# Effect: Mantel-Haenszel rate ratio, strata protein x residue type (N, S/T). RR = sum(a E0/E) / sum(b E1/E); a, b = counts at glycosite and background positions;
#   E1, E0 = cases x glycosite positions and cases x background positions; E = E1 + E0. The case count is the same in every stratum of a scope, so it cancels (the
#   weights are positions). Reported as log2 RR. A stratum enters only if it has at least one position of each class.
# Bootstrap: 2,000 resamples of the eligible proteins (at least one usable stratum of the residue group) with replacement; the strata of a protein move together; the same
#   resample is used for all scopes (pooled and 33 cancer types) of a tier x residue group. Resamples with a zero top sum are -inf, with a zero bottom sum +inf (kept and
#   counted); resamples with both zero are dropped and counted. Percentile CIs (order statistic k = ceil(q x n), no interpolation): 95% and 90%.
#   p-value (decision 2A): two-sided bootstrap p = min(1, 2 min[(1 + #{log2 RR* <= 0}) / (n + 1), (1 + #{log2 RR* >= 0}) / (n + 1)]), n = kept resamples.
# Tested unit (fixed before any data; label-independent, so the same in every shuffle): at least 20 informative proteins (eligible proteins with at least one counted
#   mutation in the usable strata of the scope) AND an expected glycosite count under the null of at least 10 (sum over counted mutations of glycosite positions / all
#   positions of their stratum). Otherwise "Not tested (too few data)". [Detail fixed while writing, before any data: the expected-count rule replaces "both sums > 0",
#   which would depend on the labels.] A tested unit with a zero observed top or bottom sum is "Inconclusive (zero count)"; a tested unit with more than 1% of resamples
#   dropped is "Inconclusive (unstable bootstrap)".
# Labels: Enriched / Depleted = BH FDR < 0.05 AND |log2 RR| >= 0.32 (direction from the sign); Equivalent = 90% CI inside [-0.14, 0.14]; Inconclusive = neither.
# Units (decision 3A, 144 units): main model pooled, tiers T1/T2/T3 x residue groups N, S/T, both (9); main model per cancer type, T1 x 3 groups x 33 types (99); and four
#   set-ups with the 9 pooled units each: strict (strata also split by pLDDT >= 70 or < 70 and relative exposure >= 0.25 or < 0.25 from P2_07_residues.tsv.gz, proteins with
#   AlphaFold status 'use' only), nohyper (cases with more than 1,000 missense rows in the chosen file left out), union (rows of a case's other primary-tumour files added,
#   each case-position counted once), norescue (only chosen-file rows with only_rescued = 0).
# BH families (within each set-up, level and repeat): the primary unit (pooled, T1, N + S/T together) alone; the other 8 pooled units; per cancer type one family of 33 for
#   each residue group.
# Control C2 (natural null): in each stratum (strict: also within each pLDDT x exposure bin) the tier-t label is moved at random to the same number of positions drawn from
#   that stratum's tier-t glycosite and background positions; each tier is shuffled on its own; then the full test. 20 repeats. PASSED for a set-up if Changed calls are
#   <= 5% of its tested units in every repeat. Reported, not judged: median log2 RR under the null and how often the 90% and 95% CIs contain 0.
# Control C3 (planted change, decision 4A): on the same shuffled data, each counted mutation at a pseudo-glycosite of tier t is kept with probability r (r < 1) or gets one
#   duplicate with probability r - 1 (r > 1); one uniform number per mutation and tier, shared by the levels. Levels: strong 1.5 and 2/3, weak 1.25 and 0.80; r = 1 is C2.
#   Per unit: "Test can answer" if each strong level is labelled Changed in the right direction in >= 16 of 20 repeats; otherwise the unit's real result can only be
#   "Inconclusive, test not sensitive". Equivalent allowed only if Equivalent in >= 16 of 20 repeats at r = 1 and in <= 1 of 20 at each strong level. Weak levels are
#   reported, not judged (1.25 lies on the Changed threshold, so about half or fewer are expected to be called Changed).
#   [Detail fixed while writing, before any data: a unit's capability needs the rule met in run1 AND run3; run2 is the identity check of run1.]
#   If C2 fails for a set-up, all its units are labelled "Inconclusive, test not sensitive (C2 failed)" and any amendment is dated.
# Seeds: base 31400000 for run1 and run2 (tables must be byte-identical), base 31500000 for run3; repeat seed = base + repeat number, split by NumPy SeedSequence into
#   streams for the main shuffle, the strict shuffle, the planting and the bootstrap (always drawn in the same order). All bases at least 100,000 apart from each other
#   and from 31000000 to 31300000 (checked in the code).
# Order: section 1 structure check and inputs (md5, position classes, res_role cross-check, saved arrays; no rate ratio); section 2 C2 and C3 (runs 1 to 3; a restart
#   continues from the saved repeats); section 3 verdict and the per-unit capability table that SqnPAN_09 reads.
# Usage (from the repository folder):  python scripts/SqnPAN_08_selection_controls.py 1
#                          python scripts/SqnPAN_08_selection_controls.py 2 [workers=N] [run=1|2|3]   (default: all three runs, 6 workers)
#                          python scripts/SqnPAN_08_selection_controls.py 3
#   Options for quick test runs only: limit=N (first N proteins), test=1 (3 repeats, 100 resamples).
# Logs: results\logs\P2_08_section<k>_<name>.txt. Outputs: data\processed\P2_08_inputs.npz, P2_08_inputs_ok.txt, P2_08_parts\ (one file per run and repeat),
#   P2_08_controls_verdict.txt; results\tables\P2_08_C2C3_run1.tsv to run3.tsv and P2_08_unit_capability.tsv.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")   # one BLAS thread per process: identical sums in every run
import sys, re, csv, gzip, math, time, collections
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import SqnPAN_05_helpers as H

ARGS = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
LIMIT = int(ARGS.get("limit", "0")) or None
TEST = ARGS.get("test", "0") == "1"
WORKERS = int(ARGS.get("workers", "6"))
RUN_ARG = ARGS.get("run", "all")

N_REP = 3 if TEST else 20
B = 100 if TEST else 2000
CHUNK = 250
SEED_USED = (31000000, 31100000, 31200000, 31300000)
RUN_BASE = {1: 31400000, 2: 31400000, 3: 31500000}
LEVELS = [("1.0", 1.0), ("1.5", 1.5), ("0.667", 2.0 / 3.0), ("1.25", 1.25), ("0.8", 0.8)]
STRONG = {"1.5": "Enriched", "0.667": "Depleted"}; WEAK = {"1.25": "Enriched", "0.8": "Depleted"}
CONFIGS = ["main", "strict", "nohyper", "union", "norescue"]
TIERS = [1, 2, 3]; GROUPS = ["N", "ST", "both"]
MIN_INF = 20; MIN_EXP = 10.0
FDR_CUT = 0.05; CHG = 0.32; EQ = 0.14; NAN_MAX = 0.01; C2_MAX = 0.05
PASS_N = math.ceil(0.8 * N_REP); EQ_STRONG_MAX = math.floor(0.05 * N_REP)
PLDDT_CUT = 70.0; RSA_CUT = 0.25

SITES = os.path.join(H.PROC, "P2_06_sites.tsv"); PROTS = os.path.join(H.PROC, "P2_06_proteins.tsv"); CASES = os.path.join(H.PROC, "P2_06_cases.tsv")
MUTS = os.path.join(H.PROC, "P2_06_mutations.tsv.gz"); RESID = os.path.join(H.PROC, "P2_07_residues.tsv.gz")
MD5 = {H.UNIPROT_MAIN: "af9af22f1e6eea80e96751c93c9da2b3", SITES: "34e0a114e0e7a6955418d837e9a5a245", PROTS: "f77f5ffabf9bdce910e4aa71f8099eb4",
       CASES: "1a131193724b29368975e7cc8d02c396", MUTS: "3be28f9265312957c13ab2f0faa2b2b2", RESID: "9681c55c6f177dd128ad9ab64fec12fd",
       H.USE_TSV: "bc009789438ac60310524d4ef3e6f901"}
NPZ = os.path.join(H.PROC_OUT, "P2_08_inputs.npz"); NPZ_OK = os.path.join(H.PROC_OUT, "P2_08_inputs_ok.txt")
PARTS = os.path.join(H.PROC_OUT, "P2_08_parts"); VERDICT = os.path.join(H.PROC_OUT, "P2_08_controls_verdict.txt")
CAPAB = os.path.join(H.TABLES, "P2_08_unit_capability.tsv")
def RUN_TABLE(k): return os.path.join(H.TABLES, f"P2_08_C2C3_run{k}.tsv")
COLS = ["repeat", "level", "config", "tier", "group", "scope", "family", "tested", "n_informative", "expected_glyco", "eligible_proteins", "num", "den",
        "log2rr", "lo95", "hi95", "lo90", "hi90", "p", "fdr", "n_inf", "n_dropped", "label"]

# ======================= Section 1: structure check and inputs =======================
def section1(log):
    log("Pre-specification: header of this script (fixed 2026-10-24, decisions 1A, 2A, 3A, 4A).")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__} | test={int(TEST)} limit={LIMIT}")
    for path, rec in MD5.items():
        got = H.md5_file(path); same = got == rec
        log(f"md5 {os.path.basename(path)}: {got} ({'same as recorded' if same else 'DIFFERENT, recorded ' + rec})")
        if not same: raise RuntimeError("an input md5 differs from the recorded value: stop and decide")
    U = H.load_uniprot_main()
    accs = sorted(U)
    if LIMIT: accs = accs[:LIMIT]; log(f"TEST RUN: first {LIMIT} proteins (sorted accessions) only")
    aidx = {a: i for i, a in enumerate(accs)}
    by_acc = collections.defaultdict(dict); n_tier = collections.Counter()
    with open(SITES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["accession"] in aidx:
                t = int(r["tier"][1]); by_acc[r["accession"]][int(r["N_position"])] = t; n_tier[t] += 1
    prot_no = {}
    with open(PROTS, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["accession"] in aidx: prot_no[r["accession"]] = (int(r["n_N_in_no_sequon"]), int(r["n_ST_in_no_sequon"]))
    log(f"proteins {len(accs)} | sequons by tier in the site table: T1 {n_tier[1]}, T2 {n_tier[2]}, T3 {n_tier[3]}")
    offs = np.zeros(len(accs) + 1, np.int64)
    for i, a in enumerate(accs): offs[i + 1] = offs[i] + len(U[a]["seq"])
    posmap = np.full(int(offs[-1]), -1, np.int32)
    P_prot, P_res, P_role, P_posn = [], [], [], []
    c = collections.Counter(); seq_mismatch = 0; noseq_mismatch = 0
    for i, a in enumerate(accs):
        s = U[a]["seq"]; sq = by_acc.get(a, {})
        if set(H.sequon_positions(s)[0]) != set(sq): seq_mismatch += 1
        role = {}
        for p, t in sq.items(): role[p] = t; role[p + 2] = t
        xpos = {p + 1 for p in sq}
        nb = [0, 0]
        for m in re.finditer("[NST]", s):
            j = m.start() + 1; res = 0 if m.group() == "N" else 1
            if j in role: code = role[j]
            elif j in xpos:
                c["excluded " + ("N" if res == 0 else "S/T") + " at the middle position X of a sequon (decision 1A)"] += 1; nb[res] += 1; continue
            else: code = 0; nb[res] += 1
            posmap[offs[i] + j - 1] = len(P_prot); P_prot.append(i); P_res.append(res); P_role.append(code); P_posn.append(j)
            c[("N" if res == 0 else "S/T") + (" background" if code == 0 else f" glycosite T{code}")] += 1
        if a in prot_no and tuple(nb) != prot_no[a]: noseq_mismatch += 1
    log("--- position classes (sequence facts, no mutation counts) ---")
    for k in sorted(c): log(f"  {k}: {c[k]}")
    log(f"proteins whose N-X-S/T positions differ from the site table: {seq_mismatch} (must be 0)")
    log(f"proteins whose background + excluded counts differ from n_N_in_no_sequon / n_ST_in_no_sequon of P2_06_proteins.tsv: {noseq_mismatch} (must be 0)")
    bad_tier = [t for t in TIERS if c[f"N glycosite T{t}"] != n_tier[t] or c[f"S/T glycosite T{t}"] != n_tier[t]]
    log(f"glycosite N and S/T positions equal the sequon count of each tier: {'yes' if not bad_tier else 'NO for tiers ' + str(bad_tier)}")
    if seq_mismatch or noseq_mismatch or bad_tier: raise RuntimeError("position classes do not agree with the SqnPAN_06 tables: stop and decide")
    P_prot = np.array(P_prot, np.int32); P_res = np.array(P_res, np.int8); P_role_a = np.array(P_role, np.int8); npos = len(P_prot)
    use = set()
    with open(H.USE_TSV, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["status"] == "use": use.add(r["accession"])
    P_usable = np.zeros(npos, bool); P_pb = np.full(npos, -1, np.int8); P_rb = np.full(npos, -1, np.int8)
    got = 0; aa_bad = 0; no_val = 0; not_use = 0
    with gzip.open(RESID, "rt", newline="", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t"); hdr = next(rd); ix = {k: hdr.index(k) for k in ("accession", "position", "aa", "plddt", "rsa")}
        for row in rd:
            aa = row[ix["aa"]]
            if aa not in ("N", "S", "T"): continue
            i = aidx.get(row[ix["accession"]])
            if i is None: continue
            pos = int(row[ix["position"]]); k = int(posmap[offs[i] + pos - 1])
            if k < 0: continue
            if U[accs[i]]["seq"][pos - 1] != aa: aa_bad += 1; continue
            pl = row[ix["plddt"]]; rs = row[ix["rsa"]]
            if pl == "" or rs == "": no_val += 1; continue
            if accs[i] not in use: not_use += 1; continue
            P_pb[k] = 1 if float(pl) >= PLDDT_CUT else 0; P_rb[k] = 1 if float(rs) >= RSA_CUT else 0; P_usable[k] = True; got += 1
    exp_use = int(sum(1 for k in range(npos) if accs[P_prot[k]] in use))
    log(f"strict model: positions of usable proteins {exp_use}, with pLDDT and exposure {got}; residue letter differs {aa_bad}; blank value {no_val}; rows of proteins not 'use' {not_use} (all three must be 0, and the first two numbers equal)")
    if got != exp_use or aa_bad or no_val or not_use: raise RuntimeError("residue table does not cover the usable positions: stop and decide")
    cases = {}
    with open(CASES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"): cases[r["case_id"]] = (r["project"], r["hypermutated"])
    projects = sorted({p for p, _ in cases.values()}); pj = {p: i for i, p in enumerate(projects)}
    if len(projects) != 33: raise RuntimeError(f"expected 33 projects, found {len(projects)}")
    ev = {}; g = collections.Counter(); g2 = collections.Counter()
    with gzip.open(MUTS, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ref = r["ref"]
            if ref not in ("N", "S", "T"): continue
            a = r["accession"]; i = aidx.get(a)
            if i is None: continue
            pos = int(r["position"])
            if U[a]["seq"][pos - 1] != ref: g["reference residue differs from UniProt"] += 1; continue
            k = int(posmap[offs[i] + pos - 1])
            if k < 0: want = "nonseqN" if ref == "N" else "nonseqST"
            else:
                code = P_role[k]
                want = ("nonseqN" if ref == "N" else "nonseqST") if code == 0 else (f"seqN_T{code}" if ref == "N" else f"seqST_T{code}")
            if r["res_role"] != want: g["res_role differs from the rebuilt class"] += 1
            proj, hyp = cases[r["case_id"]]
            if r["project"] != proj or r["hypermutated"] != hyp: g["project or hypermutated flag differs from the case table"] += 1
            if k < 0: g2["rows at excluded middle positions X (pooled over everything, not counted)"] += 1; continue
            if ref != "N" and r["alt"] in ("S", "T"): g2["S/T to S/T rows (not counted, rule of section 6)"] += 1; continue
            chosen = r["in_chosen_file"] == "1"; norescue = chosen and r["only_rescued"] == "0"
            key = (r["case_id"], k); e = ev.get(key)
            if e is None: ev[key] = [chosen, norescue]
            else: e[0] = e[0] or chosen; e[1] = e[1] or norescue
    log("--- mutation table cross-checks (must all be 0) ---")
    for k in ("reference residue differs from UniProt", "res_role differs from the rebuilt class", "project or hypermutated flag differs from the case table"): log(f"  {k}: {g[k]}")
    for k in sorted(g2): log(f"  {k}: {g2[k]}")
    if any(g.values()): raise RuntimeError("mutation table cross-check failed: stop and decide")
    keys = sorted(ev)
    ev_pos = np.array([k for _, k in keys], np.int32)
    ev_cancer = np.array([pj[cases[cid][0]] for cid, _ in keys], np.int8)
    ev_hyper = np.array([cases[cid][1] == "1" for cid, _ in keys], bool)
    ev_main = np.array([ev[k][0] for k in keys], bool); ev_norescue = np.array([ev[k][1] for k in keys], bool); ev_union = np.ones(len(keys), bool)
    log(f"counted case-position events, pooled over all position classes and tiers: union {len(keys)}, main {int(ev_main.sum())}, nohyper {int((ev_main & ~ev_hyper).sum())}, norescue {int(ev_norescue.sum())}")
    log(f"strict model usable positions: {int(P_usable.sum())} of {npos}")
    os.makedirs(H.PROC_OUT, exist_ok=True)
    tmp = NPZ + ".tmp.npz"
    np.savez_compressed(tmp, accs=np.array(accs), projects=np.array(projects), P_prot=P_prot, P_res=P_res, P_role=P_role_a, P_usable=P_usable, P_pb=P_pb, P_rb=P_rb,
                        ev_pos=ev_pos, ev_cancer=ev_cancer, ev_hyper=ev_hyper, ev_main=ev_main, ev_norescue=ev_norescue, ev_union=ev_union)
    os.replace(tmp, NPZ); md = H.md5_file(NPZ)
    with open(NPZ_OK, "w") as o: o.write(f"OK {md} test={int(TEST)} limit={LIMIT}\n")
    log(f"saved {NPZ} | md5 {md}")
    log("No rate ratio was computed in this section.")

# ======================= Section 2: controls C2 and C3 =======================
G = {}

def init_worker(npz_path):
    d = np.load(npz_path)
    prot = d["P_prot"].astype(np.int64); res = d["P_res"].astype(np.int64); role = d["P_role"]; usable = d["P_usable"]
    pb = d["P_pb"].astype(np.int64); rb = d["P_rb"].astype(np.int64)
    npos = len(prot); nprot = len(d["accs"]); ncan = len(d["projects"])
    ev_pos = d["ev_pos"].astype(np.int64); ev_cancer = d["ev_cancer"].astype(np.int64)
    G.update(projects=[str(x) for x in d["projects"]], nev=len(ev_pos), ncan=ncan)
    K = {}
    for kind in ("main", "strict"):
        if kind == "main": valid = np.ones(npos, bool); strat = prot * 2 + res
        else: valid = usable & (pb >= 0) & (rb >= 0); strat = ((prot * 2 + res) * 2 + np.maximum(pb, 0)) * 2 + np.maximum(rb, 0)
        for t in TIERS:
            pool = np.flatnonzero(valid & ((role == t) | (role == 0)))
            uniq, st = np.unique(strat[pool], return_inverse=True); st = st.ravel()
            isg = role[pool] == t
            n1 = np.bincount(st, weights=isg.astype(np.float64), minlength=len(uniq)).round().astype(np.int64)
            n0 = np.bincount(st, minlength=len(uniq)).astype(np.int64) - n1
            tot = (n1 + n0).astype(np.float64)
            if kind == "main": s_prot = uniq // 2; s_res = uniq % 2
            else: s_prot = uniq // 8; s_res = (uniq // 4) % 2
            loc = np.full(npos, -1, np.int64); loc[pool] = np.arange(len(pool))
            K[(kind, t)] = dict(st=st, n1=n1, w0=n0 / tot, w1=n1 / tot, elig=(n1 > 0) & (n0 > 0), s_prot=s_prot, s_res=s_res, loc=loc, npool=len(pool))
    EVM = {"main": d["ev_main"], "strict": d["ev_main"], "nohyper": d["ev_main"] & ~d["ev_hyper"], "union": d["ev_union"], "norescue": d["ev_norescue"]}
    S = {}
    for cfg in CONFIGS:
        kind = "strict" if cfg == "strict" else "main"
        for t in TIERS:
            Kt = K[(kind, t)]
            ei = np.flatnonzero(EVM[cfg]); j = Kt["loc"][ev_pos[ei]]; ok = j >= 0
            ei, j = ei[ok], j[ok]; s = Kt["st"][j]; ok = Kt["elig"][s]; ei, j, s = ei[ok], j[ok], s[ok]
            for g in GROUPS:
                if g == "both": sel_s = Kt["elig"]; ok = np.ones(len(ei), bool)
                else:
                    rc = 0 if g == "N" else 1; sel_s = Kt["elig"] & (Kt["s_res"] == rc); ok = Kt["s_res"][s] == rc
                eprot = np.unique(Kt["s_prot"][sel_s]); m = len(eprot)
                pmap = np.full(nprot, -1, np.int64); pmap[eprot] = np.arange(m)
                gi, gj, gs = ei[ok], j[ok], s[ok]; lp = pmap[Kt["s_prot"][gs]]; can = ev_cancer[gi]
                ns = 1 + ncan if (cfg == "main" and t == 1) else 1
                w1e = Kt["w1"][gs]
                inf = [len(np.unique(lp))]; exp = [float(w1e.sum())]
                if ns > 1:
                    for cc in range(ncan):
                        mk = can == cc; inf.append(len(np.unique(lp[mk]))); exp.append(float(w1e[mk].sum()))
                tested = [(a >= MIN_INF and b >= MIN_EXP) for a, b in zip(inf, exp)]
                S[(cfg, t, g)] = dict(ei=gi, j=gj, s=gs, lp=lp, can=can, m=m, ns=ns, inf=inf, exp=exp, tested=tested)
    G["K"] = K; G["S"] = S

def boot_product(rng, m, M):
    out = np.zeros((B, M.shape[1]))
    if m == 0: return out
    for a in range(0, B, CHUNK):
        b = min(B, a + CHUNK); n = b - a
        idx = rng.integers(0, m, size=(n, m)) + (np.arange(n, dtype=np.int64)[:, None] * m)
        W = np.bincount(idx.ravel(), minlength=n * m).reshape(n, m).astype(np.float64)
        out[a:b] = W @ M
    return out

def unit_stats(no, do, nb, db):
    inf = float("inf"); nan = float("nan")
    if no > 0 and do > 0: est = float(np.log2(no / do))
    elif no == 0 and do > 0: est = -inf
    elif do == 0 and no > 0: est = inf
    else: est = nan
    valid = ~((nb == 0) & (db == 0)); nv = int(valid.sum()); nbv = nb[valid]; dbv = db[valid]
    r = dict(est=est, num=float(no), den=float(do), n_dropped=B - nv)
    if nv == 0:
        r.update(lo95=nan, hi95=nan, lo90=nan, hi90=nan, p=nan, n_inf=0); return r
    with np.errstate(divide="ignore", invalid="ignore"): lb = np.log2(nbv / dbv)
    le = int((nbv <= dbv).sum()); ge = int((nbv >= dbv).sum())
    srt = np.sort(lb)
    def q(a, b):
        k = -(-a * nv // b); return float(srt[max(k - 1, 0)])
    r.update(lo95=q(25, 1000), hi95=q(975, 1000), lo90=q(50, 1000), hi90=q(950, 1000),
             p=min(1.0, 2.0 * min((1 + le) / (nv + 1), (1 + ge) / (nv + 1))), n_inf=int(np.isinf(lb).sum()))
    return r

def bh(p):
    p = np.asarray(p, float); n = len(p)
    if n == 0: return p
    o = np.argsort(p, kind="mergesort"); q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(q, 1.0); return out

def fmt(x, p=False):
    if isinstance(x, str): return x
    if x is None: return ""
    if isinstance(x, (bool, np.bool_)): return str(int(x))
    if isinstance(x, (int, np.integer)): return str(int(x))
    if math.isnan(x): return "nan"
    if math.isinf(x): return "inf" if x > 0 else "-inf"
    return f"{x:.6g}" if p else f"{x:.6f}"

def run_repeat(job):
    run, base, rep = job
    out = os.path.join(PARTS, f"run{run}_rep{rep:02d}.tsv")
    if os.path.exists(out): return rep, 0.0, True
    t0 = time.time()
    r_main, r_strict, r_plant, r_boot = [np.random.default_rng(x) for x in np.random.SeedSequence(base + rep).spawn(4)]
    lab = {}
    for kind, rng in (("main", r_main), ("strict", r_strict)):
        for t in TIERS:
            Kt = G["K"][(kind, t)]; st = Kt["st"]
            u = rng.random(Kt["npool"]); order = np.lexsort((u, st)); sts = st[order]
            start = np.searchsorted(sts, sts, side="left"); rank = np.arange(len(sts)) - start
            L = np.empty(len(sts), bool); L[order] = rank < Kt["n1"][sts]
            lab[(kind, t)] = L
    upl = r_plant.random((len(TIERS), G["nev"]))
    lines = []
    for cfg in CONFIGS:
        kind = "strict" if cfg == "strict" else "main"
        recs = []
        for t in TIERS:
            Kt = G["K"][(kind, t)]; L = lab[(kind, t)]
            for g in GROUPS:
                Sx = G["S"][(cfg, t, g)]; m = Sx["m"]; ns = Sx["ns"]
                isg = L[Sx["j"]]; u = upl[t - 1][Sx["ei"]]; w0 = Kt["w0"][Sx["s"]]; w1 = Kt["w1"][Sx["s"]]
                key0 = Sx["lp"] * ns; key1 = key0 + 1 + Sx["can"] if ns > 1 else None
                def agg(wt):
                    a = np.bincount(key0, weights=wt, minlength=m * ns)
                    if ns > 1: a = a + np.bincount(key1, weights=wt, minlength=m * ns)
                    return a.reshape(m, ns)
                cols = []
                for _name, rr in LEVELS:
                    if rr == 1.0: cnt = np.ones(len(u))
                    elif rr < 1.0: cnt = (u < rr).astype(np.float64)
                    else: cnt = 1.0 + (u < (rr - 1.0))
                    cols.append(agg(np.where(isg, cnt * w0, 0.0)))
                cols.append(agg(np.where(isg, 0.0, w1)))
                M = np.concatenate(cols, axis=1); BS = boot_product(r_boot, m, M); obs = M.sum(axis=0)
                for li, (name, _rr) in enumerate(LEVELS):
                    for sc in range(ns):
                        cn = li * ns + sc; cd = len(LEVELS) * ns + sc
                        st = unit_stats(obs[cn], obs[cd], BS[:, cn], BS[:, cd])
                        fam = "primary" if (t == 1 and g == "both" and sc == 0) else ("pooled_other" if sc == 0 else f"cancer_T1_{g}")
                        recs.append(dict(level=name, tier=t, group=g, scope="pooled" if sc == 0 else G["projects"][sc - 1], family=fam,
                                         tested=Sx["tested"][sc], n_inf_prot=Sx["inf"][sc], expected=Sx["exp"][sc], m=m, **st))
        for name, _rr in LEVELS:
            byfam = collections.defaultdict(list)
            for r in recs:
                if r["level"] != name: continue
                r["fdr"] = None
                if r["tested"] and r["num"] > 0 and r["den"] > 0 and not math.isnan(r["p"]): byfam[r["family"]].append(r)
            for fam, rl in byfam.items():
                for r, q in zip(rl, bh([r["p"] for r in rl])): r["fdr"] = float(q)
        for r in recs:
            if not r["tested"]: lab_ = "Not tested (too few data)"
            elif r["num"] == 0 or r["den"] == 0: lab_ = "Inconclusive (zero count)"
            elif r["n_dropped"] > NAN_MAX * B: lab_ = "Inconclusive (unstable bootstrap)"
            elif r["fdr"] is not None and r["fdr"] < FDR_CUT and abs(r["est"]) >= CHG: lab_ = "Enriched" if r["est"] > 0 else "Depleted"
            elif r["lo90"] >= -EQ and r["hi90"] <= EQ: lab_ = "Equivalent"
            else: lab_ = "Inconclusive"
            lines.append("\t".join(fmt(x) for x in [rep, r["level"], cfg, f"T{r['tier']}", r["group"], r["scope"], r["family"], bool(r["tested"]), r["n_inf_prot"]])
                         + "\t" + "\t".join([fmt(r["expected"]), fmt(r["m"]), fmt(r["num"]), fmt(r["den"]), fmt(r["est"]), fmt(r["lo95"]), fmt(r["hi95"]),
                                              fmt(r["lo90"]), fmt(r["hi90"]), fmt(r["p"], True), fmt(r["fdr"], True), fmt(r["n_inf"]), fmt(r["n_dropped"]), lab_]))
    tmp = out + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as o: o.write("\n".join(lines) + "\n")
    os.replace(tmp, out)
    return rep, time.time() - t0, False

def check_seeds(log):
    allb = list(SEED_USED) + sorted(set(RUN_BASE.values()))
    for i in range(len(allb)):
        for j in range(i + 1, len(allb)):
            if abs(allb[i] - allb[j]) < 100000: raise RuntimeError(f"seed bases {allb[i]} and {allb[j]} are closer than 100,000")
    if N_REP >= 100000: raise RuntimeError("repeat seeds would leave their base")
    log(f"seed bases {sorted(set(RUN_BASE.values()))} checked against {list(SEED_USED)}: all at least 100,000 apart; repeat seeds base + 1 to base + {N_REP}")

def c2_summary(log, run):
    rows = list(csv.DictReader(open(RUN_TABLE(run), newline="", encoding="utf-8"), delimiter="\t"))
    for cfg in CONFIGS:
        sh = []
        for rep in range(1, N_REP + 1):
            u = [r for r in rows if r["config"] == cfg and r["level"] == "1.0" and r["repeat"] == str(rep) and r["tested"] == "1"]
            ch = sum(1 for r in u if r["label"] in ("Enriched", "Depleted")); sh.append(f"{ch}/{len(u)}")
        log(f"  run{run} C2 {cfg}: Changed / tested units per repeat: {' '.join(sh)}")

def section2(log):
    if not os.path.exists(NPZ_OK) or not os.path.exists(NPZ): raise RuntimeError("section 1 has not finished (P2_08_inputs_ok.txt missing)")
    rec = open(NPZ_OK).read().split()[1]; md = H.md5_file(NPZ)
    log(f"inputs {NPZ} md5 {md} ({'same as section 1' if md == rec else 'DIFFERENT from section 1'})")
    if md != rec: raise RuntimeError("P2_08_inputs.npz changed after section 1")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__} | repeats {N_REP} | bootstrap resamples {B} | workers {WORKERS} | test={int(TEST)}")
    check_seeds(log)
    runs = [1, 2, 3] if RUN_ARG == "all" else [int(RUN_ARG)]
    os.makedirs(PARTS, exist_ok=True); os.makedirs(H.TABLES, exist_ok=True)
    pool = None
    if WORKERS > 1: pool = mp.get_context("spawn").Pool(WORKERS, initializer=init_worker, initargs=(NPZ,))
    else: init_worker(NPZ)
    try:
        for run in runs:
            base = RUN_BASE[run]; jobs = [(run, base, r) for r in range(1, N_REP + 1)]; t0 = time.time(); done = 0
            log(f"--- run{run}: base {base}, {N_REP} repeats ---")
            it = pool.imap_unordered(run_repeat, jobs) if pool else map(run_repeat, jobs)
            for rep, sec, skipped in it:
                done += 1
                log(f"  run{run} repeat {rep}: {'already saved, skipped' if skipped else f'{sec:.1f} s'} | {done}/{N_REP} done | elapsed {time.time() - t0:.0f} s")
                if done == 1 and not skipped:
                    log(f"  measured: first repeat {sec:.1f} s in one worker; projection for this run about {sec * N_REP / max(1, min(WORKERS, N_REP)) / 60:.1f} min (a guess)")
            body = []
            for r in range(1, N_REP + 1):
                with open(os.path.join(PARTS, f"run{run}_rep{r:02d}.tsv"), encoding="utf-8") as f: body.append(f.read())
            with open(RUN_TABLE(run), "w", newline="", encoding="utf-8") as o: o.write("\t".join(COLS) + "\n" + "".join(body))
            log(f"run{run} table {RUN_TABLE(run)} | md5 {H.md5_file(RUN_TABLE(run))} | seconds {time.time() - t0:.0f}")
            c2_summary(log, run)
    finally:
        if pool: pool.close(); pool.join()

# ======================= Section 3: verdict and capability table =======================
def section3(log):
    for k in (1, 2, 3):
        if not os.path.exists(RUN_TABLE(k)): raise RuntimeError(f"run{k} table missing: run section 2 first")
    m1, m2, m3 = (H.md5_file(RUN_TABLE(k)) for k in (1, 2, 3))
    same12 = m1 == m2
    log(f"run1 md5 {m1} | run2 md5 {m2} | identical: {'yes' if same12 else 'NO'} | run3 md5 {m3}")
    T = {}
    for k in (1, 3):
        T[k] = collections.defaultdict(list)
        for r in csv.DictReader(open(RUN_TABLE(k), newline="", encoding="utf-8"), delimiter="\t"):
            T[k][(r["config"], r["tier"], r["group"], r["scope"])].append(r)
    c2ok = {}
    for k in (1, 3):
        for cfg in CONFIGS:
            worst = 0.0; nfail = 0; fam = collections.Counter(); ntest = None
            for rep in range(1, N_REP + 1):
                u = [r for key, rl in T[k].items() if key[0] == cfg for r in rl if r["level"] == "1.0" and r["repeat"] == str(rep) and r["tested"] == "1"]
                ch = [r for r in u if r["label"] in ("Enriched", "Depleted")]
                for r in ch: fam[r["family"]] += 1
                share = len(ch) / len(u) if u else 0.0; worst = max(worst, share); ntest = len(u)
                if share > C2_MAX: nfail += 1
            ok = nfail == 0; c2ok[(k, cfg)] = ok
            log(f"C2 run{k} {cfg}: tested units {ntest} | largest Changed share {worst:.3f} (limit {C2_MAX}) | repeats over the limit {nfail} | Changed calls by family over all repeats {dict(fam)} | {'PASSED' if ok else 'FAILED'}")
    out = []; lab_count = collections.Counter(); primary = None
    for key in sorted(T[1]):
        cfg, tier, group, scope = key
        rl1 = T[1][key]; rl3 = T[3][key]
        tested = rl1[0]["tested"] == "1"
        if (rl3[0]["tested"] == "1") != tested: raise RuntimeError(f"tested status differs between run1 and run3 for {key}")
        rowv = [cfg, tier, group, scope, rl1[0]["family"], int(tested), rl1[0]["n_informative"], rl1[0]["expected_glyco"]]
        pw = []; eqa = []
        for k, rl in ((1, rl1), (3, rl3)):
            def n_lab(level, labs): return sum(1 for r in rl if r["level"] == level and r["label"] in labs)
            right = {lv: n_lab(lv, (d,)) for lv, d in list(STRONG.items()) + list(WEAK.items())}
            eq1 = n_lab("1.0", ("Equivalent",)); eqs = [n_lab(lv, ("Equivalent",)) for lv in STRONG]
            null = [float(r["log2rr"]) for r in rl if r["level"] == "1.0" and r["log2rr"] not in ("nan", "inf", "-inf")]
            cov90 = sum(1 for r in rl if r["level"] == "1.0" and r["lo90"] != "nan" and float(r["lo90"]) <= 0 <= float(r["hi90"]))
            cov95 = sum(1 for r in rl if r["level"] == "1.0" and r["lo95"] != "nan" and float(r["lo95"]) <= 0 <= float(r["hi95"]))
            med = float(np.median(null)) if null else float("nan")
            pw.append(all(right[lv] >= PASS_N for lv in STRONG)); eqa.append(eq1 >= PASS_N and all(x <= EQ_STRONG_MAX for x in eqs))
            rowv += [right["1.5"], right["0.667"], right["1.25"], right["0.8"], eq1, eqs[0], eqs[1], fmt(med), cov90, cov95]
        if not tested: label = "Not tested (too few data)"; eq_ok = False
        elif not (c2ok[(1, cfg)] and c2ok[(3, cfg)]): label = "Inconclusive, test not sensitive (C2 failed)"; eq_ok = False
        elif all(pw): label = "Test can answer"; eq_ok = all(eqa)
        else: label = "Inconclusive, test not sensitive"; eq_ok = False
        rowv += [label, "yes" if eq_ok else "no"]; out.append(rowv); lab_count[(cfg, "pooled" if scope == "pooled" else "cancer type", label, eq_ok)] += 1
        if key == ("main", "T1", "both", "pooled"): primary = (label, eq_ok, rowv)
    hdr = ["config", "tier", "group", "scope", "family", "tested", "n_informative", "expected_glyco"]
    for k in (1, 3): hdr += [f"run{k}_{x}" for x in ("changed_right_1.5", "changed_right_0.667", "changed_right_1.25", "changed_right_0.8", "equivalent_1.0",
                                                    "equivalent_1.5", "equivalent_0.667", "median_null_log2rr", "null_cover90", "null_cover95")]
    hdr += ["capability", "equivalent_allowed"]
    os.makedirs(H.TABLES, exist_ok=True)
    with open(CAPAB, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(hdr); w.writerows(out)
    log(f"capability table {CAPAB}: {len(out)} units | md5 {H.md5_file(CAPAB)}")
    log("--- units by set-up, scope, capability, Equivalent allowed ---")
    for k in sorted(lab_count): log(f"  {k}: {lab_count[k]}")
    log(f"--- pooled units (C3 counts of {N_REP}: right-direction Changed at 1.5, 0.667, 1.25, 0.8; Equivalent at 1.0, 1.5, 0.667; run1 then run3) ---")
    for rowv in out:
        if rowv[3] == "pooled": log("  " + " | ".join(str(x) for x in rowv[:4]) + " | run1 " + " ".join(str(x) for x in rowv[8:15]) + " | run3 " + " ".join(str(x) for x in rowv[18:25]) + f" | {rowv[-2]} | Eq {rowv[-1]}")
    c2_all = all(c2ok.values()); passed = same12 and c2_all
    lines = ["PASSED" if passed else "FAILED",
             f"run1 and run2 identical: {'yes' if same12 else 'no'}",
             *[f"C2 run{k} {cfg}: {'PASSED' if c2ok[(k, cfg)] else 'FAILED'}" for k in (1, 3) for cfg in CONFIGS],
             f"primary unit (main, T1, both, pooled): {primary[0] if primary else 'missing'}; Equivalent allowed: {'yes' if primary and primary[1] else 'no'}",
             f"written {time.strftime('%Y-%m-%d %H:%M:%S')} test={int(TEST)}"]
    with open(VERDICT, "w") as o: o.write("\n".join(lines) + "\n")
    for x in lines: log("VERDICT " + x)

SECTIONS = {"1": ("structure_inputs", section1), "2": ("controls_C2_C3", section2), "3": ("verdict", section3)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "1"
    for k in (SECTIONS if which == "all" else [which]):
        name, fn = SECTIONS[k]; H.run_section("P2_08", k, name, fn)
# ---- END OF SCRIPT ----
