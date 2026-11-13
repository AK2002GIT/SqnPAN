# SqnPAN_09_selection_test.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-10-24; decisions 1A and 2A;
#   the analysis plan sections 6, 10 and 11 with Amendments 1 and 2; SqnPAN_08_selection_controls.py header = the test definition).
# This script computes the REAL Q2 result (glycosite residues against the same residue type elsewhere in the same protein) with the test, units, thresholds and BH families
#   of SqnPAN_08, attaches the capability labels of SqnPAN_08 (results\tables\P2_08_unit_capability.tsv), and then builds the Q1 catalogue tables. Nothing about the test is changed:
#   position classes, Mantel-Haenszel rate ratio over protein x residue-type strata, 2,000 bootstrap resamples over proteins, p-values from the same resamples, labels
#   Enriched / Depleted (BH FDR < 0.05 and |log2 RR| >= 0.32), Equivalent (90% CI inside +-0.14), Inconclusive, tested-unit rule (>= 20 informative proteins and expected glycosite
#   count >= 10), the five set-ups (main, strict, nohyper, union, norescue) and the 144 units. Any change after the real result is read is a dated amendment, called post hoc.
# Code reuse: init_worker, boot_product, unit_stats, bh and fmt are imported from SqnPAN_08 (importing it does nothing). The labelling block of run_repeat is inline in SqnPAN_08 and cannot
#   be imported, so evaluate() below is an exact copy of it, generalised only in the list of planted levels. Control C-S (section 1) proves the copy: on label-shuffled data with the
#   SqnPAN_08 seeds (base 31400000, repeats 1 and 2) evaluate() must reproduce every line of P2_08_C2C3_run1.tsv for those repeats (all five levels, all 144 units, character by character).
#   If C-S fails, section 1 writes FAILED and no later section computes anything.
# Real labels: glycosite = position 1 or 3 of a tier-t sequon, background = N, S or T in no sequon and not at the middle position X (SqnPAN_08 decision 1A), rebuilt from
#   data\processed\P2_08_inputs.npz (md5 checked). Check: the real glycosite count of every stratum equals the count of SqnPAN_08 (all strata, all tiers, main and strict).
# Seeds (bootstrap only; the real labels are not random): base 31600000 for run 1 and run 2 (tables must be identical apart from the run column, compared by an md5 of the content without that column), base 31700000 for run 3. SeedSequence(base) is the only
#   stream. All bases at least 100,000 apart from each other and from 31000000 to 31500000 (checked in the code). Labels come from run 1 and are fixed now. Run 3 is reported
#   beside it (units whose raw or final label differs, largest difference in CI ends and FDR) and changes no label.
# Final labels (decision 1A; same rule as the companion LUAD study, final_label_F): the raw label of the unit stays in its own column.
#   - Enriched / Depleted stay as they are, also in a unit whose capability is not "Test can answer"; such a unit gets the flag "low-sensitivity unit" (C2 keeps false Changed calls
#     at <= 5%, so a Changed call is trusted; a low-sensitivity unit misses real effects and its estimates are inflated by selection: state it).
#   - Equivalent stays only where SqnPAN_08 allowed Equivalent; otherwise "Inconclusive, Equivalent not allowed".
#   - Inconclusive (also zero count, unstable bootstrap) in a unit whose capability is not "Test can answer": "Inconclusive, test not sensitive".
#   - Not tested units: "Not tested (too few data)".
# Set-ups are shown side by side (main model, strict, nohyper, union, norescue). The primary unit is main, T1, N + S/T together, pooled over cancer types.
# Q1 catalogue (decision 2A, descriptive only, no test, computed after the Q2 labels are written): chosen-file mutations (as main data of Q2, rescued rows included) and the same
#   without hypermutated cases. Tables: (a) counts by cancer type (33 + pooled) of mutations, cases and unique variants for loss of a sequon by tier (loss_T1, loss_T2, loss_T3, loss_any),
#   gain by plausibility state of Amendment 2 (gain_plausible, gain_unknown, gain_implausible, gain_plausible_or_unknown, gain_any) and all missense rows; a loss_and_gain mutation
#   counts in both a loss and a gain category; (b) recurrent variants (3 or more cancer types, or the top 30 by cases, in each category); (c) the position of the lost sequon hit
#   (N = position 1, X to P = position 2, S/T = position 3) by tier and cancer type. Guards: class counts and gain-state counts equal those of the SqnPAN_06 run (loss 11,735; gain 18,994;
#   loss_and_gain 187; none 1,429,876; plausible 3,654, unknown 13,630, implausible 1,897). Raw counts are not rates (they depend on protein length and mutation load); the Q2 test is the
#   selection test. The Q1 top variants and the position-3 excess are NOT blind (Fan et al. 2018); no comparison with Fan et al. is made here (SqnPAN_12).
# Order: section 1 inputs, final-label unit tests, control C-S, real-label check; section 2 the three real runs (no values printed); section 3 labels, results, set-ups side by side,
#   seed comparison; section 4 Q1 catalogue. "all" chains the four sections and stops at the first one that fails; section 2 and later refuse to run unless the C-S verdict says PASSED.
# Usage (from the repository folder):  python scripts/SqnPAN_09_selection_test.py all
#                          (single sections: 1, 2, 3 or 4 instead of all). Options for quick test runs only: test=1 (100 resamples, C-S not judged).
# Logs: results\logs\P2_09_section<k>_<name>.txt. Outputs: data\processed\P2_09_selfcheck_verdict.txt; results\tables\P2_09_Q2_run1.tsv to run3.tsv, P2_09_Q2_results.tsv,
#   P2_09_Q1_counts.tsv, P2_09_Q1_recurrent_variants.tsv, P2_09_Q1_loss_position_in_sequon.tsv.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")   # one BLAS thread: identical sums in every run
import sys, csv, gzip, math, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_08_selection_controls as S8

RUNS = {1: 31600000, 2: 31600000, 3: 31700000}
LEVEL1 = [("1.0", 1.0)]
REAL_TABLES = os.path.join(H.PAPER2, "results", "tables")         # P2_08 results are always read from the real folder
NPZ_REAL = os.path.join(H.PROC, "P2_08_inputs.npz")
P208 = {os.path.join(HERE, "SqnPAN_08_selection_controls.py"): "e871f68b6fa8e010f2d2a8465a06d1bf",
        NPZ_REAL: "c663941bae9e31dae5f68590a7624a4d",
        os.path.join(REAL_TABLES, "P2_08_unit_capability.tsv"): "de9fa49e84400ff5cf275e0d2220d887",
        os.path.join(REAL_TABLES, "P2_08_C2C3_run1.tsv"): "c872da79c028e3010c46d255b7ff540f",
        os.path.join(REAL_TABLES, "P2_08_C2C3_run2.tsv"): "c872da79c028e3010c46d255b7ff540f",
        os.path.join(REAL_TABLES, "P2_08_C2C3_run3.tsv"): "845a48d311da92c26fc0f7b5c8afecb1",
        os.path.join(H.PROC, "P2_08_controls_verdict.txt"): "4248da95cfb4e43c3323f273edbf792c"}
CAP_TSV = os.path.join(REAL_TABLES, "P2_08_unit_capability.tsv"); RUN1_P208 = os.path.join(REAL_TABLES, "P2_08_C2C3_run1.tsv")
VERD8 = os.path.join(H.PROC, "P2_08_controls_verdict.txt")
VERD9 = os.path.join(H.PROC_OUT, "P2_09_selfcheck_verdict.txt")
def RUN_REAL(k): return os.path.join(H.TABLES, f"P2_09_Q2_run{k}.tsv")
RESULTS = os.path.join(H.TABLES, "P2_09_Q2_results.tsv")
Q1_COUNTS = os.path.join(H.TABLES, "P2_09_Q1_counts.tsv"); Q1_VARS = os.path.join(H.TABLES, "P2_09_Q1_recurrent_variants.tsv")
Q1_POS = os.path.join(H.TABLES, "P2_09_Q1_loss_position_in_sequon.tsv")
COLS_REAL = ["run", "config", "tier", "group", "scope", "family", "tested", "n_informative", "expected_glyco", "eligible_proteins", "num", "den", "log2rr", "lo95", "hi95",
             "lo90", "hi90", "p", "fdr", "n_inf", "n_dropped", "raw_label", "n_events_glyco", "n_events_bg", "n_positions_glyco", "n_positions_bg"]
NOT_TESTED = "Not tested (too few data)"; CAN = "Test can answer"; NOTSENS = "Inconclusive, test not sensitive"; FLAG_LOW = "low-sensitivity unit"
GORDER = {"N": 0, "ST": 1, "both": 2}

# ======================= shared pieces =======================
def final_label(raw, cap, eq_ok):
    """Decision 1A. Returns (final label, flag)."""
    if raw == NOT_TESTED: return NOT_TESTED, ""
    sens = cap == CAN
    if raw in ("Enriched", "Depleted"): return raw, ("" if sens else FLAG_LOW)
    if raw == "Equivalent": return ("Equivalent", "") if eq_ok else ("Inconclusive, Equivalent not allowed", "")
    return (raw, "") if sens else (NOTSENS, "")

def check_inputs(log, with_muts=False):
    pairs = list(P208.items())
    if with_muts: pairs += [(S8.MUTS, S8.MD5[S8.MUTS]), (S8.CASES, S8.MD5[S8.CASES])]
    for path, rec in pairs:
        got = H.md5_file(path); same = got == rec
        log(f"md5 {os.path.basename(path)}: {got} ({'same as recorded' if same else 'DIFFERENT, recorded ' + rec})")
        if not same: raise RuntimeError("an input md5 differs from the recorded value: stop and decide")
    with open(VERD8) as f: first = f.readline().strip()
    log(f"SqnPAN_08 verdict first line: {first}")
    if first != "PASSED": raise RuntimeError("SqnPAN_08 controls verdict is not PASSED: no real result may be computed")

def check_seeds(log):
    allb = list(S8.SEED_USED) + sorted(set(S8.RUN_BASE.values())) + sorted(set(RUNS.values()))
    for i in range(len(allb)):
        for j in range(i + 1, len(allb)):
            if allb[i] != allb[j] and abs(allb[i] - allb[j]) < 100000: raise RuntimeError(f"seed bases {allb[i]} and {allb[j]} are closer than 100,000")
    log(f"seed bases of this script {sorted(set(RUNS.values()))} (run 1 and run 2 share the first) checked against {sorted(set(allb) - set(RUNS.values()))}: all at least 100,000 apart")

def require_verdict():
    if not os.path.exists(VERD9): raise RuntimeError("section 1 has not finished (P2_09_selfcheck_verdict.txt missing)")
    lines = open(VERD9).read().splitlines()
    ok = lines[0] == "PASSED" or (S8.TEST and lines[0] == "NOT JUDGED (test run)")
    if not ok: raise RuntimeError(f"control C-S verdict is '{lines[0]}': no real result may be computed")
    if not S8.TEST and "test=0" not in lines[-1]: raise RuntimeError("the C-S verdict comes from a test run")

def real_labels(d):
    role = d["P_role"]; usable = d["P_usable"]; pb = d["P_pb"].astype(np.int64); rb = d["P_rb"].astype(np.int64); npos = len(role); Lmap = {}
    for kind in ("main", "strict"):
        valid = np.ones(npos, bool) if kind == "main" else (usable & (pb >= 0) & (rb >= 0))
        for t in S8.TIERS:
            pool = np.flatnonzero(valid & ((role == t) | (role == 0)))
            Lmap[(kind, t)] = role[pool] == t
    return Lmap

def check_real_labels(log, Lmap):
    bad = 0
    for (kind, t), L in sorted(Lmap.items()):
        Kt = S8.G["K"][(kind, t)]
        if len(L) != Kt["npool"]: bad += 1; log(f"  {kind} T{t}: pool length {len(L)} differs from SqnPAN_08 {Kt['npool']}"); continue
        n1 = np.bincount(Kt["st"], weights=L.astype(np.float64), minlength=len(Kt["n1"])).round().astype(np.int64)
        same = bool(np.array_equal(n1, Kt["n1"])); bad += (not same)
        log(f"  {kind} T{t}: pool {len(L)} positions, glycosite positions {int(L.sum())}, strata {len(Kt['n1'])}; real glycosite count equals the SqnPAN_08 count in every stratum: {'yes' if same else 'NO'}")
    if bad: raise RuntimeError("real labels do not reproduce the glycosite counts of SqnPAN_08: stop and decide")

_TOT = {}
def stratum_totals(kind, t):
    if (kind, t) not in _TOT:
        Kt = S8.G["K"][(kind, t)]; _TOT[(kind, t)] = np.bincount(Kt["st"], minlength=len(Kt["n1"])).astype(np.int64)
    return _TOT[(kind, t)]

def evaluate(Lmap, r_boot, upl, levels, extras=False):
    """Exact copy of the labelling block of SqnPAN_08 run_repeat (levels: list of (name, rate ratio); upl = planting uniforms or None when only level 1.0 is used)."""
    G = S8.G; B = S8.B; out = []
    for cfg in S8.CONFIGS:
        kind = "strict" if cfg == "strict" else "main"
        recs = []
        for t in S8.TIERS:
            Kt = G["K"][(kind, t)]; L = Lmap[(kind, t)]
            for g in S8.GROUPS:
                Sx = G["S"][(cfg, t, g)]; m = Sx["m"]; ns = Sx["ns"]
                isg = L[Sx["j"]]; w0 = Kt["w0"][Sx["s"]]; w1 = Kt["w1"][Sx["s"]]
                u = upl[t - 1][Sx["ei"]] if upl is not None else None
                key0 = Sx["lp"] * ns; key1 = key0 + 1 + Sx["can"] if ns > 1 else None
                def agg(wt):
                    a = np.bincount(key0, weights=wt, minlength=m * ns)
                    if ns > 1: a = a + np.bincount(key1, weights=wt, minlength=m * ns)
                    return a.reshape(m, ns)
                cols = []
                for _name, rr in levels:
                    if rr == 1.0: cnt = np.ones(len(Sx["ei"]))
                    elif rr < 1.0: cnt = (u < rr).astype(np.float64)
                    else: cnt = 1.0 + (u < (rr - 1.0))
                    cols.append(agg(np.where(isg, cnt * w0, 0.0)))
                cols.append(agg(np.where(isg, 0.0, w1)))
                M = np.concatenate(cols, axis=1); BS = S8.boot_product(r_boot, m, M); obs = M.sum(axis=0)
                if extras:
                    if g == "both": sel = Kt["elig"]
                    else: sel = Kt["elig"] & (Kt["s_res"] == (0 if g == "N" else 1))
                    tot = stratum_totals(kind, t); pos_g = int(Kt["n1"][sel].sum()); pos_b = int((tot - Kt["n1"])[sel].sum())
                for li, (name, _rr) in enumerate(levels):
                    for sc in range(ns):
                        cn = li * ns + sc; cd = len(levels) * ns + sc
                        st = S8.unit_stats(obs[cn], obs[cd], BS[:, cn], BS[:, cd])
                        fam = "primary" if (t == 1 and g == "both" and sc == 0) else ("pooled_other" if sc == 0 else f"cancer_T1_{g}")
                        rec = dict(level=name, cfg=cfg, tier=t, group=g, scope="pooled" if sc == 0 else G["projects"][sc - 1], family=fam,
                                   tested=Sx["tested"][sc], n_inf_prot=Sx["inf"][sc], expected=Sx["exp"][sc], m=m, **st)
                        if extras:
                            mk = np.ones(len(isg), bool) if sc == 0 else (Sx["can"] == sc - 1)
                            rec.update(ev_g=int(isg[mk].sum()), ev_b=int((~isg[mk]).sum()), pos_g=pos_g, pos_b=pos_b)
                        recs.append(rec)
        for name, _rr in levels:
            byfam = collections.defaultdict(list)
            for r in recs:
                if r["level"] != name: continue
                r["fdr"] = None
                if r["tested"] and r["num"] > 0 and r["den"] > 0 and not math.isnan(r["p"]): byfam[r["family"]].append(r)
            for fam, rl in byfam.items():
                for r, q in zip(rl, S8.bh([r["p"] for r in rl])): r["fdr"] = float(q)
        for r in recs:
            if not r["tested"]: lab_ = NOT_TESTED
            elif r["num"] == 0 or r["den"] == 0: lab_ = "Inconclusive (zero count)"
            elif r["n_dropped"] > S8.NAN_MAX * B: lab_ = "Inconclusive (unstable bootstrap)"
            elif r["fdr"] is not None and r["fdr"] < S8.FDR_CUT and abs(r["est"]) >= S8.CHG: lab_ = "Enriched" if r["est"] > 0 else "Depleted"
            elif r["lo90"] >= -S8.EQ and r["hi90"] <= S8.EQ: lab_ = "Equivalent"
            else: lab_ = "Inconclusive"
            r["label"] = lab_
            out.append(r)
    return out

fmt = S8.fmt
def p208_line(rep, r):
    return ("\t".join(fmt(x) for x in [rep, r["level"], r["cfg"], f"T{r['tier']}", r["group"], r["scope"], r["family"], bool(r["tested"]), r["n_inf_prot"]])
            + "\t" + "\t".join([fmt(r["expected"]), fmt(r["m"]), fmt(r["num"]), fmt(r["den"]), fmt(r["est"]), fmt(r["lo95"]), fmt(r["hi95"]), fmt(r["lo90"]), fmt(r["hi90"]),
                                 fmt(r["p"], True), fmt(r["fdr"], True), fmt(r["n_inf"]), fmt(r["n_dropped"]), r["label"]]))

def real_line(run, r):
    return "\t".join([fmt(run), r["cfg"], f"T{r['tier']}", r["group"], r["scope"], r["family"], fmt(bool(r["tested"])), fmt(r["n_inf_prot"]), fmt(r["expected"]), fmt(r["m"]),
                      fmt(r["num"]), fmt(r["den"]), fmt(r["est"]), fmt(r["lo95"]), fmt(r["hi95"]), fmt(r["lo90"]), fmt(r["hi90"]), fmt(r["p"], True), fmt(r["fdr"], True),
                      fmt(r["n_inf"]), fmt(r["n_dropped"]), r["label"], fmt(r["ev_g"]), fmt(r["ev_b"]), fmt(r["pos_g"]), fmt(r["pos_b"])])

# ======================= Section 1: inputs, unit tests, control C-S, real-label check =======================
def section1(log):
    log("Pre-specification: header of this script (fixed 2026-10-24, decisions 1A and 2A).")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__} | test={int(S8.TEST)} | bootstrap resamples {S8.B}")
    check_inputs(log); check_seeds(log)
    cases = [("Enriched", CAN, True, ("Enriched", "")), ("Depleted", NOTSENS, False, ("Depleted", FLAG_LOW)), ("Enriched", NOTSENS, False, ("Enriched", FLAG_LOW)),
             ("Equivalent", CAN, True, ("Equivalent", "")), ("Equivalent", CAN, False, ("Inconclusive, Equivalent not allowed", "")),
             ("Equivalent", NOTSENS, False, ("Inconclusive, Equivalent not allowed", "")), ("Inconclusive", CAN, False, ("Inconclusive", "")),
             ("Inconclusive", NOTSENS, False, (NOTSENS, "")), ("Inconclusive (zero count)", CAN, False, ("Inconclusive (zero count)", "")),
             ("Inconclusive (unstable bootstrap)", NOTSENS, False, (NOTSENS, "")), (NOT_TESTED, NOT_TESTED, False, (NOT_TESTED, "")),
             ("Enriched", "Inconclusive, test not sensitive (C2 failed)", False, ("Enriched", FLAG_LOW))]
    bad = [c for c in cases if final_label(c[0], c[1], c[2]) != c[3]]
    log(f"final_label unit tests: {len(cases) - len(bad)} of {len(cases)} correct")
    if bad: raise RuntimeError(f"final_label unit test failed: {bad}")
    d = np.load(NPZ_REAL); S8.init_worker(NPZ_REAL); G = S8.G
    log(f"inputs loaded: {len(d['accs'])} proteins, {len(d['P_prot'])} classified positions, {G['nev']} case-position events (union), {G['ncan']} cancer types")
    Lmap = real_labels(d)
    log("--- real-label check (glycosite counts per stratum against SqnPAN_08) ---")
    check_real_labels(log, Lmap)
    log("--- control C-S: evaluate() on shuffled labels against P2_08_C2C3_run1.tsv ---")
    ref = collections.defaultdict(list)
    with open(RUN1_P208, encoding="utf-8") as f:
        f.readline()
        for line in f: ref[int(line.split("\t", 1)[0])].append(line.rstrip("\n"))
    ident = 0; total = 0; first_bad = None
    for rep in (1, 2):
        t0 = time.time(); base = S8.RUN_BASE[1]
        r_main, r_strict, r_plant, r_boot = [np.random.default_rng(x) for x in np.random.SeedSequence(base + rep).spawn(4)]
        lab = {}
        for kind, rng in (("main", r_main), ("strict", r_strict)):
            for t in S8.TIERS:
                Kt = G["K"][(kind, t)]; st = Kt["st"]
                u = rng.random(Kt["npool"]); order = np.lexsort((u, st)); sts = st[order]
                start = np.searchsorted(sts, sts, side="left"); rank = np.arange(len(sts)) - start
                L = np.empty(len(sts), bool); L[order] = rank < Kt["n1"][sts]
                lab[(kind, t)] = L
        upl = r_plant.random((len(S8.TIERS), G["nev"]))
        recs = evaluate(lab, r_boot, upl, S8.LEVELS)
        lines = [p208_line(rep, r) for r in recs]; want = ref[rep]
        same = sum(1 for a, b in zip(lines, want) if a == b); ident += same; total += len(want)
        if len(lines) != len(want): first_bad = first_bad or f"repeat {rep}: {len(lines)} lines against {len(want)}"
        for a, b in zip(lines, want):
            if a != b: first_bad = first_bad or f"repeat {rep}:\n   mine {a}\n   P2_08 {b}"; break
        log(f"  repeat {rep} (seed {base + rep}): {same} of {len(want)} lines identical to P2_08_C2C3_run1.tsv | {time.time() - t0:.0f} s")
    passed = ident == total and first_bad is None and total > 0
    if S8.TEST: verdict = "NOT JUDGED (test run)"; log("TEST RUN: 100 resamples, so the comparison with the 2,000-resample SqnPAN_08 table is not judged")
    else: verdict = "PASSED" if passed else "FAILED"
    if first_bad and not S8.TEST: log("first difference: " + first_bad)
    os.makedirs(H.PROC_OUT, exist_ok=True)
    with open(VERD9, "w") as o: o.write(verdict + f"\nlines identical {ident} of {total}\nwritten {time.strftime('%Y-%m-%d %H:%M:%S')} B={S8.B} test={int(S8.TEST)}\n")
    log(f"VERDICT C-S: {verdict} ({ident} of {total} lines identical)")
    if verdict == "FAILED": raise RuntimeError("control C-S failed: no real result was computed")
    log("No real rate ratio was computed in this section.")

# ======================= Section 2: the real runs =======================
def section2(log):
    require_verdict(); check_inputs(log); check_seeds(log)
    d = np.load(NPZ_REAL); S8.init_worker(NPZ_REAL)
    Lmap = real_labels(d); check_real_labels(log, Lmap)
    os.makedirs(H.TABLES, exist_ok=True)
    log(f"bootstrap resamples {S8.B} | test={int(S8.TEST)} | no values are printed in this section")
    for run in (1, 2, 3):
        base = RUNS[run]; t0 = time.time()
        recs = evaluate(Lmap, np.random.default_rng(np.random.SeedSequence(base)), None, LEVEL1, extras=True)
        with open(RUN_REAL(run), "w", newline="", encoding="utf-8") as o:
            o.write("\t".join(COLS_REAL) + "\n" + "".join(real_line(run, r) + "\n" for r in recs))
        log(f"run{run} (seed base {base}): {len(recs)} units | {time.time() - t0:.0f} s | {RUN_REAL(run)} | md5 {H.md5_file(RUN_REAL(run))} | content md5 without the run column {body_md5(run)}")

# ======================= Section 3: labels, results, set-ups side by side, seed comparison =======================
def pf(s): return float(s)
def ff(x):
    x = float(x)
    if math.isnan(x): return "nan"
    if math.isinf(x): return "inf" if x > 0 else "-inf"
    return f"{x:.3f}"
def fp(x):
    x = float(x)
    if math.isnan(x): return "nan"
    return f"{x:.2g}"
def ukey(k): return (S8.CONFIGS.index(k[0]), k[1], GORDER[k[2]], 0 if k[3] == "pooled" else 1, k[3])

def body_md5(k):
    """md5 of a real-run table without its first column (the run number), so that run 1 and run 2 can be compared."""
    import hashlib
    h = hashlib.md5()
    with open(RUN_REAL(k), encoding="utf-8", newline="") as f:
        for line in f: h.update(line.split("\t", 1)[1].encode("utf-8"))
    return h.hexdigest()

def read_run(k):
    with open(RUN_REAL(k), newline="", encoding="utf-8") as f:
        return {(r["config"], r["tier"], r["group"], r["scope"]): r for r in csv.DictReader(f, delimiter="\t")}

def section3(log):
    require_verdict()
    runs = {k: read_run(k) for k in (1, 2, 3)}
    m = {k: body_md5(k) for k in (1, 2, 3)}
    log(f"table content md5 without the run column: run1 {m[1]} | run2 {m[2]} | identical: {'yes' if m[1] == m[2] else 'NO'} | run3 {m[3]}")
    if m[1] != m[2]: log("WARNING: run1 and run2 differ although the seed is the same")
    with open(CAP_TSV, newline="", encoding="utf-8") as f: cap = {(r["config"], r["tier"], r["group"], r["scope"]): r for r in csv.DictReader(f, delimiter="\t")}
    if set(cap) != set(runs[1]) or set(cap) != set(runs[3]): raise RuntimeError("the unit lists of the real runs and of the SqnPAN_08 capability table differ")
    bad = [k for k in cap if runs[1][k]["tested"] != cap[k]["tested"] or runs[1][k]["n_informative"] != cap[k]["n_informative"]
           or runs[1][k]["expected_glyco"] != cap[k]["expected_glyco"] or runs[3][k]["tested"] != cap[k]["tested"]]
    log(f"units: {len(cap)} | tested flag, informative proteins and expected glycosite count equal the SqnPAN_08 values in every unit: {'yes' if not bad else 'NO (' + str(len(bad)) + ')'}")
    if bad: raise RuntimeError("tested units differ from SqnPAN_08: stop and decide")
    keys = sorted(cap, key=ukey); out = []; fin = {1: {}, 3: {}}
    for key in keys:
        c = cap[key]; eq_ok = c["equivalent_allowed"] == "yes"
        for k in (1, 3): fin[k][key] = final_label(runs[k][key]["raw_label"], c["capability"], eq_ok)
    hdr = ["config", "tier", "group", "scope", "family", "tested", "n_informative", "expected_glyco", "n_events_glyco", "n_events_bg", "n_positions_glyco", "n_positions_bg",
           "log2rr", "rate_ratio", "lo95", "hi95", "lo90", "hi90", "p", "fdr", "raw_label", "capability", "equivalent_allowed", "final_label", "flag",
           "run3_log2rr", "run3_lo95", "run3_hi95", "run3_fdr", "run3_raw_label", "run3_final_label"]
    for key in keys:
        r = runs[1][key]; r3 = runs[3][key]
        out.append([key[0], key[1], key[2], key[3], r["family"], r["tested"], r["n_informative"], r["expected_glyco"], r["n_events_glyco"], r["n_events_bg"], r["n_positions_glyco"],
                    r["n_positions_bg"], r["log2rr"], fmt(2.0 ** pf(r["log2rr"])), r["lo95"], r["hi95"], r["lo90"], r["hi90"], r["p"], r["fdr"], r["raw_label"], cap[key]["capability"],
                    cap[key]["equivalent_allowed"], fin[1][key][0], fin[1][key][1], r3["log2rr"], r3["lo95"], r3["hi95"], r3["fdr"], r3["raw_label"], fin[3][key][0]])
    with open(RESULTS, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(hdr); w.writerows(out)
    log(f"results table {RESULTS} | {len(out)} units | md5 {H.md5_file(RESULTS)}")
    log("Labels are those of run 1 (seed base 31600000), fixed before the result was read; run 3 is shown beside it and changes no label.")
    def line(key, k=1):
        r = runs[k][key]; fl, flag = fin[k][key]
        return (f"log2RR {ff(r['log2rr'])} (RR {ff(2.0 ** pf(r['log2rr']))}; 95% CI {ff(r['lo95'])} to {ff(r['hi95'])}; 90% CI {ff(r['lo90'])} to {ff(r['hi90'])}) p {fp(r['p'])} FDR {fp(r['fdr']) if r['fdr'] else '-'}"
                f" | events glyco {r['n_events_glyco']} bg {r['n_events_bg']} | raw {r['raw_label']} | capability {cap[key]['capability']} | FINAL {fl}" + (f" [{flag}]" if flag else ""))
    pk = ("main", "T1", "both", "pooled")
    log("=== PRIMARY UNIT (main model, T1, N + S/T together, pooled over cancer types) ===")
    log("  " + line(pk))
    log("=== POOLED UNITS, all set-ups (log2 rate ratio of glycosite residues against the same residue type elsewhere in the same protein; > 0 = more mutated) ===")
    for cfg in S8.CONFIGS:
        for t in ("T1", "T2", "T3"):
            for g in ("N", "ST", "both"):
                key = (cfg, t, g, "pooled"); log(f"  {cfg:9s} {t} {g:4s} " + line(key))
    log("=== SET-UPS SIDE BY SIDE (pooled; log2RR and final label) ===")
    ab = {"Enriched": "ENR", "Depleted": "DEP", "Equivalent": "EQV", "Inconclusive": "INC", NOT_TESTED: "NT", NOTSENS: "INS", "Inconclusive, Equivalent not allowed": "ENA"}
    for t in ("T1", "T2", "T3"):
        for g in ("N", "ST", "both"):
            parts = []
            for cfg in S8.CONFIGS:
                key = (cfg, t, g, "pooled"); parts.append(f"{cfg} {ff(runs[1][key]['log2rr'])} {ab.get(fin[1][key][0], fin[1][key][0][:3].upper())}")
            log(f"  {t} {g:4s} | " + " | ".join(parts))
    log("=== PER CANCER TYPE (T1, main model): final labels by residue group ===")
    for g in ("N", "ST", "both"):
        cnt = collections.Counter(fin[1][k][0] for k in keys if k[0] == "main" and k[1] == "T1" and k[2] == g and k[3] != "pooled")
        log(f"  {g:4s}: " + "; ".join(f"{a} {b}" for a, b in sorted(cnt.items())))
    log("  units with a final label Enriched, Depleted or Equivalent (flag = low-sensitivity unit):")
    n_list = 0
    for key in keys:
        if key[0] == "main" and key[1] == "T1" and key[3] != "pooled" and fin[1][key][0] in ("Enriched", "Depleted", "Equivalent"):
            log(f"   {key[3]} {key[2]:4s} " + line(key)); n_list += 1
    if not n_list: log("   none")
    n_oth = sum(1 for key in keys if key[0] != "main" and key[3] == "pooled" and fin[1][key][0] in ("Enriched", "Depleted"))
    log(f"=== pooled units of strict, nohyper, union and norescue with final label Enriched or Depleted: {n_oth} of 36 (listed in the pooled block above) ===")
    log("=== SEED COMPARISON (run 3, seed base 31700000, against run 1, seed base 31600000) ===")
    tested = [k for k in keys if cap[k]["tested"] == "1"]
    dr = [k for k in tested if runs[1][k]["raw_label"] != runs[3][k]["raw_label"]]
    df = [k for k in tested if fin[1][k][0] != fin[3][k][0]]
    est_diff = max([abs(pf(runs[1][k]["log2rr"]) - pf(runs[3][k]["log2rr"])) for k in tested if math.isfinite(pf(runs[1][k]["log2rr"]))] or [0.0])
    def maxdiff(col): return max([abs(pf(runs[1][k][col]) - pf(runs[3][k][col])) for k in tested if math.isfinite(pf(runs[1][k][col])) and math.isfinite(pf(runs[3][k][col]))] or [0.0])
    log(f"  tested units {len(tested)} | raw label differs {len(dr)} | final label differs {len(df)} | largest difference in log2RR {est_diff:.6f} (must be 0: the estimate does not use the bootstrap)")
    log(f"  largest difference in lo95 {maxdiff('lo95'):.4f}, hi95 {maxdiff('hi95'):.4f}, lo90 {maxdiff('lo90'):.4f}, hi90 {maxdiff('hi90'):.4f}")
    fd = [abs(pf(runs[1][k]["fdr"]) - pf(runs[3][k]["fdr"])) for k in tested if runs[1][k]["fdr"] and runs[3][k]["fdr"]]
    log(f"  largest difference in FDR {max(fd) if fd else 0.0:.4f}")
    for k in df: log(f"  differs: {k[0]} {k[1]} {k[2]} {k[3]}: run1 {fin[1][k][0]} -> run3 {fin[3][k][0]}")
    log("Placement (analysis plan section 11, rule 8): units with 'Test can answer' go into Methods and Results; units that could not answer get one line in the main text, methods in the Supplement.")

# ======================= Section 4: Q1 catalogue =======================
CATS = ["loss_T1", "loss_T2", "loss_T3", "loss_any", "gain_plausible", "gain_unknown", "gain_implausible", "gain_any"]
VARCATS = ["loss_T1", "loss_T2", "loss_T3", "gain_plausible", "gain_unknown"]
EXPECT_CLASS = {"loss": 11735, "gain": 18994, "loss_and_gain": 187, "none": 1429876}
EXPECT_GAIN = {"plausible": 3654, "unknown": 13630, "implausible": 1897}
RELNAME = {0: "1 (N)", 1: "2 (X to P)", 2: "3 (S/T)"}

def section4(log):
    require_verdict(); check_inputs(log, with_muts=True)
    cases = {}; cases_total = collections.Counter(); projects = set()
    with open(S8.CASES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cases[r["case_id"]] = (r["project"], r["hypermutated"] == "1"); projects.add(r["project"])
            cases_total[("main", r["project"])] += 1; cases_total[("main", "pooled")] += 1
            if r["hypermutated"] != "1": cases_total[("nohyper", r["project"])] += 1; cases_total[("nohyper", "pooled")] += 1
    rows_all = collections.Counter(); case_all = collections.defaultdict(set)
    rows_ct = collections.Counter(); case_ct = collections.defaultdict(set); var_ct = collections.defaultdict(set)
    vd = {}; pos_ct = collections.Counter(); cls = collections.Counter(); gst = collections.Counter(); guard = collections.Counter()
    n_chosen = 0
    with gzip.open(S8.MUTS, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["in_chosen_file"] != "1": continue
            n_chosen += 1; cls[r["class"]] += 1
            proj = r["project"]; hyp = r["hypermutated"] == "1"
            sets = ["main"] if hyp else ["main", "nohyper"]
            cats = []; lt = r["lost_tiers"]; gs = r["gain_state_best"]
            if lt: cats += ["loss_" + t for t in sorted(set(lt.split(";")))] + ["loss_any"]
            if gs: cats += ["gain_" + gs, "gain_any"]; gst[gs] += 1
            variant = (r["accession"], r["position"], r["ref"], r["alt"]); cid = r["case_id"]
            for s in sets:
                for scope in ("pooled", proj):
                    k = (s, scope); rows_all[k] += 1; case_all[k].add(cid)
                    for c in cats:
                        kc = (s, scope, c); rows_ct[kc] += 1; case_ct[kc].add(cid); var_ct[kc].add(variant)
            for c in cats:
                if c in VARCATS:
                    e = vd.get((c, variant))
                    if e is None: e = vd[(c, variant)] = [r["gene"], 0, 0, collections.Counter()]
                    e[1] += 1; e[2] += (0 if hyp else 1); e[3][proj] += 1
            if lt:
                ns = r["lost_sequon_N"].split(";"); ts = lt.split(";")
                if len(ns) != len(ts): guard["lost_sequon_N and lost_tiers have different numbers of entries"] += 1; continue
                for n_, t_ in zip(ns, ts):
                    rel = int(r["position"]) - int(n_)
                    if rel not in (0, 1, 2): guard["lost sequon position outside 0 to 2 from the N"] += 1; continue
                    for s in sets:
                        for scope in ("pooled", proj): pos_ct[(s, scope, t_, rel)] += 1
    log(f"chosen-file mutation rows read: {n_chosen} | classes {dict(cls)} | gain states of gain and loss_and_gain rows {dict(gst)}")
    log(f"guards (must be 0): {dict(guard) if guard else 0}")
    if guard: raise RuntimeError("Q1 guard failed")
    if dict(cls) != EXPECT_CLASS or dict(gst) != EXPECT_GAIN:
        raise RuntimeError(f"class or gain-state counts differ from the SqnPAN_06 run (expected {EXPECT_CLASS} and {EXPECT_GAIN}): stop and decide")
    log("class counts and gain-state counts equal those of the SqnPAN_06 run: yes")
    scopes = ["pooled"] + sorted(projects)
    os.makedirs(H.TABLES, exist_ok=True)
    with open(Q1_COUNTS, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["set", "scope", "category", "n_mutations", "n_cases", "n_unique_variants", "scope_mutations_total", "scope_cases_with_mutations", "scope_cases_total"])
        for s in ("main", "nohyper"):
            for scope in scopes:
                k = (s, scope); base = [rows_all[k], len(case_all[k]), cases_total[k]]
                w.writerow([s, scope, "all_missense", rows_all[k], len(case_all[k]), "", *base])
                for c in CATS: kc = (s, scope, c); w.writerow([s, scope, c, rows_ct[kc], len(case_ct[kc]), len(var_ct[kc]), *base])
                ku = [(s, scope, "gain_plausible"), (s, scope, "gain_unknown")]
                w.writerow([s, scope, "gain_plausible_or_unknown", sum(rows_ct[x] for x in ku), len(case_ct[ku[0]] | case_ct[ku[1]]), len(var_ct[ku[0]] | var_ct[ku[1]]), *base])
    log(f"Q1 counts table {Q1_COUNTS} | md5 {H.md5_file(Q1_COUNTS)}")
    with open(Q1_POS, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["set", "scope", "tier", "sequon_position", "n_sequon_losses"])
        for s in ("main", "nohyper"):
            for scope in scopes:
                for t in ("T1", "T2", "T3"):
                    for rel in (0, 1, 2): w.writerow([s, scope, t, RELNAME[rel], pos_ct[(s, scope, t, rel)]])
    log(f"Q1 loss-position table {Q1_POS} | md5 {H.md5_file(Q1_POS)}")
    nrec = collections.Counter(); keep = []
    for c in VARCATS:
        items = [(k, v) for k, v in vd.items() if k[0] == c]
        items.sort(key=lambda kv: (-kv[1][1], -len(kv[1][3]), kv[0][1]))
        for rank, (k, v) in enumerate(items):
            if len(v[3]) >= 3 or rank < 30: keep.append((c, k[1], v, rank))
        nrec[c] = sum(1 for k, v in items if len(v[3]) >= 3)
    with open(Q1_VARS, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["category", "accession", "gene", "position", "ref", "alt", "n_cases", "n_cases_without_hypermutated", "n_cancer_types", "cancer_types_with_cases", "rank_in_category"])
        for c, var, v, rank in keep:
            w.writerow([c, var[0], v[0], var[1], var[2], var[3], v[1], v[2], len(v[3]), ";".join(f"{p}:{n}" for p, n in sorted(v[3].items(), key=lambda x: (-x[1], x[0]))), rank + 1])
    log(f"Q1 variant table {Q1_VARS} | {len(keep)} rows | md5 {H.md5_file(Q1_VARS)}")
    log("=== Q1 pooled counts (mutations / cases / unique variants) ===")
    for s in ("main", "nohyper"):
        k = (s, "pooled"); log(f"  {s}: all missense rows {rows_all[k]} in {len(case_all[k])} cases (cases in the case table {cases_total[k]})")
        for c in CATS: kc = (s, "pooled", c); log(f"    {c:17s} {rows_ct[kc]:7d} / {len(case_ct[kc]):6d} / {len(var_ct[kc]):6d}")
    log("=== Q1 variants seen in 3 or more cancer types (main set) ===")
    for c in VARCATS: log(f"  {c}: {nrec[c]} variants")
    log("=== Q1 top 10 variants by cases, loss_T1 and gain_plausible (main set; descriptive, not blind: Fan et al. 2018) ===")
    for c in ("loss_T1", "gain_plausible"):
        for cc, var, v, rank in [x for x in keep if x[0] == c][:10]:
            log(f"  {c} #{rank + 1}: {v[0]} {var[2]}{var[1]}{var[3]} ({var[0]}) cases {v[1]} (without hypermutated {v[2]}) cancer types {len(v[3])}")
    log("=== Q1 position of the lost sequon hit, pooled (sequon losses; 1 = N, 2 = X to P, 3 = S/T) ===")
    for s in ("main", "nohyper"):
        for t in ("T1", "T2", "T3"): log(f"  {s} {t}: " + " | ".join(f"{RELNAME[rel]} {pos_ct[(s, 'pooled', t, rel)]}" for rel in (0, 1, 2)))
    log("Raw counts are not rates: they depend on protein length and mutation load. The Q2 selection test is the test; this catalogue is descriptive.")

SECTIONS = {"1": ("selfcheck_inputs", section1), "2": ("real_Q2_runs", section2), "3": ("labels_results_seeds", section3), "4": ("Q1_catalogue", section4)}

def run_one(k):
    name, fn = SECTIONS[k]; ok = {"v": False}
    def wrapped(log):
        fn(log); ok["v"] = True
    H.run_section("P2_09", k, name, wrapped)
    return ok["v"]

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for k in (list(SECTIONS) if which == "all" else [which]):
        if not run_one(k): print(f"SECTION {k} DID NOT FINISH: the chain stops here (see the log)", flush=True); break
# ---- END OF SCRIPT ----
