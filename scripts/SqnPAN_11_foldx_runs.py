# SqnPAN_11_foldx_runs.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-10-25; D1 A, D2 A, D3 A, D4 B and the fixed details;
#   the analysis plan sections 8 and 10 with Amendment 3; set-up and timing in SqnPAN_10_foldx_setup.py, verdict PASSED 2026-10-25).
# Question (analysis plan Q3, FoldX part): do T1 sequon-loss substitutions change the folding stability of the protein (FoldX ddG, glycan not modelled) more or less than matched
#   control substitutions in the same protein?
# FoldX: version 5.1 (foldx_1_20270131.exe, md5 c66ffd383d4d1d274e940bf99979bf00, expires 2026-12-31), all settings at the program defaults (decision D3 A): the physics
#   values are passed explicitly (298 K, pH 7, ionic strength 0.05, -IGNORE, vdwDesign 2, moveNeighbours true; SqnPAN_10 check D proved them equal to the defaults), --out-pdb is not passed
#   (default, mutant PDBs written and deleted after each call; SqnPAN_10 amendment, option B); RepairPDB once per model (protein x fragment), then one BuildModel call per model with all its
#   mutations, numberOfRuns 3; ddG of a mutation = mean of its 3 runs (checked against the Average table within 0.01). Functions from SqnPAN_10_foldx_helpers.py and SqnPAN_10_foldx_setup.py.
# FoldX set (analysis plan section 8): unique T1 loss substitutions of the chosen files, protein usable (rule 1A), pLDDT >= 70 at the position, model fragment by rule 2A.
#   D3 A: the PRIMARY comparison uses losses at the N (SqnPAN_06 res_role seqN_T1, reference N) or the S/T (seqST_T1, reference S or T) of a T1 sequon; the other losses (X to P, and rows whose
#   SqnPAN_06 res_role is not seqN_T1 or seqST_T1) are built and reported one by one, descriptively.
# Matched controls (D2 A): observed TCGA missense substitutions of the chosen files (unique substitutions, class none, so no sequon change) in the same protein and the same model fragment,
#   with the same wild-type amino acid exactly (N with N, S with S, T with T), at a position in no N-X-S/T (X not P) sequon of the UniProt 2026_03 sequence and not at its middle position X
#   (as the Q2 background, SqnPAN_08 decision 1A), pLDDT >= 70; for S and T only substitutions to a residue that is not S or T. Two controls per loss, drawn without replacement within the protein
#   (losses taken in sorted order); a loss with 1 control keeps it, a loss with 0 controls is counted and left out of the paired comparison.
#   Stricter model side by side: controls must also have the same mutant residue (separate draw; losses with at least 1 such control).
# Comparison (D1 A): for each primary loss with at least 1 control, d = ddG(loss) - mean ddG(its controls). Statistic: median of d. CIs (95% and 90%, percentile) from 2,000 bootstrap
#   resamples over proteins (all pairs of a protein move together). Labels (kcal/mol; > 0 = losses more destabilising): Higher = 95% CI entirely above +0.5; Lower = 95% CI entirely
#   below -0.5; Equivalent = 90% CI inside -0.5 to +0.5; Inconclusive = otherwise. A unit with fewer than 20 proteins is "Not tested (too few data)" (rule fixed now, before any ddG).
#   Units: (1) PRIMARY: N + S/T losses, draw 1, all models; (2) N losses; (3) S/T losses; (4) without fragment models; (5) stricter matching (same mutant residue); (6) draw 2 (seed check of
#   the control draw). Also reported per unit: median ddG of the losses and of the controls, and the shares above 1 kcal/mol (destabilising) and above 3 kcal/mol (highly destabilising;
#   Tokuriki et al. 2007, analysis plan Amendment 3) among losses and among the unique control substitutions.
# Control C4 (analysis plan section 10, run and judged before any loss or control mutation is built): 50 usable proteins drawn at random among those with at least one buried hydrophobic residue
#   (I, L, V, F or M; pLDDT >= 90; relative exposure <= 0.05; not terminal) and at least one exposed lysine (relative exposure >= 0.5; pLDDT >= 70; not terminal), the SqnPAN_10 test protein P14174
#   left out; one buried hydrophobic residue changed to D and one exposed lysine changed to R per protein. PASSED if >= 80% of the buried group and <= 20% of the lysine group have
#   ddG > 1 kcal/mol, judged only if at least 45 of each group have a ddG (rule fixed now). If C4 fails, the comparison is labelled "Inconclusive, test not sensitive" and section 3 does not run.
# Gains (D4 B; descriptive only, after the loss result is written): sequon-gain substitutions of the chosen files in proteins with at least one T1 sequon, usable, pLDDT >= 70; plausible
#   first, then unknown (Amendment 2; implausible not built). Reported for plausible and for plausible plus unknown: number, median ddG, shares above 1 and 3 kcal/mol.
# Seeds: control draw 1 SeedSequence(31900000 + 1), drawn twice (must be identical); stricter draw 31900000 + 2; C4 draw 31900000 + 3; control draw 2 SeedSequence(32000000 + 1);
#   bootstrap SeedSequence([32100000, unit]) for run 1 and run 2 (tables must be identical apart from the run column) and [32200000, unit] for run 3. Labels come from run 1; run 3 and draw 2
#   are reported beside it and change no label. All bases at least 100,000 apart from 31000000 to 31800000 (checked).
# Parallel work: 8 FoldX processes at once (threads, one FoldX process each), largest models first. Each finished model writes one result file (written last, so a model is done only with
#   it); a stopped run continues where it stopped when the same section is started again; failed models are tried once more at the end of the section.
#   Measured time (SqnPAN_10 section 2, a projection): core about 15 h, plus plausible gains about 11 h, plus unknown gains about 15 h with 8 processes.
#   PERFORMANCE AMENDMENT (2026-10-25, to use more of the CPU; made during section 3 before any ddG was read; no rule, label, threshold or seed changes):
#     the number of FoldX processes at once is now 16 (WORKERS; the environment variable P2_WORKERS overrides it). The section 3 log and Task Manager showed 8 processes using
#     8 of the 24 CPU threads (37% CPU, 59% memory, FoldX about 280 MB each). A ddG does not depend on how many processes run at once (SqnPAN_10 checks C and D: separate FoldX calls gave
#     identical tables). The "rough time left" line is now weighted by model size (fitted seconds from the SqnPAN_10 timing fit); the count-based line said 43 h because the largest
#     models run first. Finished models are kept, so the section continues where the 8-process run stopped.
# Order: section 1 inputs, set, control draws, C4 draw, plan table and time projection (no FoldX); section 2 control C4 (FoldX) and verdict; section 3 core FoldX runs (losses and all controls;
#   no value printed); section 4 comparison, labels, runs side by side; section 5 gains (FoldX) and the descriptive gain summary. "all" chains the five sections and stops at the first one
#   that fails; section 3 refuses to run unless C4 PASSED; sections 2 to 5 check the md5 of the plan table written by section 1.
# Usage (from the repository folder):  python scripts/SqnPAN_11_foldx_runs.py all      (or one section: 1, 2, 3, 4 or 5)
# Logs: results\logs\P2_11_section<k>_<name>.txt. Outputs: data\processed\P2_11_foldx_plan.tsv and P2_11_plan_md5.txt, data\processed\P2_11_C4_verdict.txt,
#   data\processed\P2_10_foldx\P2_11_results\<section>\ (one file per model, with copies of the FoldX Dif and Average tables), data\processed\P2_10_foldx\repaired\ (repaired models),
#   results\tables\P2_11_foldx_ddg.tsv, P2_11_comparison_run1.tsv to run3.tsv, P2_11_comparison_results.tsv, P2_11_foldx_gains_ddg.tsv, P2_11_foldx_gains_summary.tsv.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, csv, gzip, math, time, shutil, hashlib, datetime, collections, threading
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_10_foldx_helpers as FX
import SqnPAN_10_foldx_setup as S10

DRAW1, STRICT, C4SEED, DRAW2 = 31900000 + 1, 31900000 + 2, 31900000 + 3, 32000000 + 1
BOOT = {1: 32100000, 2: 32100000, 3: 32200000}
NEW_BASES = [31900000, 32000000, 32100000, 32200000]
USED_BASES = [31000000, 31100000, 31200000, 31300000, 31400000, 31500000, 31600000, 31700000, 31800000]
WORKERS = int(os.environ.get("P2_WORKERS", "16"))             # 8 until the performance amendment of 2026-10-25 (see the header)
N_BOOT = int(os.environ.get("P2_NBOOT", "2000"))             # P2_NBOOT is used only for quick test runs
MARGIN = 0.5; DESTAB = 1.0; HIGH = 3.0; MIN_PROT = 20
TEST_PROTEIN = "P14174"
REAL_TABLES = os.path.join(H.PAPER2, "results", "tables")
EXTRA_INPUTS = {os.path.join(HERE, "SqnPAN_10_foldx_helpers.py"): "4cf976f15c5ef96c7a7e9edbc12fb05a",
                os.path.join(HERE, "SqnPAN_10_foldx_setup.py"): "197110e63d53cde430d494ec63e81708",
                os.path.join(REAL_TABLES, "P2_10_timing.tsv"): "7371feab33c984d637ab1ddf46f2e636"}
P210_VERDICT = os.path.join(H.PROC, "P2_10_setup_verdict.txt")
PLAN = os.path.join(H.PROC_OUT, "P2_11_foldx_plan.tsv"); PLAN_MD5 = os.path.join(H.PROC_OUT, "P2_11_plan_md5.txt")
C4_VERDICT = os.path.join(H.PROC_OUT, "P2_11_C4_verdict.txt")
RES_ROOT = os.path.join(FX.FOLDX_WORK, "P2_11_results"); WORK_ROOT = os.path.join(FX.FOLDX_WORK, "P2_11_work")
DDG_TSV = os.path.join(H.TABLES, "P2_11_foldx_ddg.tsv")
RUN_TSV = os.path.join(H.TABLES, "P2_11_comparison_run{}.tsv"); CMP_TSV = os.path.join(H.TABLES, "P2_11_comparison_results.tsv")
GAIN_DDG = os.path.join(H.TABLES, "P2_11_foldx_gains_ddg.tsv"); GAIN_SUM = os.path.join(H.TABLES, "P2_11_foldx_gains_summary.tsv")
PLAN_COLS = ["set", "key", "accession", "gene", "fragment", "position", "model_number", "wt", "mt", "role", "plddt", "rsa", "fragment_model", "loss_key"]
CORE_SETS = ("loss", "ctrl_d1", "ctrl_d2", "ctrl_strict"); C4_SETS = ("c4_buried", "c4_lysine")
UNITS = [(1, "primary: N + S/T losses, draw 1, all models", "ctrl_d1", ("N", "ST"), True),
         (2, "N losses, draw 1", "ctrl_d1", ("N",), True),
         (3, "S/T losses, draw 1", "ctrl_d1", ("ST",), True),
         (4, "N + S/T losses, draw 1, without fragment models", "ctrl_d1", ("N", "ST"), False),
         (5, "stricter matching (same mutant residue)", "ctrl_strict", ("N", "ST"), True),
         (6, "N + S/T losses, draw 2 (seed check of the control draw)", "ctrl_d2", ("N", "ST"), True)]

def mkey(a, fr, p, wt, mt): return f"{a}|F{fr}|{p}|{wt}{mt}"

# ------------------------------------------------------------------ shared ------------------------------------------------------------------
def check_plan(log):
    if not os.path.exists(PLAN) or not os.path.exists(PLAN_MD5): raise RuntimeError("no plan table: run section 1 first")
    want = open(PLAN_MD5).read().split()[0]; got = H.md5_file(PLAN)
    if want != got: raise RuntimeError(f"plan table md5 {got} differs from the md5 written by section 1 ({want}): stop")
    with open(PLAN, newline="", encoding="utf-8") as f: rows = list(csv.DictReader(f, delimiter="\t"))
    log(f"plan table {PLAN}: {len(rows)} rows, md5 {got} (as written by section 1)")
    return rows

def result_file(section, m): return os.path.join(RES_ROOT, section, f"{m[0]}_F{m[1]}.tsv")

def read_result(path):
    """{key: (mean ddG, sd, run texts)} of one model result file (empty if absent)."""
    out = {}
    if not os.path.exists(path): return out
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"): out[r["key"]] = (float(r["ddg_mean"]), r["sd_average"], (r["run0"], r["run1"], r["run2"]))
    return out

def load_results(sections):
    out = {}
    for s in sections:
        d = os.path.join(RES_ROOT, s)
        if not os.path.isdir(d): continue
        for fn in sorted(os.listdir(d)):
            if fn.endswith(".tsv"): out.update(read_result(os.path.join(d, fn)))
    return out

def run_model(section, m, muts, seqs, ranges):
    """RepairPDB (or the kept repaired model) + one BuildModel call with all mutations of model m. Writes the result file last. Returns (m, ok, info)."""
    acc, fr = m; t0 = time.time()
    wd = os.path.join(WORK_ROOT, section, f"{acc}_F{fr}")
    try:
        if os.path.isdir(wd): shutil.rmtree(wd)
        os.makedirs(wd)
        rep_keep = os.path.join(FX.REPAIRED_DIR, f"{acc}_F{fr}_Repair.pdb")
        if not os.path.exists(rep_keep):
            path = S10.write_model(acc, fr, seqs, ranges, wd)[0]
            r1 = FX.repair(path, wd)
            if not (r1["rc"] == 0 and r1["exists"]): return m, False, "RepairPDB failed: " + " / ".join(FX.stdout_tail(r1["stdout"], 3))
            os.makedirs(FX.REPAIRED_DIR, exist_ok=True); tmp = rep_keep + ".tmp"
            shutil.copy2(r1["repaired"], tmp); os.replace(tmp, rep_keep)
        r2 = FX.build(rep_keep, [(x["wt"], "A", int(x["model_number"]), x["mt"]) for x in muts], wd, tag="BuildModel")
        if r2["rc"] != 0 or not os.path.exists(r2["dif"]): return m, False, f"BuildModel rc {r2['rc']}: " + " / ".join(FX.stdout_tail(r2["stdout"], 3))
        dif = FX.dif_ddg(r2["dif"]); avg = FX.average_ddg(r2["average"])
        rows = []
        for i, x in enumerate(muts, start=1):
            runs = dif.get(i, {})
            if sorted(runs) != list(range(FX.N_RUNS)) or i not in avg: return m, False, f"mutation {i} ({x['key']}) incomplete in the FoldX tables"
            mean = float(np.mean([float(runs[r]) for r in range(FX.N_RUNS)]))
            if abs(mean - float(avg[i][1])) > 0.01: return m, False, f"mutation {i}: Dif mean {mean} differs from Average {avg[i][1]}"
            rows.append([x["key"], acc, fr, x["position"], x["model_number"], x["wt"], x["mt"]] + [runs[r] for r in range(FX.N_RUNS)] + [f"{mean:.6f}", avg[i][0], avg[i][1]])
        resdir = os.path.join(RES_ROOT, section); os.makedirs(os.path.join(resdir, "fxout"), exist_ok=True)
        for p, lab in ((r2["dif"], "Dif"), (r2["average"], "Average")): shutil.copy2(p, os.path.join(resdir, "fxout", f"{acc}_F{fr}_{lab}.fxout"))
        out = result_file(section, m); tmp = out + ".tmp"
        with open(tmp, "w", newline="", encoding="utf-8") as o:
            w = csv.writer(o, delimiter="\t", lineterminator="\n")
            w.writerow(["key", "accession", "fragment", "position", "model_number", "wt", "mt", "run0", "run1", "run2", "ddg_mean", "sd_average", "mean_average"])
            w.writerows(rows)
        os.replace(tmp, out)
        shutil.rmtree(wd, ignore_errors=True)
        return m, True, f"{len(muts)} mutations, {time.time() - t0:.0f} s"
    except Exception as ex:
        return m, False, f"{type(ex).__name__}: {ex}"

def run_jobs(log, section, rows, seqs, ranges, skip_keys=frozenset()):
    """Runs every model of the rows (plan dicts) not already done in this section. Two passes: the second retries failed models once."""
    jobs = {}
    for r in rows:
        if r["key"] in skip_keys: continue
        jobs.setdefault((r["accession"], int(r["fragment"])), {})[r["key"]] = r
    mlen = lambda m: [e - b + 1 for f, b, e in ranges[m[0]] if f == m[1]][0]
    wt = lambda m: 0.373 * mlen(m) + len(jobs[m]) * 0.338 * mlen(m) ** 0.48      # fitted seconds (P2_10 timing fit), used only for the progress estimate
    todo = []
    for m in sorted(jobs, key=lambda m: (-mlen(m), m)):
        have = read_result(result_file(section, m))
        if not set(jobs[m]) <= set(have): todo.append(m)
    n_mut = sum(len(v) for v in jobs.values())
    log(f"[{section}] models {len(jobs)} with {n_mut} unique mutations; already done {len(jobs) - len(todo)}; to run {len(todo)}; {WORKERS} FoldX processes at once")
    failed = []
    for attempt in (1, 2):
        if not todo: break
        if attempt == 2: log(f"[{section}] second try for {len(todo)} failed model(s)")
        t0 = time.time(); failed = []; ndone = 0; last = t0; wdone = 0.0; wtot = sum(wt(m) for m in todo)
        with ThreadPoolExecutor(max_workers=WORKERS) as ex:
            fut = {ex.submit(run_model, section, m, [jobs[m][k] for k in sorted(jobs[m])], seqs, ranges): m for m in todo}
            pending = set(fut)
            while pending:
                done, pending = wait(pending, timeout=60, return_when=FIRST_COMPLETED)
                for f in done:
                    m, ok, info = f.result(); ndone += 1; wdone += wt(m)
                    if not ok: failed.append(m); log(f"  FAILED {m[0]} F{m[1]}: {info}")
                if time.time() - last >= 600 or not pending:
                    el = time.time() - t0; last = time.time()
                    eta = el / wdone * (wtot - wdone) if wdone else float("nan")
                    log(f"  [{section}] {time.strftime('%H:%M')} done {ndone} of {len(todo)} models ({100 * wdone / wtot:.0f}% of the work by model size), failed {len(failed)}, "
                        f"elapsed {el / 3600:.2f} h, rough time left {eta / 3600:.1f} h (weighted by model size)")
        todo = failed
    log(f"[{section}] finished; models still failed after two tries: {len(failed)}" + ("" if not failed else " (" + ", ".join(f"{a} F{fr}" for a, fr in failed[:20]) + ")"))
    return failed

# ------------------------------------------------------------------ section 1 ------------------------------------------------------------------
def section1(log):
    S10.check_inputs(log)
    for b in NEW_BASES:
        for u in USED_BASES + [x for x in NEW_BASES if x != b]:
            if abs(b - u) < 100000: raise ValueError(f"seed base {b} is closer than 100,000 to {u}")
    log(f"seed bases {NEW_BASES}: at least 100,000 apart from each other and from {USED_BASES} (checked)")
    bad = 0
    for p, m in EXTRA_INPUTS.items():
        got = H.md5_file(p) if os.path.exists(p) else "MISSING"; bad += got != m
        log(f"  {'OK ' if got == m else 'BAD'} {p} md5 {got}" + ("" if got == m else f" (expected {m})"))
    if bad: raise ValueError("SqnPAN_10 files differ from the recorded md5: stop")
    if not (os.path.exists(P210_VERDICT) and open(P210_VERDICT).read().startswith("PASSED")): raise RuntimeError("SqnPAN_10 set-up verdict is not PASSED")
    log("SqnPAN_10 set-up verdict: " + open(P210_VERDICT).readline().strip())
    use, ranges, seqs = S10.load_tables()
    usable = {a for a, r in use.items() if r["status"] == "use"}
    nT1 = {}
    with open(os.path.join(H.PROC, "P2_06_proteins.tsv"), newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"): nT1[r["accession"]] = int(r["n_T1_sequons"])
    loss = {}; gain = {}; pool = {}
    with gzip.open(os.path.join(H.PROC, "P2_06_mutations.tsv.gz"), "rt", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["in_chosen_file"] != "1": continue
            k = (r["accession"], int(r["position"]), r["ref"], r["alt"])
            if "T1" in r["lost_tiers"].split(";"): loss[k] = r["res_role"]
            if r["class"] in ("gain", "loss_and_gain"): gain[k] = r["gain_state_best"]
            if r["class"] == "none" and r["ref"] in ("N", "S", "T"): pool[k] = r["res_role"]
    loss_acc = {k[0] for k in loss if k[0] in usable}
    need = collections.defaultdict(set)
    for k in loss:
        if k[0] in usable: need[k[0]].add(k[1])
    for k in pool:
        if k[0] in loss_acc: need[k[0]].add(k[1])
    for k, st in gain.items():
        if k[0] in usable and nT1.get(k[0], 0) > 0 and st in ("plausible", "unknown"): need[k[0]].add(k[1])
    feat = {}; c4b = collections.defaultdict(list); c4k = collections.defaultdict(list)
    with gzip.open(os.path.join(H.PROC, "P2_07_residues.tsv.gz"), "rt", newline="") as f:
        rd = csv.reader(f, delimiter="\t"); hdr = next(rd)
        if hdr != ["accession", "position", "aa", "plddt", "asa", "rsa", "fragment", "n_fragments", "terminal", "cut_distance"]: raise ValueError(f"residue table header {hdr}")
        for a, p, aa, pl, asa, rsa, fr, nf, term, cut in rd:
            p = int(p); pl = float(pl); rs = float(rsa) if rsa else None; fr = int(fr)
            s = need.get(a)
            if s is not None and p in s: feat[(a, p)] = (aa, pl, rs, fr, term == "1")
            if term == "1" or rs is None or a not in usable: continue
            if aa in "ILVFM" and pl >= 90 and rs <= 0.05: c4b[a].append((p, aa, fr, pl, rs))
            elif aa == "K" and rs >= 0.5 and pl >= 70: c4k[a].append((p, aa, fr, pl, rs))
    log(f"residue features read for {len(feat)} mutated positions; C4 candidates: {len(c4b)} proteins with a buried hydrophobic residue, {len(c4k)} with an exposed lysine")
    def begin(a, fr): return [b for f2, b, e in ranges[a] if f2 == fr][0]
    def row(st, a, p, wt, mt, role="", loss_key=""):
        aa, pl, rs, fr, term = feat[(a, p)] if (a, p) in feat else (None, None, None, None, None)
        if aa != wt or seqs[a][p - 1] != wt: raise ValueError(f"wild type {wt} at {a} {p} differs from the residue table ({aa}) or the sequence ({seqs[a][p - 1]})")
        return {"set": st, "key": mkey(a, fr, p, wt, mt), "accession": a, "gene": use[a]["gene"], "fragment": fr, "position": p, "model_number": p - begin(a, fr) + 1,
                "wt": wt, "mt": mt, "role": role, "plddt": f"{pl:.2f}", "rsa": "" if rs is None else f"{rs:.3f}", "fragment_model": "yes" if use[a]["n_fragments"] != "1" else "no",
                "loss_key": loss_key}
    # losses
    plan = []; cnt = collections.Counter(); losses = []
    for k, role in sorted(loss.items()):
        a, p, ref, alt = k
        if a not in usable: cnt["protein not usable"] += 1; continue
        if feat[(a, p)][1] < 70: cnt["pLDDT < 70"] += 1; continue
        kind = "N" if (role == "seqN_T1" and ref == "N") else ("ST" if (role == "seqST_T1" and ref in "ST") else "odd")
        r = row("loss", a, p, ref, alt, role=kind); losses.append(r); plan.append(r); cnt[f"eligible {kind}"] += 1
    log(f"T1 losses (unique substitutions, chosen files): {len(loss)}; {dict(cnt)}")
    odd = [r for r in losses if r["role"] == "odd"]
    log(f"  descriptive only (D3 A): {len(odd)} = " + ", ".join(f"{r['gene']} {r['wt']}{r['position']}{r['mt']} (SqnPAN_06 res_role {loss[(r['accession'], int(r['position']), r['wt'], r['mt'])]})" for r in odd))
    # control pools
    excl = {}
    for a in loss_acc:
        e = set()
        for n in H.sequon_positions(seqs[a])[0]: e.update((n, n + 1, n + 2))
        excl[a] = e
    pools = collections.defaultdict(list); disagree = 0
    for k, role in sorted(pool.items()):
        a, p, ref, alt = k
        if a not in loss_acc: continue
        insq = p in excl[a]
        if (role in ("nonseqN", "nonseqST")) == insq: disagree += 1
        if insq: continue
        x = feat[(a, p)]
        if x[1] < 70: continue
        if ref in "ST" and alt in "ST": continue
        pools[(a, x[3], ref)].append((p, alt))
    log(f"control pools: {sum(len(v) for v in pools.values())} candidate substitutions in {len({k[0] for k in pools})} proteins; rows where the SqnPAN_06 res_role (nonseqN/nonseqST) and our "
        f"sequon check disagree: {disagree} (middle-X positions are background in res_role but excluded here, as in SqnPAN_08 decision 1A)")
    prim = [r for r in losses if r["role"] in ("N", "ST")]
    def draw(seed, strict):
        rng = np.random.default_rng(np.random.SeedSequence(seed)); used = collections.defaultdict(set); out = {}
        for L in prim:
            a = L["accession"]
            cand = sorted(c for c in pools.get((a, L["fragment"], L["wt"]), []) if c not in used[a] and (not strict or c[1] == L["mt"]))
            k = min(2, len(cand)); pick = [cand[i] for i in sorted(rng.choice(len(cand), size=k, replace=False))] if k else []
            used[a].update(pick); out[L["key"]] = pick
        return out
    d1 = draw(DRAW1, False); d1b = draw(DRAW1, False); d2 = draw(DRAW2, False); ds = draw(STRICT, True)
    log(f"control draw 1 repeated with the same seed is identical: {d1 == d1b}")
    if d1 != d1b: raise RuntimeError("control draw 1 is not reproducible")
    for name, d, st in (("draw 1", d1, "ctrl_d1"), ("draw 2", d2, "ctrl_d2"), ("stricter", ds, "ctrl_strict")):
        dist = collections.Counter(len(v) for v in d.values())
        log(f"  {name}: losses with 2 / 1 / 0 controls: {dist[2]} / {dist[1]} / {dist[0]}; unique control substitutions {len({(L.split('|')[0], c) for L, v in d.items() for c in v})}")
        for L in prim:
            for p, alt in d[L["key"]]: plan.append(row(st, L["accession"], p, L["wt"], alt, role=L["role"], loss_key=L["key"]))
    s1 = {(L, c) for L, v in d1.items() for c in v}; s2 = {(L, c) for L, v in d2.items() for c in v}
    log(f"  draw 1 and draw 2 share {len(s1 & s2)} of {len(s1)} loss-control pairs")
    # C4
    rng = np.random.default_rng(np.random.SeedSequence(C4SEED))
    elig = sorted(a for a in usable if c4b.get(a) and c4k.get(a) and a != TEST_PROTEIN)
    chosen = sorted(elig[i] for i in rng.choice(len(elig), size=50, replace=False))
    for a in chosen:
        pb, aab, frb, plb, rsb = c4b[a][int(rng.integers(len(c4b[a])))]; pk, aak, frk, plk, rsk = c4k[a][int(rng.integers(len(c4k[a])))]
        feat[(a, pb)] = (aab, plb, rsb, frb, False); feat[(a, pk)] = (aak, plk, rsk, frk, False)
        rb = row("c4_buried", a, pb, aab, "D"); rk = row("c4_lysine", a, pk, "K", "R")
        plan.extend([rb, rk])
    log(f"control C4: {len(elig)} eligible proteins, 50 drawn; {len({(r['accession'], r['fragment']) for r in plan if r['set'] in C4_SETS})} models")
    # gains
    gc = collections.Counter()
    for k, st in sorted(gain.items()):
        a, p, ref, alt = k
        if st not in ("plausible", "unknown") or a not in usable or nT1.get(a, 0) == 0: continue
        if feat[(a, p)][1] < 70: continue
        plan.append(row("gain_" + st, a, p, ref, alt)); gc[st] += 1
    log(f"gains in T1 glycoproteins, usable, pLDDT >= 70: {dict(gc)}")
    with open(PLAN, "w", newline="", encoding="utf-8") as o:
        w = csv.DictWriter(o, fieldnames=PLAN_COLS, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(plan)
    md5 = H.md5_file(PLAN)
    with open(PLAN_MD5, "w") as o: o.write(md5 + "\n")
    log(f"plan table {PLAN}: {len(plan)} rows, md5 {md5}")
    # time projection from the SqnPAN_10 timing table
    with open(os.path.join(REAL_TABLES, "P2_10_timing.tsv"), newline="") as f: tt = [r for r in csv.DictReader(f, delimiter="\t") if not r["error"]]
    Ls = np.log([float(r["model_length"]) for r in tt])
    fr_ = np.polyfit(Ls, np.log([float(r["repair_seconds"]) for r in tt]), 1); fb_ = np.polyfit(Ls, np.log([float(r["build_seconds_3_mutations"]) / 3 for r in tt]), 1)
    def mlen(m): return [e - b + 1 for f2, b, e in ranges[m[0]] if f2 == m[1]][0]
    seen = {f.replace("_Repair.pdb", "") for f in os.listdir(FX.REPAIRED_DIR)} if os.path.isdir(FX.REPAIRED_DIR) else set()
    done_keys = set()
    for lab, sets in (("section 2 (C4)", C4_SETS), ("section 3 (core)", CORE_SETS), ("section 5 (plausible gains)", ("gain_plausible",)), ("section 5 (unknown gains)", ("gain_unknown",))):
        rs = [r for r in plan if r["set"] in sets and r["key"] not in done_keys]
        ms = collections.Counter((r["accession"], int(r["fragment"])) for r in {r["key"]: r for r in rs}.values())
        sec = 0.0; longest = 0.0
        for m, c in ms.items():
            t = (0 if f"{m[0]}_F{m[1]}" in seen else math.exp(fr_[1] + fr_[0] * math.log(mlen(m)))) + c * math.exp(fb_[1] + fb_[0] * math.log(mlen(m)))
            sec += t; longest = max(longest, t); seen.add(f"{m[0]}_F{m[1]}")
        done_keys |= {r["key"] for r in rs}
        log(f"  projection {lab}: {sum(ms.values())} mutations in {len(ms)} models, {sec / 3600:.1f} CPU h, about {max(sec / WORKERS, longest) / 3600:.1f} h with {WORKERS} processes (a guess)")

# ------------------------------------------------------------------ section 2 ------------------------------------------------------------------
def section2(log):
    plan = check_plan(log); use, ranges, seqs = S10.load_tables()
    S10.check_inputs(log)
    rows = [r for r in plan if r["set"] in C4_SETS]
    failed = run_jobs(log, "c4", rows, seqs, ranges)
    dd = load_results(["c4"])
    res = {}
    for st in C4_SETS:
        v = [dd[r["key"]][0] for r in rows if r["set"] == st and r["key"] in dd]
        res[st] = v
        log(f"C4 {st}: {len(v)} of 50 with a ddG; median {np.median(v):.2f} kcal/mol; share > {DESTAB} kcal/mol {np.mean(np.array(v) > DESTAB):.2f}; share > {HIGH} {np.mean(np.array(v) > HIGH):.2f}" if v else f"C4 {st}: no ddG")
    enough = all(len(res[s]) >= 45 for s in C4_SETS)
    fb = float(np.mean(np.array(res["c4_buried"]) > DESTAB)) if res["c4_buried"] else 0.0
    fk = float(np.mean(np.array(res["c4_lysine"]) > DESTAB)) if res["c4_lysine"] else 1.0
    ok = enough and fb >= 0.80 and fk <= 0.20
    with open(C4_VERDICT, "w", newline="\n") as o:
        o.write(("PASSED" if ok else "FAILED") + f" {time.strftime('%Y-%m-%d %H:%M:%S')} P2_11 control C4\n")
        o.write(f"buried hydrophobic to D: {len(res['c4_buried'])} with ddG, share > 1 kcal/mol {fb:.3f} (needed >= 0.80)\n")
        o.write(f"exposed lysine to R: {len(res['c4_lysine'])} with ddG, share > 1 kcal/mol {fk:.3f} (needed <= 0.20)\n")
    log(f"C4 VERDICT: {'PASSED' if ok else 'FAILED'} (buried share {fb:.3f} >= 0.80; lysine share {fk:.3f} <= 0.20; at least 45 of each with a ddG: {enough}) -> {C4_VERDICT}")
    if not ok: raise RuntimeError("control C4 FAILED: the FoldX comparison is 'Inconclusive, test not sensitive'; section 3 will not run")

# ------------------------------------------------------------------ section 3 ------------------------------------------------------------------
def section3(log):
    plan = check_plan(log)
    if not (os.path.exists(C4_VERDICT) and open(C4_VERDICT).read().startswith("PASSED")): raise RuntimeError("control C4 has not PASSED: section 3 does not run")
    log("C4 verdict: " + open(C4_VERDICT).readline().strip())
    S10.check_inputs(log); use, ranges, seqs = S10.load_tables()
    rows = [r for r in plan if r["set"] in CORE_SETS]
    run_jobs(log, "core", rows, seqs, ranges)
    dd = load_results(["core"]); keys = {r["key"] for r in rows}
    log(f"core mutations with a ddG: {len(keys & set(dd))} of {len(keys)} (no ddG value is printed in this section)")

# ------------------------------------------------------------------ section 4 ------------------------------------------------------------------
def boot_unit(pairs, seed, B):
    by = collections.defaultdict(list)
    for a, d in pairs: by[a].append(d)
    accs = sorted(by); arrs = [np.array(by[a]) for a in accs]
    est = float(np.median(np.concatenate(arrs)))
    rng = np.random.default_rng(np.random.SeedSequence(seed)); bs = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, len(accs), len(accs)); bs[b] = np.median(np.concatenate([arrs[i] for i in idx]))
    return est, np.percentile(bs, [2.5, 97.5]), np.percentile(bs, [5, 95])

def label(c95, c90):
    if c95[0] > MARGIN: return "Higher"
    if c95[1] < -MARGIN: return "Lower"
    if c90[0] > -MARGIN and c90[1] < MARGIN: return "Equivalent"
    return "Inconclusive"

RUN_COLS = ["run", "unit", "description", "control_set", "n_pairs", "n_proteins", "n_losses_2_controls", "median_d", "ci95_lo", "ci95_hi", "ci90_lo", "ci90_hi", "label",
            "median_ddg_losses", "median_ddg_controls", "losses_gt1", "losses_gt3", "controls_gt1", "controls_gt3", "n_unique_controls"]

def section4(log):
    plan = check_plan(log)
    if not (os.path.exists(C4_VERDICT) and open(C4_VERDICT).read().startswith("PASSED")): raise RuntimeError("control C4 has not PASSED")
    log("C4 verdict: " + open(C4_VERDICT).readline().strip())
    dd = load_results(["c4", "core"])
    core = [r for r in plan if r["set"] in CORE_SETS]
    miss = {r["key"] for r in core if r["key"] not in dd}
    log(f"core mutations without a ddG (failed models): {len(miss)} of {len({r['key'] for r in core})}")
    with open(DDG_TSV, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(PLAN_COLS + ["ddg_mean", "ddg_sd", "run0", "run1", "run2"])
        for r in plan:
            if r["set"] in CORE_SETS or r["set"] in C4_SETS:
                v = dd.get(r["key"]); w.writerow([r[c] for c in PLAN_COLS] + ([f"{v[0]:.4f}", v[1]] + list(v[2]) if v else ["", "", "", "", ""]))
    log(f"per-mutation table {DDG_TSV} md5 {H.md5_file(DDG_TSV)}")
    losses = {r["key"]: r for r in plan if r["set"] == "loss"}
    ctrl = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in plan:
        if r["set"].startswith("ctrl_"): ctrl[r["set"]][r["loss_key"]].append(r["key"])
    tables = {}
    for run in (1, 2, 3):
        out = []
        for u, desc, cs, kinds, allm in UNITS:
            pairs = []; lv = []; cv = set(); n2 = 0
            for k, L in losses.items():
                if L["role"] not in kinds or (not allm and L["fragment_model"] == "yes") or k not in dd: continue
                cs_ = [c for c in ctrl[cs].get(k, []) if c in dd]
                if not cs_: continue
                pairs.append((L["accession"], dd[k][0] - float(np.mean([dd[c][0] for c in cs_])))); lv.append(dd[k][0]); cv.update(cs_); n2 += len(cs_) == 2
            nprot = len({a for a, d in pairs}); cvv = np.array([dd[c][0] for c in sorted(cv)]); lvv = np.array(lv)
            base = [run, u, desc, cs, len(pairs), nprot, n2]
            if nprot < MIN_PROT:
                out.append(base + ["", "", "", "", "", "Not tested (too few data)"] + [""] * 7); continue
            est, c95, c90 = boot_unit(pairs, [BOOT[run], u], N_BOOT)
            out.append(base + [f"{est:.4f}", f"{c95[0]:.4f}", f"{c95[1]:.4f}", f"{c90[0]:.4f}", f"{c90[1]:.4f}", label(c95, c90),
                               f"{np.median(lvv):.4f}", f"{np.median(cvv):.4f}", f"{np.mean(lvv > DESTAB):.4f}", f"{np.mean(lvv > HIGH):.4f}",
                               f"{np.mean(cvv > DESTAB):.4f}", f"{np.mean(cvv > HIGH):.4f}", len(cvv)])
        tables[run] = out
        with open(RUN_TSV.format(run), "w", newline="", encoding="utf-8") as o:
            w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(RUN_COLS); w.writerows(out)
        log(f"run {run} (bootstrap base {BOOT[run]}): {RUN_TSV.format(run)} md5 {H.md5_file(RUN_TSV.format(run))}")
    cm = lambda t: hashlib.md5("\n".join("\t".join(str(x) for x in row[1:]) for row in t).encode()).hexdigest()
    log(f"run 1 and run 2 identical apart from the run column: {cm(tables[1]) == cm(tables[2])} (content md5 {cm(tables[1])}, {cm(tables[2])}); run 3 content md5 {cm(tables[3])}")
    with open(CMP_TSV, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(RUN_COLS[1:] + ["label_run3", "ci95_lo_run3", "ci95_hi_run3", "ci90_lo_run3", "ci90_hi_run3"])
        for a, b in zip(tables[1], tables[3]): w.writerow(a[1:] + [b[12], b[8], b[9], b[10], b[11]])
    log(f"results table {CMP_TSV} md5 {H.md5_file(CMP_TSV)}")
    log("=== FoldX comparison (labels from run 1; d = ddG(loss) - mean ddG(controls), kcal/mol, > 0 = losses more destabilising) ===")
    for a, b in zip(tables[1], tables[3]):
        if a[12].startswith("Not tested"): log(f"  unit {a[1]} {a[2]}: {a[4]} pairs, {a[5]} proteins: {a[12]}"); continue
        log(f"  unit {a[1]} {a[2]}: {a[4]} pairs in {a[5]} proteins ({a[6]} with 2 controls); median d {a[7]} (95% CI {a[8]} to {a[9]}; 90% CI {a[10]} to {a[11]}) -> {a[12]}"
            f" | run 3: {b[12]}, 95% CI {b[8]} to {b[9]}")
        log(f"      median ddG losses {a[13]}, controls {a[14]} ({a[19]} unique controls); share > 1 kcal/mol: losses {a[15]}, controls {a[17]}; > 3 kcal/mol: losses {a[16]}, controls {a[18]}")
    d3 = sum(1 for a, b in zip(tables[1], tables[3]) if a[12] != b[12]); log(f"units whose label differs in run 3: {d3}")
    log("=== losses reported one by one (D3 A: not in the paired comparison) ===")
    for k, L in sorted(losses.items()):
        if L["role"] == "odd": log(f"  {L['gene']} {L['wt']}{L['position']}{L['mt']} ({L['accession']} F{L['fragment']}, pLDDT {L['plddt']}): ddG {dd[k][0]:.2f}" if k in dd else f"  {k}: no ddG")

# ------------------------------------------------------------------ section 5 ------------------------------------------------------------------
def section5(log):
    plan = check_plan(log)
    if not os.path.exists(CMP_TSV): raise RuntimeError("the loss comparison (section 4) is not written yet: gains run after it")
    S10.check_inputs(log); use, ranges, seqs = S10.load_tables()
    before = set(load_results(["c4", "core"]))
    for st in ("gain_plausible", "gain_unknown"):
        rows = [r for r in plan if r["set"] == st]
        run_jobs(log, st, rows, seqs, ranges, skip_keys=frozenset(before))
        before |= set(load_results([st]))
    dd = load_results(["c4", "core", "gain_plausible", "gain_unknown"])
    gr = [r for r in plan if r["set"] in ("gain_plausible", "gain_unknown")]
    with open(GAIN_DDG, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(PLAN_COLS + ["ddg_mean", "ddg_sd", "run0", "run1", "run2"])
        for r in gr:
            v = dd.get(r["key"]); w.writerow([r[c] for c in PLAN_COLS] + ([f"{v[0]:.4f}", v[1]] + list(v[2]) if v else ["", "", "", "", ""]))
    log(f"gain table {GAIN_DDG} md5 {H.md5_file(GAIN_DDG)}")
    out = []
    for lab, sets in (("plausible", ("gain_plausible",)), ("plausible_or_unknown", ("gain_plausible", "gain_unknown"))):
        for allm in (True, False):
            v = np.array([dd[r["key"]][0] for r in gr if r["set"] in sets and r["key"] in dd and (allm or r["fragment_model"] == "no")])
            n_plan = sum(1 for r in gr if r["set"] in sets and (allm or r["fragment_model"] == "no"))
            row = [lab, "all models" if allm else "without fragment models", n_plan, len(v)] + ([f"{np.median(v):.4f}", f"{np.percentile(v, 25):.4f}", f"{np.percentile(v, 75):.4f}",
                   f"{np.mean(v > DESTAB):.4f}", f"{np.mean(v > HIGH):.4f}"] if len(v) else [""] * 5)
            out.append(row); log(f"  gains {lab}, {row[1]}: {len(v)} of {n_plan} with a ddG; median {row[4]} (IQR {row[5]} to {row[6]}); share > 1 kcal/mol {row[7]}; > 3 kcal/mol {row[8]}")
    with open(GAIN_SUM, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["gain_set", "models", "n_planned", "n_with_ddg", "median_ddg", "q25", "q75", "share_gt1", "share_gt3"]); w.writerows(out)
    log(f"gain summary {GAIN_SUM} md5 {H.md5_file(GAIN_SUM)} (descriptive, no test)")

SECTIONS = {"1": ("plan", section1), "2": ("control_C4", section2), "3": ("core_runs", section3), "4": ("comparison", section4), "5": ("gains", section5)}

def run_one(k):
    name, fn = SECTIONS[k]; ok = {"v": False}
    def wrapped(log):
        fn(log); ok["v"] = True
    H.run_section("P2_11", k, name, wrapped)
    return ok["v"]

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "1"
    for k in (list(SECTIONS) if which == "all" else [which]):
        if not run_one(k): print(f"SECTION {k} DID NOT FINISH: the chain stops here (see the log)", flush=True); break
# ---- END OF SCRIPT ----
