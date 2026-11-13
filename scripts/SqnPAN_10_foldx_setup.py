# SqnPAN_10_foldx_setup.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# STRUCTURE CHECK AND TIMING for the FoldX part of SqnPAN (analysis plan sections 8 and 10 with Amendment 3).
# No FoldX result of interest is computed here: no T1 loss substitution and no control mutation is built. The FoldX run itself (control C4, loss substitutions,
#   matched controls, gains, comparison) is the next script, SqnPAN_11, with its own pre-specification fixed after this script's timing was read.
# Section 1 (format test, a few minutes): one small, confident AlphaFold model of a protein with no T1 loss mutation (single model, 100 to 160 residues, highest mean
#   pLDDT among those with a buried hydrophobic residue and an exposed lysine; ties to the lower accession). Tests:
#   A RepairPDB on the PDB written from the packed tar: return code 0, output <base>_Repair.pdb exists [name FROM MEMORY], same residues and numbering as the input.
#   B BuildModel with the mutation list (--mutant-file=individual_list.txt, line format WT+chain+number+mutant+";" [FROM MEMORY]), 3 mutations x numberOfRuns 3, mutant PDBs
#     written: Dif table has 3 x 3 rows named <base>_<mutation>_<run>.pdb [FROM MEMORY], Average table 3 rows, mean of the Dif runs equals the Average value within 0.01,
#     and the mutant PDB of run 0 carries the mutant residue at the given number (proves the numbering).
#   C Same command as B in a fresh folder, mutant PDBs deleted after the call: identical ddG text in every row (two separate calls give the same values) and no
#     mutant PDB left in the folder.
#   D Same as B without any physics argument: identical ddG text to B (proves that the explicit values 298 K, pH 7, ionic strength 0.05, -IGNORE, vdwDesign 2,
#     moveNeighbours true equal the program defaults; decision D3 A).
#   POST HOC AMENDMENT (2026-10-25, option B; made after the first two runs of section 1, before any FoldX result of interest). The first version passed
#     --out-pdb=0 in C and D and required C to equal B (--out-pdb=1) to the last digit. Run 2 gave FAILED on that check only: one mutation differed in the 6th
#     significant digit (largest 0.000002 kcal/mol); C and D (two separate calls) were identical, and FoldX refused a wrong wild-type residue (E). Decision: BuildModel runs
#     at the program default for --out-pdb (the argument is not passed, mutant PDBs are written), the PDBs are deleted after each call, and C and D are redefined as above.
#     Run 1 stopped on a parser error: the Average table names its rows <base>_<mutation> without .pdb (fixed in the helpers).
#   E Mutation list with a wrong wild-type letter (our own wild-type check switched off for this call only): FoldX's behaviour is reported, not judged. Our code always
#     checks the wild type itself before writing a list (P2_10_foldx_helpers.build).
#   PASSED if A to D all pass; written to data\processed\P2_10_setup_verdict.txt. The three test mutations (buried hydrophobic to A, exposed residue to A, exposed K to R) are a
#   format test on one protein, not part of control C4; that protein is left out of the C4 draw in SqnPAN_11.
# Section 2 (timing, runs until the largest model is done; can take hours): refuses to run unless section 1 PASSED. Builds the FoldX set as defined in analysis plan section 8 (unique
#   T1 loss substitutions in chosen files, protein usable by rule 1A, pLDDT >= 70 at the position; fragment by rule 2A) and prints its counts (no FoldX value). Timing proteins
#   from that set: the single-model proteins at the 10%, 50%, 75% and 90% points of model length and the largest single model, plus the fragment-model protein with the most
#   eligible losses (its fragment holding most of them). Each gets RepairPDB and BuildModel with 3 timing mutations (random residues with pLDDT >= 70 that are not T1 loss
#   positions, not G, P or A, not terminal, changed to A; seed SeedSequence(31800000 + 2)); ddG values of timing mutations are not printed. All six run at the same time
#   (threads, one FoldX process each), so the times include that load. Repaired models are kept in data\processed\P2_10_foldx\repaired for reuse in SqnPAN_11.
#   Projection (a guess from six timed proteins): log-log fit of seconds on model length for RepairPDB and for BuildModel per mutation, applied to every
#   model and mutation count of the planned run (losses x 3 for the loss plus two matched controls; gains in T1 glycoproteins with pLDDT >= 70, plausible and plausible
#   plus unknown; 100 C4 mutations at the median model length), printed as CPU hours and as wall hours for 4, 6, 8 and 12 parallel FoldX processes, against the 5-day rule
#   of analysis plan section 8 and the build expiry 2026-12-31.
# Seeds: base 31800000 (section 1 stream 31800001, section 2 stream 31800002); at least 100,000 apart from 31000000 to 31700000 (checked).
# Usage (from the repository folder):  python scripts/SqnPAN_10_foldx_setup.py 1      (section 1)
#                          python scripts/SqnPAN_10_foldx_setup.py 2      (section 2; or "all" for both)
# Logs: results\logs\P2_10_section1_format_test.txt and P2_10_section2_timing.txt. Outputs: data\processed\P2_10_setup_verdict.txt,
#   data\processed\P2_10_foldx\ (working folders and repaired models), results\tables\P2_10_timing.tsv.
import os, sys, csv, gzip, glob, math, time, shutil, datetime, collections, threading
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_04_alphafold_helpers as AF
import SqnPAN_10_foldx_helpers as FX

SEED_BASE = 31800000
USED_BASES = [31000000, 31100000, 31200000, 31300000, 31400000, 31500000, 31600000, 31700000]
INPUTS = {FX.FOLDX_EXE: FX.FOLDX_MD5,
          os.path.join(H.PROC, "P2_04_alphafold_index.tsv"): "b674fa1d5c67da8bbfa5145f30ebadf1",
          os.path.join(H.PROC, "P2_04_alphafold_use.tsv"): "bc009789438ac60310524d4ef3e6f901",
          os.path.join(H.PROC, "P2_06_mutations.tsv.gz"): "3be28f9265312957c13ab2f0faa2b2b2",
          os.path.join(H.PROC, "P2_06_proteins.tsv"): "f77f5ffabf9bdce910e4aa71f8099eb4",
          os.path.join(H.PROC, "P2_07_residues.tsv.gz"): "9681c55c6f177dd128ad9ab64fec12fd"}
PRINT_ONLY = [os.path.join(H.PROC, "P2_04_alphafold_models.tsv"), os.path.join(H.PROC, "P2_04_alphafold_sequences.tsv")]
VERDICT = os.path.join(H.PROC_OUT, "P2_10_setup_verdict.txt")
TIMING_TSV = os.path.join(H.TABLES, "P2_10_timing.tsv")
SETUP_DIR = os.path.join(FX.FOLDX_WORK, "setup")
HYDROPHOBIC = set("ILVFM")
_IDX = {}; _IDX_LOCK = threading.Lock()

def check_inputs(log):
    for b in USED_BASES:
        if abs(SEED_BASE - b) < 100000: raise ValueError(f"seed base {SEED_BASE} is closer than 100,000 to {b}")
    log(f"seed base {SEED_BASE}: at least 100,000 from {USED_BASES} (checked)")
    bad = 0
    for p, m in INPUTS.items():
        got = H.md5_file(p) if os.path.exists(p) else "MISSING"
        log(f"  {'OK ' if got == m else 'BAD'} {p} md5 {got}" + ("" if got == m else f" (expected {m})")); bad += got != m
    for p in PRINT_ONLY: log(f"  (printed only) {p} md5 {H.md5_file(p)}")
    if bad: raise ValueError(f"{bad} input file(s) with a wrong md5: stop")
    today = datetime.date.today(); exp = datetime.date.fromisoformat(FX.FOLDX_EXPIRY)
    log(f"today {today}; FoldX build expires {exp}: {(exp - today).days} days left")
    if today >= exp: raise ValueError("the FoldX build has expired: download a new build and record it")
    log(f"python {sys.version.split()[0]} | numpy {np.__version__} | cpu count {os.cpu_count()}")
    try:
        import psutil; log(f"psutil {psutil.__version__}: peak memory of each FoldX process is measured")
    except ImportError:
        log("psutil not installed: peak memory of FoldX is not measured (Task Manager can be watched instead)")

def load_tables():
    use = {}
    with open(H.USE_TSV, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"): use[r["accession"]] = r
    ranges = collections.defaultdict(list)
    with open(os.path.join(H.PROC, "P2_04_alphafold_models.tsv"), newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"): ranges[r["accession"]].append((int(r["fragment"]), int(r["align_begin"]), int(r["align_end"])))
    seqs = {}
    with open(os.path.join(H.PROC, "P2_04_alphafold_sequences.tsv"), newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"): seqs[r["accession"]] = r["sequence"]
    return use, ranges, seqs

def load_mutations():
    """Unique substitutions of the chosen files: T1 losses {(acc, pos, ref, alt): res_role} and gains {(acc, pos, ref, alt): gain_state_best}."""
    loss = {}; gain = {}
    with gzip.open(os.path.join(H.PROC, "P2_06_mutations.tsv.gz"), "rt", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["in_chosen_file"] != "1": continue
            k = (r["accession"], int(r["position"]), r["ref"], r["alt"])
            if "T1" in r["lost_tiers"].split(";"): loss[k] = r["res_role"]
            if r["class"] in ("gain", "loss_and_gain"): gain[k] = r["gain_state_best"]
    return loss, gain

def scan_residues(accs):
    """{acc: {position: (aa, plddt, rsa or None, fragment, terminal)}} for the given accessions."""
    out = collections.defaultdict(dict)
    with gzip.open(os.path.join(H.PROC, "P2_07_residues.tsv.gz"), "rt", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            a = r["accession"]
            if a in accs:
                out[a][int(r["position"])] = (r["aa"], float(r["plddt"]), float(r["rsa"]) if r["rsa"] else None, int(r["fragment"]), int(r["terminal"]))
    return out

def write_model(acc, fr, seqs, ranges, folder, log=None):
    """Writes <acc>_F<fr>.pdb from the packed tar into folder and checks it against UniProt 2026_03. Returns (path, begin, end, n_res, numbering_ok, record types)."""
    with _IDX_LOCK:
        if "idx" not in _IDX: _IDX["idx"] = AF.load_index()
    idx = _IDX["idx"]; fh = AF.open_archive()
    try:
        text = FX.model_pdb_text(fh, idx, acc, fr)
    finally:
        fh.close()
    b, e = [(bb, ee) for f2, bb, ee in ranges[acc] if f2 == fr][0]
    n, mism, chains, num_ok = FX.check_model_sequence(text, seqs[acc], b)
    if mism or n != e - b + 1 or chains != ["A"] or not num_ok:
        raise ValueError(f"{acc} F{fr}: model check failed (residues {n}, expected {e - b + 1}, mismatches {mism}, chains {chains}, numbering 1..n {num_ok})")
    os.makedirs(folder, exist_ok=True); path = os.path.join(folder, f"{acc}_F{fr}.pdb")
    with open(path, "w", newline="\n") as o: o.write(text)
    return path, b, e, n, num_ok, FX.record_types(text)

# ------------------------------------------------------------------ section 1 ------------------------------------------------------------------
def section1(log):
    check_inputs(log)
    use, ranges, seqs = load_tables()
    loss, gain = load_mutations()
    loss_prot = {k[0] for k in loss}
    cand = {a for a, r in use.items() if r["status"] == "use" and r["n_fragments"] == "1" and 100 <= int(r["alphafold_length"]) <= 160 and a not in loss_prot}
    log(f"candidate test proteins (usable, single model, 100-160 residues, no T1 loss substitution): {len(cand)}")
    feat = scan_residues(cand); best = None
    for a in sorted(feat):
        v = feat[a]
        bur = [p for p, x in v.items() if x[0] in HYDROPHOBIC and x[2] is not None and x[2] <= 0.05 and x[1] >= 90 and not x[4]]
        lys = [p for p, x in v.items() if x[0] == "K" and x[2] is not None and x[2] >= 0.5 and x[1] >= 70 and not x[4]]
        oth = [p for p, x in v.items() if x[0] not in "GPAK" and x[2] is not None and x[2] >= 0.5 and x[1] >= 70 and not x[4]]
        if not (bur and lys and oth): continue
        m = float(np.mean([x[1] for x in v.values()]))
        if best is None or m > best[0]: best = (m, a, sorted(bur), sorted(lys), sorted(oth))
    m, acc, bur, lys, oth = best
    log(f"test protein {acc} ({use[acc]['gene']}), {use[acc]['alphafold_length']} residues, mean pLDDT {m:.2f}")
    rng = np.random.default_rng(np.random.SeedSequence(SEED_BASE + 1))
    p1 = int(rng.choice(bur)); p2 = int(rng.choice([p for p in oth if p != p1])); p3 = int(rng.choice([p for p in lys if p not in (p1, p2)]))
    v = feat[acc]
    muts = [(v[p1][0], "A", p1, "A"), (v[p2][0], "A", p2, "A"), ("K", "A", p3, "R")]
    log("test mutations (model number = UniProt position for this single model): " + ", ".join(
        f"{FX.mutation_line(*mu)} (pLDDT {v[mu[2]][1]:.1f}, rsa {v[mu[2]][2]:.3f})" for mu in muts))
    if os.path.isdir(os.path.join(SETUP_DIR, "s1")): shutil.rmtree(os.path.join(SETUP_DIR, "s1"))
    d0 = os.path.join(SETUP_DIR, "s1", "model")
    path, b, e, n, num_ok, rec = write_model(acc, 1, seqs, ranges, d0)
    log(f"PDB written {path}: {n} residues, chain A, numbered 1..{n}, all residues equal UniProt 2026_03; record types {rec}")
    log("first 6 lines of the PDB: "); [log("   | " + l) for l in FX.file_head(path, 6)]
    checks = {}
    # A RepairPDB
    dA = os.path.join(SETUP_DIR, "s1", "A_repair"); rA = FX.repair(path, dA)
    log(f"A RepairPDB: rc {rA['rc']}, {rA['seconds']:.1f} s, peak MB {rA['peak_mb']}, output {os.path.basename(rA['repaired'])} exists {rA['exists']}")
    log(f"  args: {rA['args']}"); log("  files in the folder: " + ", ".join(sorted(os.listdir(dA))))
    log("  last stdout lines:"); [log("   | " + l) for l in FX.stdout_tail(rA["stdout"], 15)]
    with open(rA["stdout"], encoding="utf-8", errors="replace") as f:
        flagged = [l.strip() for l in f if any(w in l.lower() for w in ("rotabase", "error", "warning", "not found", "cannot"))]
    log(f"  stdout lines mentioning rotabase/error/warning/not found/cannot: {len(flagged)}"); [log("   ! " + l) for l in flagged[:20]]
    okA = rA["rc"] == 0 and rA["exists"]
    if okA:
        with open(path) as f: r_in = [(c, nn, nm) for c, nn, ic, nm in FX.pdb_ca_residues(f.read())]
        with open(rA["repaired"]) as f: r_out = [(c, nn, nm) for c, nn, ic, nm in FX.pdb_ca_residues(f.read())]
        okA = r_in == r_out
        log(f"  repaired model has the same residues and numbering as the input: {okA}")
    checks["A RepairPDB runs, output name, same residues"] = okA
    if not okA:
        finish1(log, checks); return
    rep = rA["repaired"]
    # B BuildModel with mutant PDBs
    dB = os.path.join(SETUP_DIR, "s1", "B_build_pdb"); rB = FX.build(rep, muts, dB, tag="BuildModel_B", keep_pdbs=True)
    log(f"B BuildModel (--out-pdb not passed = default): rc {rB['rc']}, {rB['seconds']:.1f} s, peak MB {rB['peak_mb']}, mutant and WT PDBs written {rB['pdbs_written']}")
    log(f"  args: {rB['args']}")
    log("  mutation list written: " + " | ".join(FX.file_head(os.path.join(dB, FX.MUTANT_FILE), 5)))
    log("  files in the folder: " + ", ".join(sorted(os.listdir(dB))))
    log("  last stdout lines:"); [log("   | " + l) for l in FX.stdout_tail(rB["stdout"], 12)]
    okB = False; difB = None
    if rB["rc"] == 0 and os.path.exists(rB["dif"]) and os.path.exists(rB["average"]):
        for p in (rB["dif"], rB["average"]):
            log(f"  first lines of {os.path.basename(p)}:"); [log("   | " + l) for l in FX.file_head(p, 14)]
        difB = FX.dif_ddg(rB["dif"]); avgB = FX.average_ddg(rB["average"])
        shape = sorted(difB) == [1, 2, 3] and all(sorted(difB[i]) == [0, 1, 2] for i in difB) and sorted(avgB) == [1, 2, 3]
        meanok = shape and all(abs(np.mean([float(x) for x in difB[i].values()]) - float(avgB[i][1])) <= 0.01 for i in difB)
        pdbok = True
        for i, (wt, ch, num, mt) in enumerate(muts, start=1):
            mp = os.path.join(dB, f"{rB['base']}_{i}_0.pdb")
            if not os.path.exists(mp): pdbok = False; log(f"  mutant PDB {os.path.basename(mp)} not found"); continue
            with open(mp) as f: got = {(c, nn): nm for c, nn, ic, nm in FX.pdb_ca_residues(f.read())}.get((ch, num))
            log(f"  mutant PDB {os.path.basename(mp)}: residue at {ch}{num} = {got} (expected {FX.ONE[mt]})"); pdbok &= got == FX.ONE[mt]
        log(f"  Dif rows 3 mutations x 3 runs and Average 3 rows: {shape}; Dif mean = Average within 0.01: {meanok}; mutant residue at the given number: {pdbok}")
        log("  ddG of the three format-test mutations (kcal/mol; one protein, format test only, not control C4): " + ", ".join(
            f"{FX.mutation_line(*muts[i - 1])} runs {[difB[i][r] for r in (0, 1, 2)]} mean {avgB[i][1]} SD {avgB[i][0]}" for i in (1, 2, 3)))
        okB = shape and meanok and pdbok
    checks["B BuildModel with --mutant-file, output names, numbering"] = okB
    if not okB:
        finish1(log, checks); return
    # C same command in a fresh folder, PDBs deleted after the call; D no physics arguments (amended 2026-10-25, option B)
    dC = os.path.join(SETUP_DIR, "s1", "C_build_repeat"); rC = FX.build(rep, muts, dC, tag="BuildModel_C")
    dD = os.path.join(SETUP_DIR, "s1", "D_build_defaults"); rD = FX.build(rep, muts, dD, phys=False, tag="BuildModel_D")
    for lab, r, d in (("C", rC, dC), ("D", rD, dD)):
        log(f"{lab} BuildModel: rc {r['rc']}, {r['seconds']:.1f} s, mutant and WT PDBs written {r['pdbs_written']} (deleted after the call); args: {r['args']}")
        log("  files left in the folder: " + ", ".join(sorted(os.listdir(d))))
    difC = FX.dif_ddg(rC["dif"]) if rC["rc"] == 0 and os.path.exists(rC["dif"]) else None
    difD = FX.dif_ddg(rD["dif"]) if rD["rc"] == 0 and os.path.exists(rD["dif"]) else None
    nopdb = rC["pdbs_written"] > 0 and not FX.mutant_pdbs(dC, rC["base"])
    checks["C same command again: identical ddG to B, mutant PDBs deleted"] = difC == difB and nopdb
    checks["D no physics arguments: identical ddG to B (explicit values = defaults)"] = difD is not None and difD == difB
    log(f"C identical to B in every row: {difC == difB}; PDBs were written and none is left after deletion: {nopdb}")
    log(f"D identical to B in every row: {difD == difB}")
    # E wrong wild type
    wt0, ch0, num0, mt0 = muts[0]; wrong = "G" if wt0 != "G" else "A"
    dE = os.path.join(SETUP_DIR, "s1", "E_wrong_wt"); rE = FX.build(rep, [(wrong, ch0, num0, mt0)], dE, tag="BuildModel_E", check_wt=False)
    rowsE = None
    if os.path.exists(rE["dif"]):
        try: rowsE = FX.dif_ddg(rE["dif"])
        except Exception as ex: rowsE = f"unreadable ({ex})"
    log(f"E wrong wild type ({FX.mutation_line(wrong, ch0, num0, mt0)}, the PDB has {wt0}): rc {rE['rc']}, Dif file exists {os.path.exists(rE['dif'])}, rows {rowsE}")
    log("  last stdout lines:"); [log("   | " + l) for l in FX.stdout_tail(rE["stdout"], 8)]
    log("  (reported, not judged: our code checks the wild type itself before every list)")
    finish1(log, checks)

def finish1(log, checks):
    ok = all(checks.values()) and len(checks) == 4
    log("=== section 1 checks ===")
    for k, v in checks.items(): log(f"  {'PASS' if v else 'FAIL'}  {k}")
    if len(checks) < 4: log("  (later checks not reached)")
    os.makedirs(os.path.dirname(VERDICT), exist_ok=True)
    with open(VERDICT, "w", newline="\n") as o:
        o.write(("PASSED" if ok else "FAILED") + f" {time.strftime('%Y-%m-%d %H:%M:%S')} P2_10 section 1 format test\n")
        for k, v in checks.items(): o.write(f"{'PASS' if v else 'FAIL'}\t{k}\n")
    log(f"VERDICT: {'PASSED' if ok else 'FAILED'} (written to {VERDICT})")
    if not ok: raise RuntimeError("section 1 FAILED: section 2 will not run")

# ------------------------------------------------------------------ section 2 ------------------------------------------------------------------
def section2(log):
    if not (os.path.exists(VERDICT) and open(VERDICT).read().startswith("PASSED")): raise RuntimeError("section 1 has not PASSED: run section 1 first")
    check_inputs(log)
    use, ranges, seqs = load_tables()
    loss, gain = load_mutations()
    nT1 = {}
    with open(os.path.join(H.PROC, "P2_06_proteins.tsv"), newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"): nT1[r["accession"]] = int(r["n_T1_sequons"])
    usable = {a for a, r in use.items() if r["status"] == "use"}
    accs = {k[0] for k in loss} | {k[0] for k in gain if nT1.get(k[0], 0) > 0}
    feat = scan_residues(accs & usable)
    # FoldX set (analysis plan section 8): counts only
    cnt = collections.Counter(); elig = []
    for k, role in sorted(loss.items()):
        a, p = k[0], k[1]
        if a not in usable: cnt["protein not usable (rule 1A or no model)"] += 1; continue
        x = feat[a].get(p)
        if x is None or x[0] != k[2]: raise ValueError(f"residue table and mutation table disagree at {a} {p}")
        if x[1] < 70: cnt["pLDDT < 70"] += 1; continue
        elig.append((a, p, k[2], k[3], role, x[3]))
    log(f"FoldX set: unique T1 loss substitutions in chosen files {len(loss)}; left out: {dict(cnt)}; eligible {len(elig)}")
    log(f"  by position role: {dict(collections.Counter(e[4] for e in elig))}")
    prot = collections.Counter(e[0] for e in elig); models = collections.Counter((e[0], e[5]) for e in elig)
    fragp = {a for a in prot if use[a]["n_fragments"] != "1"}
    log(f"  proteins {len(prot)}; models to repair (protein x fragment) {len(models)}; fragment-model proteins {len(fragp)} with {sum(prot[a] for a in fragp)} losses")
    def mlen(a, fr): return [ee - bb + 1 for f2, bb, ee in ranges[a] if f2 == fr][0]
    gcnt = collections.Counter(); gmodels = {"plausible": collections.Counter(), "plausible_or_unknown": collections.Counter()}
    for k, st in gain.items():
        a, p = k[0], k[1]
        if nT1.get(a, 0) == 0 or a not in usable: continue
        x = feat[a].get(p)
        if x is None or x[1] < 70: continue
        gcnt[st] += 1
        if st == "plausible": gmodels["plausible"][(a, x[3])] += 1
        if st in ("plausible", "unknown"): gmodels["plausible_or_unknown"][(a, x[3])] += 1
    log(f"  gains in T1 glycoproteins with pLDDT >= 70 (unique substitutions, descriptive set): {dict(gcnt)}")
    for g, mc in gmodels.items(): log(f"    {g}: {sum(mc.values())} mutations in {len(mc)} models, {len(set(mc) - set(models))} of them models not in the loss set")
    # timing proteins
    single = sorted(((mlen(a, 1), a) for a in prot if use[a]["n_fragments"] == "1"))
    picks = []
    for q in (0.10, 0.50, 0.75, 0.90, 1.0):
        L, a = single[int(round(q * (len(single) - 1)))]
        if a not in [x[0] for x in picks]: picks.append((a, 1, f"single model, {int(q * 100)}% point of length"))
    fa = sorted(fragp, key=lambda a: (-prot[a], a))[0]
    ffr = sorted(((-c, fr) for (a, fr), c in models.items() if a == fa))[0][1]
    picks.append((fa, ffr, f"fragment model, most eligible losses ({prot[fa]})"))
    rng = np.random.default_rng(np.random.SeedSequence(SEED_BASE + 2))
    lossp = {(k[0], k[1]) for k in loss}; items = []
    for a, fr, why in picks:
        b, e = [(bb, ee) for f2, bb, ee in ranges[a] if f2 == fr][0]
        pool = sorted(p for p, x in feat[a].items() if b <= p <= e and x[3] == fr and x[1] >= 70 and x[0] not in "GPA" and not x[4] and (a, p) not in lossp)
        pos = sorted(int(p) for p in rng.choice(pool, size=3, replace=False))
        muts = [(feat[a][p][0], "A", p - b + 1, "A") for p in pos]
        items.append({"acc": a, "fr": fr, "why": why, "L": e - b + 1, "begin": b, "muts": muts, "n_loss_in_model": models.get((a, fr), 0)})
        log(f"timing protein {a} ({use[a]['gene']}) F{fr}: model length {e - b + 1}, {why}; timing mutations " + ", ".join(FX.mutation_line(*m) for m in muts))
    os.makedirs(FX.REPAIRED_DIR, exist_ok=True)
    def time_one(it):
        wd = os.path.join(SETUP_DIR, f"s2_{it['acc']}_F{it['fr']}")
        if os.path.isdir(wd): shutil.rmtree(wd)
        path = write_model(it["acc"], it["fr"], seqs, ranges, wd)[0]
        r1 = FX.repair(path, wd); it["repair_s"] = r1["seconds"]; it["repair_mb"] = r1["peak_mb"]; it["repair_rc"] = r1["rc"]
        if not (r1["rc"] == 0 and r1["exists"]): it["error"] = "RepairPDB failed: " + " / ".join(FX.stdout_tail(r1["stdout"], 3)); return it
        shutil.copy2(r1["repaired"], os.path.join(FX.REPAIRED_DIR, os.path.basename(r1["repaired"])))
        r2 = FX.build(r1["repaired"], it["muts"], wd, tag="BuildModel_timing")
        it["build_s"] = r2["seconds"]; it["build_mb"] = r2["peak_mb"]; it["build_rc"] = r2["rc"]
        try: it["rows"] = sum(len(v) for v in FX.dif_ddg(r2["dif"]).values())
        except Exception as ex: it["error"] = f"BuildModel output unreadable: {ex}"
        return it
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=len(items)) as ex:
        fut = {ex.submit(time_one, it): it for it in items}; pending = set(fut)
        while pending:
            done, pending = wait(pending, timeout=300, return_when=FIRST_COMPLETED)
            for f in done:
                it = f.result()
                log(f"  finished {it['acc']} F{it['fr']} (length {it['L']}) after {(time.time() - t0) / 60:.1f} min: repair {it.get('repair_s', float('nan')):.0f} s, "
                    f"build 3 mutations {it.get('build_s', float('nan')):.0f} s, rows {it.get('rows')}, peak MB {it.get('repair_mb')}/{it.get('build_mb')}"
                    + (f", ERROR {it['error']}" if "error" in it else ""))
            if pending: log(f"  ... {(time.time() - t0) / 60:.1f} min, still running: " + ", ".join(f"{fut[f]['acc']} (length {fut[f]['L']})" for f in pending))
    with open(TIMING_TSV, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "fragment", "model_length", "why", "repair_seconds", "repair_peak_mb", "build_seconds_3_mutations", "build_peak_mb", "build_rows", "error"])
        for it in items: w.writerow([it["acc"], it["fr"], it["L"], it["why"], f"{it.get('repair_s', float('nan')):.1f}", it.get("repair_mb"), f"{it.get('build_s', float('nan')):.1f}", it.get("build_mb"), it.get("rows"), it.get("error", "")])
    log(f"timing table {TIMING_TSV} md5 {H.md5_file(TIMING_TSV)}")
    good = [it for it in items if "error" not in it and it.get("rows") == 9]
    if len(good) < 3: raise RuntimeError(f"only {len(good)} timing proteins finished without error: no projection")
    Ls = np.log([it["L"] for it in good])
    fr_ = np.polyfit(Ls, np.log([it["repair_s"] for it in good]), 1); fb_ = np.polyfit(Ls, np.log([it["build_s"] / 3 for it in good]), 1)
    rep_s = lambda L: math.exp(fr_[1] + fr_[0] * math.log(L)); mut_s = lambda L: math.exp(fb_[1] + fb_[0] * math.log(L))
    log(f"fit (log seconds on log model length): RepairPDB exponent {fr_[0]:.2f}; BuildModel per mutation (3 runs) exponent {fb_[0]:.2f}")
    for it in good: log(f"  {it['acc']} L {it['L']}: repair measured {it['repair_s']:.0f} s, fitted {rep_s(it['L']):.0f} s; per mutation measured {it['build_s'] / 3:.0f} s, fitted {mut_s(it['L']):.0f} s")
    medL = float(np.median([mlen(a, fr) for (a, fr) in models]))
    jobs = {m: rep_s(mlen(*m)) + 3 * c * mut_s(mlen(*m)) for m, c in models.items()}
    core = sum(jobs.values()) + 100 * mut_s(medL)
    extra = {}
    for g, mc in gmodels.items():
        extra[g] = sum(c * mut_s(mlen(*m)) for m, c in mc.items()) + sum(rep_s(mlen(*m)) for m in set(mc) - set(models))
    longest = max(jobs.values())
    log("=== PROJECTION (a guess from the timed proteins, not a measurement of the full run) ===")
    log(f"  losses + 2 matched controls each ({3 * len(elig)} mutations) + {len(models)} repairs + 100 C4 mutations: {core / 3600:.1f} CPU hours")
    for g in extra: log(f"  gains {g}: + {extra[g] / 3600:.1f} CPU hours")
    log(f"  longest single model (repair + its mutations): {longest / 3600:.1f} h")
    days_left = (datetime.date.fromisoformat(FX.FOLDX_EXPIRY) - datetime.date.today()).days
    for W in (4, 6, 8, 12):
        for lab, tot in (("core", core), ("core + plausible gains", core + extra["plausible"]), ("core + plausible or unknown gains", core + extra["plausible_or_unknown"])):
            wall = max(tot / W, longest)
            log(f"  {W:2d} parallel FoldX processes, {lab}: about {wall / 3600:.1f} h wall ({wall / 86400:.2f} days; 5-day rule {'met' if wall <= 5 * 86400 else 'NOT met'})")
    log(f"  days left before the build expires: {days_left}. Times were measured with {len(items)} FoldX processes at once on {os.cpu_count()} logical CPUs; more processes than physical cores slow each one down.")

SECTIONS = {"1": ("format_test", section1), "2": ("timing", section2)}

def run_one(k):
    name, fn = SECTIONS[k]; ok = {"v": False}
    def wrapped(log):
        fn(log); ok["v"] = True
    H.run_section("P2_10", k, name, wrapped)
    return ok["v"]

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "1"
    for k in (list(SECTIONS) if which == "all" else [which]):
        if not run_one(k): print(f"SECTION {k} DID NOT FINISH: the chain stops here (see the log)", flush=True); break
# ---- END OF SCRIPT ----
