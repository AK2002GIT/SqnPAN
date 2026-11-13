# SqnPAN_07_structure_features.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-10-24; decisions 1A and 2A; the analysis plan section 7 with Amendments 1 and 2)
# Purpose: one table with pLDDT and relative solvent exposure for every residue of every usable AlphaFold v6 model (20,249 proteins, 23,275 models, about 11.3 million
#   residues). It feeds the stricter Q2 model, Q3, the FoldX set, control C4 and script 06 of the companion LUAD study. This script joins nothing to the mutation tables.
# Method: Shrake-Rupley, probe 1.4 A, 401 Fibonacci points per atom (DSSP's formula), DSSP atom radii (N 1.65, CA 1.87, C 1.76, O 1.40, other atoms 1.80; read from the
#   PDB-REDO/dssp source on 2026-10-24; Tien et al. used the 2011 DSSP version, not checked). Atom area = 4 pi (r+1.4)^2 x exposed points / 401; residue area = sum over its atoms.
#   Relative exposure rsa = residue area / Tien et al. 2013 theoretical maximum (Table 1, ALLOWED region; PLoS One 8(11):e80635, doi 10.1371/journal.pone.0080635, PMID 24278298).
#   Values above 1 are kept, not capped. pLDDT = B-factor of the CA atom. Heavy atoms only (the models have no hydrogens).
#   Every fragment model is computed on its own; each UniProt position takes its values from the fragment chosen by rule 2A (choose_fragment of SqnPAN_04).
# Output columns of P2_07_residues.tsv.gz (one row per UniProt position): accession, position, aa, plddt, asa, rsa, fragment, n_fragments, terminal (1 = position 1 or the last
#   position of the protein; Tien excludes chain-terminating residues), cut_distance (residues to the nearest artificial fragment cut in the chosen fragment; blank for unsplit
#   proteins; a protein terminus is not a cut). No threshold is applied here; later scripts choose their own. rsa is blank for residues without a Tien value (X, U).
# Controls (section 1; section 2 refuses to start unless the verdict file says PASSED):
#   C-A exact geometry. One atom (r 1.87): area = 4 pi (3.27)^2 to 1e-9 relative. Two equal atoms (r 1.8) at 3.0 A: analytic area per atom 4 pi R^2 - 2 pi R (R - d/2), R = 3.2;
#       error <= 0.5% with 20,001 points and <= 5% with 401 points. One atom with six neighbours at 1.5 A along +-x, +-y, +-z: area exactly 0.
#   C-B neighbour search. On 5 random models (30 to 150 residues) the grid gives the same neighbour pairs and the same exposed-point counts as an all-pairs brute force (difference 0).
#   C-C file consistency. On 200 random usable models: CA B-factors of the cif equal those of the pdb (within 0.01); an independent parser (SqnPAN_04 parse_cif) agrees with ours on residue
#       names, pLDDT and UniProt numbers (= align_begin + label_seq_id - 1); residue letters equal the UniProt sequence. Zero mismatches.
#   C-D scale and biophysics. On 200 random usable proteins, non-terminal residues with pLDDT >= 70: median rsa of each of L, I, V, F < median rsa of each of K, E, D, R (all 16 pairs);
#       share of rsa > 1 <= 1%; 99.9th percentile of rsa <= 1.15 for every residue type (the 1.15 allows a 15% scale mismatch from the DSSP version and the point count; fixed
#       before any value was seen). If C-D fails the script stops and a dated amendment is needed (empirical normalisation); the label is then "exposure on a different scale".
#   C-E point count. On 20 random proteins, rsa from 401 points against rsa from 2,001 points: 95th percentile of |difference| <= 0.03 and <= 2% of residues change side of 0.25.
#   Seeds: bases 31200000 (run twice, results must be identical) and 31300000 (once, must also pass); derived seeds are base + repeat number; all bases at least 100,000 apart and
#       from the SqnPAN_06 seeds (31000000, 31100000), checked in the code.
#   Timing (after the controls): one protein first, then 50 random proteins (seed 31200001); the reported time is this measurement, the projection is labelled as one.
# Section 3 checks: rows per protein = UniProt length (all 20,249 proteins, positions consecutive); every sequon of P2_06_sites.tsv in a usable protein has N at the N position and S/T
#   at +2 in the table (zero mismatches); 50 proteins (seed 31200002) are recomputed and identical to the table; residue letters equal the UniProt sequence. The log prints pooled
#   numbers only (rows, mean pLDDT, share of rsa > 1 over all residues): nothing by tier, gene, mutation or cancer type.
# Limits to state in the Methods: monomer exposure (interfaces and complexes missed); exposure in low-pLDDT regions is not meaningful; radii from DSSP 4, not the 2011 version.
# Usage (from the repository folder):  python scripts/SqnPAN_07_structure_features.py 1            (controls and timing)
#                          python scripts/SqnPAN_07_structure_features.py 2 [workers=N]   (table; a restart continues)
#                          python scripts/SqnPAN_07_structure_features.py 3            (checks)
#   Options for quick test runs only: limit=N (first N proteins), test=1 (small control samples).
# Logs: results\logs\P2_07_section<k>_<name>.txt. Outputs in data\processed: P2_07_controls_verdict.txt, P2_07_parts\ (part files), P2_07_residues.tsv.gz,
#   P2_07_table_checks_verdict.txt.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")   # one BLAS thread per process; the workers are the parallelism
import sys, csv, gzip, math, time, random, shutil, tracemalloc, collections
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import SqnPAN_04_alphafold_helpers as AF
import SqnPAN_05_helpers as H
import SqnPAN_07_helpers as X

LIMIT = None; WORKERS = None; TEST = False
for a in sys.argv[2:]:
    if a.startswith("limit="): LIMIT = int(a[6:])
    if a.startswith("workers="): WORKERS = int(a[8:])
    if a.startswith("test="): TEST = a[5:] == "1"
SEEDS = (31200000, 31300000)
OLD_SEEDS = (31000000, 31100000)
INDEX_MD5 = "b674fa1d5c67da8bbfa5145f30ebadf1"; SITES_MD5 = "34e0a114e0e7a6955418d837e9a5a245"
INDEX = os.path.join(H.PROC, "P2_04_alphafold_index.tsv"); SITES = os.path.join(H.PROC, "P2_06_sites.tsv")
VERDICT = os.path.join(H.PROC_OUT, "P2_07_controls_verdict.txt"); CHECKS_VERDICT = os.path.join(H.PROC_OUT, "P2_07_table_checks_verdict.txt")
PARTS = os.path.join(H.PROC_OUT, "P2_07_parts"); TABLE = os.path.join(H.PROC_OUT, "P2_07_residues.tsv.gz")
BATCH = 200
SIZE = {"bf": 5, "cc": 200, "cd": 200, "ce": 20, "time": 50} if not TEST else {"bf": 2, "cc": 10, "cd": 12, "ce": 3, "time": 5}
HYD, CHG = "LIVF", "KEDR"

def check_inputs(log):
    allseeds = list(SEEDS) + list(OLD_SEEDS)
    for i in range(len(allseeds)):
        for j in range(i + 1, len(allseeds)): assert abs(allseeds[i] - allseeds[j]) >= 100000, "seed bases closer than 100,000"
    log(f"python {sys.version.split()[0]} | numpy {np.__version__} | cpu count {os.cpu_count()} | seed bases {SEEDS} (at least 100,000 apart from each other and from {OLD_SEEDS})")
    for path, want in ((INDEX, INDEX_MD5), (SITES, SITES_MD5)):
        got = H.md5_file(path); log(f"md5 {os.path.basename(path)}: {got} ({'matches' if got == want else 'DIFFERS from ' + want})")
        if got != want: raise RuntimeError("input md5 differs from the recorded value: stop")
    log(f"md5 P2_04_alphafold_use.tsv (recorded now): {H.md5_file(H.USE_TSV)}")
    log("Tien maximum areas (A^2) as used, to compare with Table 1 of the paper: " + " ".join(f"{k}{int(v)}" for k, v in sorted(X.TIEN_MAX.items())) + f" | sum {int(sum(X.TIEN_MAX.values()))}")
    log(f"DSSP radii: {X.RADIUS_BY_ATOM} side atoms {X.RADIUS_SIDE} | probe {X.PROBE} | points {X.N_POINTS}")

# ======================= Section 1: controls and timing =======================
def control_A(log):
    ok = True; U401 = X.unit_points(401); U20k = X.unit_points(20001)
    R = np.array([1.87 + 1.4]); a = X.atom_areas(R, X.exposed_points(np.zeros((1, 3)), R, U401), 401)[0]; want = 4 * math.pi * 3.27 ** 2
    e1 = abs(a / want - 1); okA1 = e1 <= 1e-9; log(f"  C-A one atom: relative error {e1:.2e} (limit 1e-9) {'ok' if okA1 else 'FAIL'}")
    R = np.array([3.2, 3.2]); xyz = np.array([[0, 0, 0], [3.0, 0, 0]], float); want = 4 * math.pi * 3.2 ** 2 - 2 * math.pi * 3.2 * (3.2 - 1.5)
    res = []
    for U, tol in ((U20k, 0.005), (U401, 0.05)):
        err = max(abs(v / want - 1) for v in X.atom_areas(R, X.exposed_points(xyz, R, U), len(U))); res.append(err <= tol)
        log(f"  C-A two atoms, {len(U)} points: largest relative error {err:.4f} (limit {tol}) {'ok' if err <= tol else 'FAIL'}")
    xyz = np.array([[0, 0, 0], [1.5, 0, 0], [-1.5, 0, 0], [0, 1.5, 0], [0, -1.5, 0], [0, 0, 1.5], [0, 0, -1.5]], float); R = np.full(7, 3.2)
    c = int(X.exposed_points(xyz, R, U401)[0]); log(f"  C-A buried atom: exposed points {c} (must be 0) {'ok' if c == 0 else 'FAIL'}")
    return okA1 and all(res) and c == 0

def control_B(ctx, rng, log):
    cand = sorted(a for a, (g, L, nf) in ctx.use.items() if nf == 1 and 30 <= L <= 150)
    diff = 0; pairs_ok = True; n_atoms = 0
    for acc in rng.sample(cand, SIZE["bf"]):
        a = X.parse_atoms(AF.read_model(ctx.fh, ctx.idx, acc, 1, "cif").split("\n")); R = a["radii"] + X.PROBE; xyz = a["xyz"]; n = len(R); n_atoms += n
        d = np.sqrt(((xyz[:, None] - xyz[None]) ** 2).sum(2))
        brute = {(i, j) for i in range(n) for j in range(n) if i != j and d[i, j] < R[i] + R[j]}
        pairs_ok &= (X.neighbour_pair_set(xyz, R) == brute)
        diff += int(np.abs(X.exposed_points(xyz, R, ctx.U) - X.exposed_points_bruteforce(xyz, R, ctx.U)).sum())
    return {"pairs_equal": bool(pairs_ok), "exposed_difference": diff, "atoms": n_atoms}

def pdb_ca(text):
    out = []
    for l in text.split("\n"):
        if l.startswith("ATOM") and l[12:16].strip() == "CA": out.append((l[17:20], int(l[22:26]), float(l[60:66])))
    return out

def control_C(ctx, rng, log):
    models = sorted((a, f) for a, rs in ctx.ranges.items() for f, b, e in rs)
    bad = collections.Counter(); n_res = 0
    for acc, fr in rng.sample(models, SIZE["cc"]):
        ctext = AF.read_model(ctx.fh, ctx.idx, acc, fr, "cif"); ptext = AF.read_model(ctx.fh, ctx.idx, acc, fr, "pdb")
        mine = X.parse_atoms(ctext.split("\n")); ref = AF.parse_cif(ctext); pdb = pdb_ca(ptext); begin = ref["align_begin"]; seq = ctx.seqs[acc]
        if not (len(ref["ca"]) == mine["n_res"] == len(pdb)): bad["residue count differs between parsers or pdb"] += 1; continue
        for k, (sid, comp, un, ures, b) in enumerate(ref["ca"]):
            n_res += 1
            if sid != k + 1: bad["label_seq_id not 1..n"] += 1
            if comp != mine["comp"][k]: bad["residue name differs (our parser, P2_04 parser)"] += 1
            if abs(b - mine["plddt"][k]) > 1e-9: bad["pLDDT differs (our parser, P2_04 parser)"] += 1
            if un != begin + sid - 1: bad["sifts UniProt number != align_begin + label_seq_id - 1"] += 1
            if X.THREE_TO_ONE.get(comp, "X") != ures: bad["sifts UniProt residue != residue name"] += 1
            if abs(pdb[k][2] - b) > 0.01 or pdb[k][0] != comp: bad["pdb CA B-factor or residue name differs from the cif"] += 1
            if begin + k <= len(seq) and seq[begin + k - 1] != X.THREE_TO_ONE.get(comp, "X"): bad["residue letter differs from the UniProt sequence"] += 1
    return {"residues_checked": n_res, "mismatches": dict(bad)}

def parse_rows(lines):
    for l in lines:
        f = l.split("\t"); yield f[2], float(f[3]), (float(f[5]) if f[5] != "" else None), int(f[8])

def control_D(ctx, rng, log):
    accs = sorted(ctx.use); by = collections.defaultdict(list); n = 0
    for acc in rng.sample(accs, SIZE["cd"]):
        lines, _ = X.protein_rows(acc, ctx)
        for aa, pl, rsa, term in parse_rows(lines):
            if term == 0 and pl >= 70 and rsa is not None: by[aa].append(rsa); n += 1
    med = {aa: round(float(np.median(v)), 6) for aa, v in sorted(by.items())}
    p999 = {aa: round(float(np.percentile(v, 99.9)), 6) for aa, v in sorted(by.items())}
    allv = np.concatenate([np.array(v) for v in by.values()]); share = float((allv > 1.0).mean())
    pairs = [(h, c, med[h] < med[c]) for h in HYD for c in CHG]
    return {"n_residues": n, "types_with_data": len(by), "median": med, "p999": p999, "share_gt1": round(share, 6), "order_pairs_ok": sum(p[2] for p in pairs),
            "order_pairs": len(pairs), "p999_max": max(p999.values()), "min_n_type": min(len(v) for v in by.values())}

def control_E(ctx, rng, log):
    accs = sorted(ctx.use); U2 = X.unit_points(2001); d = []; cross = 0; n = 0
    for acc in rng.sample(accs, SIZE["ce"]):
        l1, _ = X.protein_rows(acc, ctx); l2, _ = X.protein_rows(acc, ctx, U=U2)
        for (a1, p1, r1, t1), (a2, p2, r2, t2) in zip(parse_rows(l1), parse_rows(l2)):
            if r1 is None or r2 is None: continue
            d.append(abs(r1 - r2)); n += 1; cross += int((r1 >= 0.25) != (r2 >= 0.25))
    return {"n_residues": n, "p95_abs_diff": round(float(np.percentile(d, 95)), 6), "share_cross_0.25": round(cross / n, 6)}

def judge(r):
    p = {}
    p["C-B"] = r["B"]["pairs_equal"] and r["B"]["exposed_difference"] == 0
    p["C-C"] = r["C"]["residues_checked"] > 0 and not r["C"]["mismatches"]
    D = r["D"]
    p["C-D order"] = D["order_pairs_ok"] == D["order_pairs"]; p["C-D share>1"] = D["share_gt1"] <= 0.01; p["C-D p99.9"] = D["p999_max"] <= 1.15
    p["C-D all 20 types present"] = D["types_with_data"] == 20
    p["C-E"] = r["E"]["p95_abs_diff"] <= 0.03 and r["E"]["share_cross_0.25"] <= 0.02
    return p

def run_random_controls(ctx, seed, label, log):
    rng = random.Random(seed); t0 = time.time()
    r = {"B": control_B(ctx, rng, log), "C": control_C(ctx, rng, log), "D": control_D(ctx, rng, log), "E": control_E(ctx, rng, log)}
    log(f"--- {label} (seed {seed}, {round(time.time() - t0)} s) ---")
    log(f"  C-B atoms {r['B']['atoms']}: neighbour pairs equal {r['B']['pairs_equal']}, exposed-point difference {r['B']['exposed_difference']}")
    log(f"  C-C residues {r['C']['residues_checked']}: mismatches {r['C']['mismatches'] if r['C']['mismatches'] else 'none'}")
    D = r["D"]; log(f"  C-D residues {D['n_residues']} (pLDDT >= 70, non-terminal), types with data {D['types_with_data']}, smallest type n {D['min_n_type']}")
    log("    median rsa: " + " ".join(f"{k}:{D['median'][k]:.3f}" for k in sorted(D["median"])))
    log("    99.9th percentile rsa: " + " ".join(f"{k}:{D['p999'][k]:.3f}" for k in sorted(D["p999"])))
    log(f"    hydrophobic < charged pairs {D['order_pairs_ok']} of {D['order_pairs']} | share rsa > 1: {100 * D['share_gt1']:.3f}% | largest 99.9th percentile {D['p999_max']:.3f}")
    log(f"  C-E residues {r['E']['n_residues']}: 95th percentile |rsa(401) - rsa(2001)| {r['E']['p95_abs_diff']:.4f} | share crossing 0.25 {100 * r['E']['share_cross_0.25']:.3f}%")
    p = judge(r); log("  " + " | ".join(f"{k} {'ok' if v else 'FAIL'}" for k, v in p.items()))
    return r, all(p.values())

def timing(ctx, log):
    accs = sorted(ctx.use); rng = random.Random(SEEDS[0] + 1); samp = rng.sample(accs, SIZE["time"]); tot_res = sum(ctx.use[a][1] for a in ctx.use)
    t = time.perf_counter(); X.protein_rows(samp[0], ctx); t1 = time.perf_counter() - t
    log(f"  timing, first protein {samp[0]} ({ctx.use[samp[0]][1]} residues, {ctx.use[samp[0]][2]} fragment(s)): {t1:.2f} s")
    t = time.perf_counter(); res = 0; sizes = {}
    for a in samp:
        t_a = time.perf_counter(); X.protein_rows(a, ctx); sizes[a] = (ctx.use[a][1], time.perf_counter() - t_a); res += ctx.use[a][1]
    el = time.perf_counter() - t; per = el / res
    log(f"  timing, {len(samp)} random proteins: {res} residues in {el:.1f} s = {1000 * per:.3f} s per 1,000 residues (one process, includes reading and formatting)")
    big = max(samp, key=lambda a: ctx.use[a][1]); tracemalloc.start(); X.protein_rows(big, ctx); cur, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
    log(f"  peak Python/NumPy allocation for the largest sampled protein {big} ({ctx.use[big][1]} residues): {peak / 1e6:.0f} MB (tracemalloc)")
    w = WORKERS or max(1, min(8, (os.cpu_count() or 2) - 2))
    log(f"  PROJECTION from this sample (not a measurement of the full run): {tot_res} residues x {per:.6f} s = {tot_res * per / 3600:.1f} h in one process; about {tot_res * per / 3600 / w:.1f} h with {w} workers")

def section1(log):
    log("SqnPAN_07 section 1: controls C-A to C-E and timing (pre-specification in the header of the script)")
    check_inputs(log)
    okA = control_A(log); log(f"C-A {'PASSED' if okA else 'FAILED'}")
    ctx = X.Ctx(); runs = {}; allok = okA
    for label, seed in (("seed 31200000 run 1", SEEDS[0]), ("seed 31200000 run 2", SEEDS[0]), ("seed 31300000", SEEDS[1])):
        runs[label], ok = run_random_controls(ctx, seed, label, log); allok &= ok
    same = runs["seed 31200000 run 1"] == runs["seed 31200000 run 2"]; allok &= same
    log(f"same seed twice identical: {same}")
    a, b = runs["seed 31200000 run 1"], runs["seed 31300000"]
    log(f"difference between the two seeds (reported, not judged): share rsa > 1 {100 * a['D']['share_gt1']:.3f}% vs {100 * b['D']['share_gt1']:.3f}% | largest 99.9th percentile {a['D']['p999_max']:.3f} vs {b['D']['p999_max']:.3f} | "
        f"C-E p95 {a['E']['p95_abs_diff']:.4f} vs {b['E']['p95_abs_diff']:.4f}")
    timing(ctx, log); ctx.close()
    verdict = "PASSED" if allok else "FAILED"; log(f"CONTROLS VERDICT: {verdict}")
    if not allok: log("A control failed: section 2 will not run. Stop here; the label is 'Inconclusive, test not sensitive' until a dated amendment is made.")
    with open(VERDICT, "w") as o: o.write(f"{verdict} {time.strftime('%Y-%m-%d %H:%M:%S')} test={int(TEST)}\n")
    log(f"verdict file: {VERDICT}")

# ======================= Section 2: the residue table =======================
def section2(log):
    log("SqnPAN_07 section 2: per-residue table (resumable)")
    if not os.path.exists(VERDICT) or not open(VERDICT).read().startswith("PASSED"): raise RuntimeError("controls have not PASSED (verdict file missing or FAILED): section 2 does not run")
    log("controls verdict file:", open(VERDICT).read().strip()); check_inputs(log)
    use = X.load_use(); accs = sorted(use)
    if LIMIT: accs = accs[:LIMIT]; log(f"TEST RUN: first {len(accs)} proteins only")
    os.makedirs(PARTS, exist_ok=True)
    batches = [(i // BATCH + 1, accs[i:i + BATCH]) for i in range(0, len(accs), BATCH)]
    todo = []
    for bid, ab in batches:
        if os.path.exists(os.path.join(PARTS, f"part_{bid:05d}.tsv.gz")) and os.path.exists(os.path.join(PARTS, f"part_{bid:05d}.stats")): continue
        todo.append((bid, ab, PARTS))
    w = WORKERS or max(1, min(8, (os.cpu_count() or 2) - 2))
    log(f"{len(accs)} proteins in {len(batches)} batches of up to {BATCH}; {len(batches) - len(todo)} batches already finished (restart), {len(todo)} to do; {w} workers (cpu count {os.cpu_count()})")
    t0 = time.time(); done = 0
    if todo:
        ctx_mp = mp.get_context("spawn")
        with ctx_mp.Pool(w, initializer=X.worker_init) as pool:
            for bid, n, tot, secs in pool.imap_unordered(X.process_batch, todo):
                done += 1; el = time.time() - t0
                log(f"  batch {bid} done ({n} proteins, {tot['rows']} rows, {secs} s) | {done} of {len(todo)} this run | {round(el)} s elapsed, about {round(el / done * (len(todo) - done))} s left")
    tot = collections.Counter(); nprot = 0
    for bid, ab in batches:
        with open(os.path.join(PARTS, f"part_{bid:05d}.stats")) as f: v = f.read().split()
        nprot += int(v[1]); tot["rows"] += int(v[2]); tot["mismatch"] += int(v[3]); tot["rsa_blank"] += int(v[4]); tot["rsa_gt1"] += int(v[5])
    want_rows = sum(use[a][1] for a in accs)
    log(f"parts: {nprot} proteins, {tot['rows']} rows (UniProt lengths sum to {want_rows}) | residue letters differing from the UniProt sequence: {tot['mismatch']} | rows without a Tien value: {tot['rsa_blank']}")
    if tot["rows"] != want_rows or tot["mismatch"] != 0 or nprot != len(accs): raise RuntimeError("row count or residue letters differ from the expected values: table not merged")
    tmp = TABLE + ".tmp"
    with open(tmp, "wb") as out:
        out.write(gzip.compress(("\t".join(X.TABLE_HEADER) + "\n").encode("utf-8")))
        for bid, ab in batches:
            with open(os.path.join(PARTS, f"part_{bid:05d}.tsv.gz"), "rb") as pf: shutil.copyfileobj(pf, out)
    os.replace(tmp, TABLE)
    log(f"written {TABLE} | size {os.path.getsize(TABLE)} | md5 {H.md5_file(TABLE)}")
    log("The part files in P2_07_parts can be kept or deleted (the table is complete without them).")

# ======================= Section 3: checks on the table =======================
def section3(log):
    log("SqnPAN_07 section 3: checks on the table")
    check_inputs(log)
    if not os.path.exists(TABLE): raise RuntimeError("table not found: run section 2 first")
    ctx = X.Ctx(); use = ctx.use
    accs = sorted(use)
    if LIMIT: accs = accs[:LIMIT]
    keep = set(accs)
    need = collections.defaultdict(dict)       # accession -> position -> list of allowed letter strings
    with open(SITES, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["accession"] in keep and r["alphafold_status"] == "use":
                n = int(r["N_position"]); need[r["accession"]].setdefault(n, []).append("N"); need[r["accession"]].setdefault(n + 2, []).append("ST")
    n_sites_needed = sum(len(v) for v in need.values()) // 2
    recheck = set(random.Random(SEEDS[0] + 2).sample(accs, min(50, len(accs)))); kept_lines = collections.defaultdict(list)
    cnt = collections.Counter(); bad = collections.Counter(); last = {}; nrow = 0; sum_pl = 0.0; gt1 = 0; nr = 0; blank = 0; term = 0; pl70 = 0
    with gzip.open(TABLE, "rt", encoding="utf-8", newline="") as f:
        header = f.readline().rstrip("\n").split("\t")
        if header != X.TABLE_HEADER: raise RuntimeError(f"header differs: {header}")
        cur = None; nd = None; seq = None
        for line in f:
            line = line.rstrip("\n"); v = line.split("\t"); acc = v[0]; pos = int(v[1]); aa = v[2]
            if acc != cur:
                if acc in last: bad["accession rows not contiguous"] += 1
                cur = acc; nd = need.get(acc); seq = ctx.seqs.get(acc)
                if seq is None: bad["accession not in the usable set"] += 1; seq = ""
                last[acc] = 0
            if pos != last[acc] + 1: bad["positions not consecutive"] += 1
            last[acc] = pos; cnt[acc] += 1; nrow += 1
            if pos > len(seq) or seq[pos - 1] != aa: bad["residue letter differs from the UniProt sequence"] += 1
            if nd is not None and pos in nd:
                for allowed in nd[pos]:
                    if aa not in allowed: bad["sequon position has the wrong residue"] += 1
            sum_pl += float(v[3]); pl70 += float(v[3]) >= 70; term += int(v[8])
            if v[5] == "": blank += 1
            else:
                nr += 1; gt1 += float(v[5]) > 1.0
            if acc in recheck: kept_lines[acc].append(line)
    for a in accs:
        if cnt[a] != use[a][1]: bad["rows per protein differ from the UniProt length"] += 1
    if set(last) != keep: bad["set of proteins differs from the usable set"] += 1
    ident = 0
    for a in sorted(recheck):
        lines, _ = X.protein_rows(a, ctx); ident += int(lines == kept_lines[a])
    log(f"rows {nrow} in {len(last)} proteins | sequon positions checked: {n_sites_needed} sites (N at N position, S/T at +2)")
    log(f"recomputed proteins identical to the table: {ident} of {len(recheck)}")
    log(f"mismatch counts: {dict(bad) if bad else 'none'}")
    log(f"pooled: mean pLDDT {sum_pl / nrow:.2f} | share pLDDT >= 70 {100 * pl70 / nrow:.1f}% | share of rsa > 1 (all residues with a Tien value) {100 * gt1 / nr:.3f}% | rows without rsa {blank} | terminal rows {term}")
    ok = not bad and ident == len(recheck)
    log(f"TABLE CHECKS VERDICT: {'PASSED' if ok else 'FAILED'}")
    with open(CHECKS_VERDICT, "w") as o: o.write(f"{'PASSED' if ok else 'FAILED'} {time.strftime('%Y-%m-%d %H:%M:%S')} rows {nrow} md5 {H.md5_file(TABLE)}\n")
    log(f"table size {os.path.getsize(TABLE)} | md5 {H.md5_file(TABLE)}")
    ctx.close()

SECTIONS = {"1": ("controls_timing", section1), "2": ("residue_table", section2), "3": ("table_checks", section3)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "1"
    for k in (SECTIONS if which == "all" else [which]):
        name, fn = SECTIONS[k]; H.run_section("P2_07", k, name, fn)
# ---- END OF SCRIPT ----
