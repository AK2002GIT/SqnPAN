# SqnPAN_17_constrained_check.py - germline test at constrained Cys and His residues, follow-up of the A2 control (Amendment 4; added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (fixed 2026-11-01). Post hoc: no primary label changes.
# Question: can the germline test of SqnPAN_15 (same Engine, same variants, same expected rates) see depletion of missense variation at residues known to be under strong
#   constraint? If yes, the missing depletion at glycosite N in SqnPAN_15 is a real finding of the test. If no, the presence-based germline counts cannot show residue-level
#   selection and the A2 control of SqnPAN_15 is reported as uninformative.
# Same as SqnPAN_15: Engine of SqnPAN_14_context_helpers.py, germline variants option A (set-ups main and uniform_restricted) and option B (sensB), transcript map of SqnPAN_14, germline
#   expected rates (results\tables\P2_15_germline_spectra.tsv, rate_A and rate_B; read from P2_15_engine.npz and compared with the table), C2/C3 and the final-label rule as SqnPAN_15.
# Positions: every Cys (slot "N" of the Engine = group Cys) and every His (slot "ST" = group His) of the proteins mapped in SqnPAN_14; counted event = any missense (OPP from CHG["any"]).
# Roles (tiers), exclusive; a position with two feature types is left out; background = Cys/His of the same protein with none of the three:
#   T1 = single-residue BINDING record (UniProt 2026_03 features table) whose ligand is a metal ion, metal cation, iron-sulfur cluster or heme (regex in t1_ligand below) = DECISION UNIT
#   T2 = ACT_SITE record at one residue; T3 = Cys in a DISULFID record (either end; also a single-position interchain record). BINDING ranges are not used.
#   (The planning counts before mapping were Cys 4,706 and His 2,602 for T1; the regex fixed here gives slightly different numbers, logged in section 1.)
# Set-ups main, uniform_restricted, sensB (no strict set-up: there are no pLDDT or exposure bins for Cys/His). Groups Cys, His, both (Engine groups N, ST, both). Tiers T1 to T3.
# Expectation stated before any result: Depleted.
# DECISION UNIT: main, T1, both (Cys + His), pooled. Reading rule: Depleted -> the test sees strong germline selection, so no detectable germline depletion at glycosite N
#   in SqnPAN_15 is a real finding. Not Depleted -> presence-based germline counts cannot show residue-level selection and the A2 control is reported as uninformative.
#   No label changes either way; the A2 suffix (_biology_control_FAILED) on the tumour Equivalent stays.
# Seeds (block 33500000 to 33599999): C2/C3 runs 1 and 2 33500000, run 3 33520000, real runs 1 and 2 33540000, run 3 33560000.
# Usage (from the repository folder):  python scripts/SqnPAN_17_constrained_check.py all [workers=6]    (run after SqnPAN_14 and SqnPAN_15 have finished)
# Expected run time about 7 min (a guess from SqnPAN_15, which had four set-ups and a 105 s variation file). Logs: results\logs\P2_17_section<k>_<name>.txt.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, re, csv, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_06_helpers as P
import SqnPAN_08_selection_controls as S8
import SqnPAN_14_context_helpers as C

ARGS = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
WORKERS = int(ARGS.get("workers", "6"))
BLOCK = 33500000
SEEDS = dict(cc1=33500000, cc3=33520000, real1=33540000, real3=33560000)
CTRL_BASES = {1: SEEDS["cc1"], 2: SEEDS["cc1"], 3: SEEDS["cc3"]}; REAL_BASES = {1: SEEDS["real1"], 2: SEEDS["real1"], 3: SEEDS["real3"]}
OK14_2 = os.path.join(H.PROC_OUT, "P2_14_section2_ok.txt"); OK15_2 = os.path.join(H.PROC_OUT, "P2_15_section2_ok.txt")
POS14 = os.path.join(H.PROC_OUT, "P2_14_positions.npz"); TMAP14 = os.path.join(H.TABLES, "P2_14_transcript_map.tsv")
ENG15 = os.path.join(H.PROC_OUT, "P2_15_engine.npz"); ENG15_MD5 = "31862cb291742fe5025675709cb27f66"
VARS15 = os.path.join(H.PROC_OUT, "P2_15_variants.npz"); SPEC15 = os.path.join(H.TABLES, "P2_15_germline_spectra.tsv")
POS_NPZ = os.path.join(H.PROC_OUT, "P2_17_positions.npz"); ENGINE_NPZ = os.path.join(H.PROC_OUT, "P2_17_engine.npz")
OK1 = os.path.join(H.PROC_OUT, "P2_17_section1_ok.txt"); OK2 = os.path.join(H.PROC_OUT, "P2_17_section2_ok.txt")
PARTS = os.path.join(H.PROC_OUT, "P2_17_parts"); CAPAB = os.path.join(H.TABLES, "P2_17_unit_capability.tsv")
VERDICT = os.path.join(H.PROC_OUT, "P2_17_controls_verdict.txt"); RESULTS = os.path.join(H.TABLES, "P2_17_constrained_results.tsv")
DEC_VERDICT = os.path.join(H.PROC_OUT, "P2_17_constrained_verdict.txt")
def CTRL_TABLE(k): return os.path.join(H.TABLES, f"P2_17_C2C3_run{k}.tsv")
def REAL_TABLE(k): return os.path.join(H.TABLES, f"P2_17_constrained_run{k}.tsv")

CONFIGS = [dict(name="main", kind="main", mask="germA", mode="ctx", spec="germA", cancer=False),
           dict(name="uniform_restricted", kind="main", mask="germA", mode="uniform", spec=None, cancer=False),
           dict(name="sensB", kind="main", mask="germB", mode="ctx", spec="germB", cancer=False)]
CFG_NAMES = [c["name"] for c in CONFIGS]
RESNAME = {"N": "Cys", "ST": "His", "both": "Cys+His"}
DIS_RE = re.compile(r"DISULFID\s+(\d+)(?:\.\.(\d+))?\s*(?:;|$)")
ACT_RE = re.compile(r"ACT_SITE\s+(\d+)\s*(?:;|$)")
BIND_SPLIT = re.compile(r"(?:(?<=; )|^)(?=BINDING\s)")
BIND_RE = re.compile(r"BINDING\s+(\d+)\s*;(.*)$", re.S)
LIG_RE = re.compile(r'/ligand="([^"]*)"')

def t1_ligand(l):
    """Metal ion (Zn(2+), Cu(+), K(+) ...; not H(+)), metal cation (Fe cation, a divalent metal cation), iron-sulfur cluster or heme."""
    return bool(re.fullmatch(r"(?!H\()[A-Z][a-z]?\(\d?\+\)", l) or re.fullmatch(r"[A-Z][a-z]? cation", l) or "metal cation" in l
                or re.search(r"\bheme\b|siroheme", l) or "cluster" in l.lower())

def need_ok(path, what):
    if not os.path.exists(path) or not open(path).read().startswith("OK"): raise RuntimeError(f"{what} has not finished ({os.path.basename(path)} missing)")

# ======================= Section 1: roles, positions, codon opportunities =======================
def section1(log):
    need_ok(OK14_2, "P2_14 section 2"); need_ok(OK15_2, "P2_15 section 2")
    log("Pre-specification: header of this script (fixed 2026-11-01). Post hoc, added after the primary analysis.")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__}")
    U = H.load_uniprot_main(); accs = sorted(U); F = H.load_uniprot_features(); t0 = time.time()
    if accs != [str(x) for x in C.load_npz(POS14)["accs"]]: raise RuntimeError("accession order differs from SqnPAN_14")
    G, _, _ = C.load_gencode(); M = {}
    with open(TMAP14, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            vl = G[r["transcript"]]; vi = [x["ver"] for x in vl].index(r["version"]); M[r["accession"]] = (r["transcript"], vi)
    bad = sum(1 for a, (tid, vi) in M.items() if C.translate(G[tid][vi]["seq"][G[tid][vi]["cs"]:G[tid][vi]["ce"]]) != U[a]["seq"])
    log(f"SqnPAN_14 transcript map: {len(M)} proteins | translation differs from UniProt: {bad} (must be 0) | GENCODE read {time.time() - t0:.0f} s")
    if bad: raise RuntimeError("transcript map does not reproduce")
    T = collections.Counter(); lig = collections.Counter()
    P_prot, P_res, P_posn, P_role = [], [], [], []
    for a_i, acc in enumerate(accs):
        if acc not in M: continue
        s = U[acc]["seq"]; fr = F.get(acc); typ = collections.defaultdict(set)
        if fr is not None:
            for m in DIS_RE.finditer(fr.get("Disulfide bond") or ""):
                for g in m.groups():
                    if g is None: continue
                    p = int(g)
                    if 1 <= p <= len(s) and s[p - 1] == "C": typ[p].add(3)
                    else: T["disulfide positions not at a Cys (left out)"] += 1
            for m in ACT_RE.finditer(fr.get("Active site") or ""):
                p = int(m.group(1))
                if 1 <= p <= len(s) and s[p - 1] in "CH": typ[p].add(2)
                else: T["active-site positions not at Cys or His (left out)"] += 1
            for part in BIND_SPLIT.split(fr.get("Binding site") or ""):
                if not part.startswith("BINDING"): continue
                m = BIND_RE.match(part)
                if not m: T["BINDING records that are ranges or uncertain (not used)"] += 1; continue
                p = int(m.group(1))
                if not (1 <= p <= len(s)) or s[p - 1] not in "CH": continue
                lm = LIG_RE.search(m.group(2)); l = lm.group(1) if lm else ""
                if t1_ligand(l): typ[p].add(1); lig[l] += 1
                else: T["single-residue BINDING at Cys/His with another ligand (stays background)"] += 1
        for mm in re.finditer("[CH]", s):
            j = mm.start() + 1; t = typ.get(j)
            role = 0 if not t else (next(iter(t)) if len(t) == 1 else -1)
            P_prot.append(a_i); P_res.append(0 if mm.group() == "C" else 1); P_posn.append(j); P_role.append(role)
    P_prot = np.array(P_prot, np.int32); P_res = np.array(P_res, np.int8); P_posn = np.array(P_posn, np.int32); P_role = np.array(P_role, np.int8)
    log(f"Cys and His positions in the {len(M)} mapped proteins: {len(P_prot)}")
    log("--- roles before the codon restriction (Cys / His) ---")
    for t, nm in ((0, "background"), (1, "T1 metal/cluster/heme binding"), (2, "T2 active site"), (3, "T3 disulfide"), (-1, "left out: two feature types")):
        log(f"  {nm}: Cys {int(((P_role == t) & (P_res == 0)).sum())} | His {int(((P_role == t) & (P_res == 1)).sum())}")
    for k in sorted(T): log(f"  {k}: {T[k]}")
    log("  most frequent T1 ligands: " + "; ".join(f"{k} {v}" for k, v in lig.most_common(12)))
    n = len(P_prot); OPP = np.full((n, 9), -1, np.int16); ctx = np.zeros(n, bool); c = collections.Counter(); last = -1; Tx = None; cs = 0; acc = None
    pl = P_prot.tolist(); rl = P_res.tolist(); ql = P_posn.tolist()
    for k in range(n):
        if pl[k] != last:
            last = pl[k]; acc = accs[last]; d = G[M[acc][0]][M[acc][1]]; Tx = d["seq"]; cs = d["cs"]; s = U[acc]["seq"]
        i = cs + 3 * (ql[k] - 1); cod = Tx[i:i + 3]
        if C.CODON.get(cod) != s[ql[k] - 1]: c["codon does not code the UniProt residue (must be 0)"] += 1; continue
        j = 0
        for kk, alt, _a in C.CHG["any"][cod]:
            idx = i + kk
            if idx < 1 or idx + 1 >= len(Tx): c["changes without both flanking bases (left out)"] += 1; continue
            ch = C.chan_tx(Tx[idx - 1:idx + 2], alt)
            if ch < 0: c["changes with a non-ACGT flank (left out)"] += 1; continue
            OPP[k, j] = ch; j += 1
        if j: ctx[k] = True; c["positions with codon opportunities"] += 1
        else: c["positions without a usable change"] += 1
    for k in sorted(c): log(f"  {k}: {c[k]}")
    if c["codon does not code the UniProt residue (must be 0)"]: raise RuntimeError("codon check failed")
    log("--- roles after the codon restriction (positions with at least one usable change; Cys / His) ---")
    for t, nm in ((0, "background"), (1, "T1"), (2, "T2"), (3, "T3")):
        log(f"  {nm}: Cys {int(((P_role == t) & (P_res == 0) & ctx).sum())} | His {int(((P_role == t) & (P_res == 1) & ctx).sum())} | proteins {len(set(P_prot[(P_role == t) & ctx].tolist()))}")
    np.savez_compressed(POS_NPZ, accs=np.array(accs), P_prot=P_prot, P_res=P_res, P_posn=P_posn, P_role=P_role, P_ctx=ctx, OPP=OPP)
    with open(OK1, "w") as o: o.write(f"OK\n{os.path.basename(POS_NPZ)}\t{H.md5_file(POS_NPZ)}\n")
    log(f"saved {POS_NPZ} | md5 {H.md5_file(POS_NPZ)} | {time.time() - t0:.0f} s")

# ======================= Section 2: events and engine inputs =======================
def section2(log):
    need_ok(OK1, "section 1")
    got = H.md5_file(ENG15); log(f"md5 {os.path.basename(ENG15)}: {got} ({'same as recorded' if got == ENG15_MD5 else 'DIFFERENT, recorded ' + ENG15_MD5})")
    if got != ENG15_MD5: raise RuntimeError("SqnPAN_15 engine file md5 differs: stop and decide")
    U = H.load_uniprot_main(); accs = sorted(U); Pz = C.load_npz(POS_NPZ); E15 = C.load_npz(ENG15); Vz = C.load_npz(VARS15)
    if accs != [str(x) for x in Pz["accs"]]: raise RuntimeError("accession order differs")
    rA = E15["R_germA"]; rB = E15["R_germB"]
    with open(SPEC15, newline="", encoding="utf-8") as f: sp = list(csv.DictReader(f, delimiter="\t"))
    if [r["channel"] for r in sp] != C.CHANNELS: raise RuntimeError("channel order of the SqnPAN_15 spectra table differs")
    dA = max(abs(float(r["rate_A"]) - rA[i, 0]) / max(rA[i, 0], 1e-12) for i, r in enumerate(sp)); dB = max(abs(float(r["rate_B"]) - rB[i, 0]) / max(rB[i, 0], 1e-12) for i, r in enumerate(sp))
    log(f"expected rates of SqnPAN_15 (engine file against spectra table): largest relative difference A {dA:.2e}, B {dB:.2e} (must be < 1e-4)")
    if dA > 1e-4 or dB > 1e-4: raise RuntimeError("expected rates differ")
    uk = Vz["uk"]; ai, ps = uk >> 25, (uk >> 9) & 65535; optA = Vz["optA"]; optB = Vz["optB"]
    offs = np.zeros(len(accs) + 1, np.int64)
    for i, a in enumerate(accs): offs[i + 1] = offs[i] + len(U[a]["seq"])
    posmap = np.full(int(offs[-1]), -1, np.int64)
    posmap[offs[Pz["P_prot"].astype(np.int64)] + Pz["P_posn"].astype(np.int64) - 1] = np.arange(len(Pz["P_posn"]))
    pidx = posmap[offs[ai] + ps - 1]; counted = (pidx >= 0) & optA; ev = np.flatnonzero(counted)
    n = len(Pz["P_prot"])
    E = dict(accs=Pz["accs"], P_prot=Pz["P_prot"], P_res=Pz["P_res"], P_role=Pz["P_role"], P_usable=np.zeros(n, bool), P_pb=np.full(n, -1, np.int8), P_rb=np.full(n, -1, np.int8),
             P_ctx=Pz["P_ctx"], OPP=Pz["OPP"], projects=np.array([], dtype="<U1"), ev_pos=pidx[ev].astype(np.int32), ev_cancer=np.zeros(len(ev), np.int8),
             ev_germA=optA[ev], ev_germB=optB[ev], R_germA=rA, R_germB=rB)
    np.savez_compressed(ENGINE_NPZ, **E)
    role = E["P_role"]; res = E["P_res"]; ctx = E["P_ctx"]; evp = E["ev_pos"]
    log(f"germline variants (option A) at Cys/His positions of the mapped proteins: {len(ev)} (option B {int(optB[ev].sum())}); variants of SqnPAN_15 {len(uk)}")
    log("--- germline events (option A; all / at positions with codon opportunities) ---")
    for t, nm in ((0, "background"), (1, "T1"), (2, "T2"), (3, "T3")):
        for rc, rn in ((0, "Cys"), (1, "His")):
            mk = (role == t) & (res == rc); log(f"  {nm} {rn}: positions {int((mk & ctx).sum())} | events {int(mk[evp].sum())} / {int((mk & ctx)[evp].sum())}")
    log(f"engine inputs {ENGINE_NPZ} | md5 {H.md5_file(ENGINE_NPZ)}")
    with open(OK2, "w") as o: o.write(f"OK\n{os.path.basename(ENGINE_NPZ)}\t{H.md5_file(ENGINE_NPZ)}\n")

# ======================= Section 3: controls C2, C3 =======================
def section3(log):
    need_ok(OK2, "section 2"); C.check_seed_blocks(log, BLOCK, SEEDS.values())
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__} | workers {WORKERS} | repeats {S8.N_REP} | bootstrap {S8.B}")
    C.run_controls(log, ENGINE_NPZ, CONFIGS, True, CTRL_BASES, PARTS, CTRL_TABLE, WORKERS)
    cap, c2ok, same12 = C.capability(log, CTRL_TABLE, CFG_NAMES, CAPAB)
    pk = ("main", "T1", "both", "pooled")
    lines = ["PASSED" if same12 else "FAILED", f"C2/C3 run1 and run2 identical: {'yes' if same12 else 'no'}",
             *[f"C2 run{k} {cfg}: {'PASSED' if c2ok[(k, cfg)] else 'FAILED'}" for k in (1, 3) for cfg in CFG_NAMES],
             f"decision unit (main, T1, both, pooled): {cap[pk][0]}", f"written {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    with open(VERDICT, "w") as o: o.write("\n".join(lines) + "\n")
    for x in lines: log("VERDICT " + x)

# ======================= Section 4: real result and decision =======================
def section4(log):
    if not os.path.exists(VERDICT) or open(VERDICT).readline().strip() != "PASSED": raise RuntimeError("controls verdict is not PASSED: no real result is computed")
    C.check_seed_blocks(log, BLOCK, SEEDS.values())
    eng = C.Engine(C.load_npz(ENGINE_NPZ), CONFIGS, restricted=True)
    C.real_runs(log, eng, REAL_BASES, REAL_TABLE); cap = C.read_cap(CAPAB)
    runs, fin = C.assemble(log, REAL_TABLE, cap, CFG_NAMES, RESULTS, extra_hdr=("residues",), extra_fn=lambda key: [RESNAME[key[2]]])
    pk = ("main", "T1", "both", "pooled")
    log("=== DECISION UNIT (germline, option A, T1 metal/cluster/heme-binding Cys + His, germline expected rates); expectation stated before the result: Depleted ===")
    log("  " + C.unit_line(runs, fin, cap, pk))
    log("=== POOLED UNITS (residues: Cys, His, Cys+His) ===")
    for cfg in CFG_NAMES:
        for t in ("T1", "T2", "T3"):
            for g in ("N", "ST", "both"):
                if (cfg, t, g, "pooled") in runs[1]: log(f"  {cfg} {t} {RESNAME[g]}: " + C.unit_line(runs, fin, cap, (cfg, t, g, "pooled")))
    passed = fin[1][pk][0] == "Depleted"
    reading = ("READING: the germline test sees strong selection at residues known to be constrained, so the missing germline depletion at glycosite N in P2_15 is a real finding of the test."
               if passed else
               "READING: the presence-based germline counts do not show depletion even at residues known to be constrained, so the A2 control of P2_15 is reported as uninformative about glycosites.")
    lines = ["PASSED" if passed else "FAILED", f"decision unit final label: {fin[1][pk][0]}" + (f" [{fin[1][pk][1]}]" if fin[1][pk][1] else ""), reading,
             "No label changes; the suffix _biology_control_FAILED of P2_15 stays on the tumour Equivalent.", f"written {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    with open(DEC_VERDICT, "w") as o: o.write("\n".join(lines) + "\n")
    for x in lines: log("CONSTRAINED-RESIDUE CHECK " + x)

SECTIONS = {"1": ("roles_positions", section1), "2": ("events_engine", section2), "3": ("controls_C2_C3", section3), "4": ("real_result_decision", section4)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    C.run_chain("P2_17", SECTIONS, "1" if which == "all" else which)
# ---- END OF SCRIPT ----
