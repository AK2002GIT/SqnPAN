# SqnPAN_15_germline_q2.py - A2 germline positive control (Amendment 4; added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-10-31; decision unit T1 N + S/T; option-B sensitivity included).
#   Post hoc: no primary label changes by itself; the reading rule below only adds a suffix to the tumour Equivalent.
# Question: does the same Q2 code find the expected depletion of germline missense variation at N-glycosite residues? If it does not, the tumour Equivalent cannot be
#   read as "selection absent", only as "the test did not separate".
# Expectation stated before any result: Depleted. Mazumder et al. 2012 (PLoS ONE 7(5):e36212, doi 10.1371/journal.pone.0036212; dbSNP via UniProt 2011_01 and PolyPhen)
#   report "slightly lower" variation at the motif N than at other N, with no test; from their Table 3 (computed here from the published counts, not printed in the paper): N 126/18,662
#   against 4,188/419,793 (ratio about 0.68), S 54/10,268 against 8,390/981,010 (about 0.62), T 90/8,394 against 7,265/609,107 (about 0.90). Not blind.
# Variants (option A): UniProt homo_sapiens_variation 2026_03 (md5 ddfc80d965ebafb900d64fe8c51d7cec), Consequence Type "missense variant",
#   AC equal to a reviewed 2026_03 accession (isoform accessions "-n" left out), reference residue equal to UniProt 2026_03 at the position (mismatches dropped, counted).
#   Germline evidence = gnomAD, TOPMed, ExAC, 1000Genomes or ESP in the Evidence column. Somatic flag = Evidence NCI-TCGA, NCI-TCGACosmic or cosmiccurated, or a Source DB
#   ID starting COSV, COSM, TCGA or NCI. Option A: germline evidence and no somatic flag in any row of the same protein variant (accession, position, alt).
#   Option B (pre-specified sensitivity): germline evidence from two or more of the five resources, same somatic rule.
#   Event = one distinct nucleotide change at a position (accession, position, alt, genomic change from "Chromosome Coordinate"; rows without an SNV coordinate keep one
#   event per protein variant). Counted at N: any missense; at S/T: missense to a residue that is not S or T (as in Q2). No allele frequency exists in the file, so no
#   frequency rule is possible (stated limit); ClinVar-only rows have no germline evidence and are out.
# Germline expected rates: positions of the proteins mapped in SqnPAN_14 (same transcripts), outside every N-X-S/T motif (N, X and S/T positions left out). Channel of a
#   variant: the single-nucleotide change of the codon that gives the alt residue and whose bases equal the genomic change (either strand); if several remain, the count
#   is split equally; if none, no channel (counted). Rate of channel k = variants in k / possible missense changes in k at the same positions. Position exposure = sum of
#   the rates of its counted changes (SqnPAN_14 opportunities). GATE (spectrum sanity): mean rate of the four CpG C>T channels >= 2 x mean rate of the other 92.
# Test: the Engine of SqnPAN_14_context_helpers.py (SqnPAN_08 / SqnPAN_09 code, proved line for line by control C-S of SqnPAN_14). Set-ups: main (option A, germline expected rates),
#   strict (option A, also pLDDT and exposure bins), uniform_restricted (option A, position counts on the same proteins), sensB (option B, its own spectrum). Tiers T1 to
#   T3, groups N, S/T, both; no cancer types. 36 units. Labels, BH families (decision unit alone; the other 8 pooled units of a set-up), tested-unit rule as SqnPAN_08.
# Controls: C2 and C3 as SqnPAN_08 (20 repeats; runs 1 and 2 identical; run 3 other seed); final-label rule 1A of SqnPAN_09.
# DECISION (germline positive control A2): unit main, T1, N + S/T. PASSED if its final label is Depleted. If not, the tumour primary Equivalent (SqnPAN_09) and the A1 primary
#   label (SqnPAN_14) are read with the suffix _biology_control_FAILED. T1 N, T1 S/T and the other set-ups are reported beside it.
# Seeds (sub-bases inside the block reserved for A2, 32800000 to 32899999): C2/C3 run1 and run2 32800000, run3 32820000; real runs 1 and 2 32840000, run 3 32860000.
# Usage (from the repository folder):  python scripts/SqnPAN_15_germline_q2.py all [workers=6]     (run after SqnPAN_14 has finished)
# Logs: results\logs\P2_15_section<k>_<name>.txt.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, re, csv, gzip, time, collections
from array import array
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_06_helpers as P
import SqnPAN_08_selection_controls as S8
import SqnPAN_14_context_helpers as C

ARGS = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
WORKERS = int(ARGS.get("workers", "6"))
BLOCK = 32800000
SEEDS = dict(cc1=32800000, cc3=32820000, real1=32840000, real3=32860000)
CTRL_BASES = {1: SEEDS["cc1"], 2: SEEDS["cc1"], 3: SEEDS["cc3"]}; REAL_BASES = {1: SEEDS["real1"], 2: SEEDS["real1"], 3: SEEDS["real3"]}
VAR = os.path.join(H.RAW, "Uniprot_variation", "homo_sapiens_variation.txt.gz"); VAR_MD5 = "ddfc80d965ebafb900d64fe8c51d7cec"
OK14_1 = os.path.join(H.PROC_OUT, "P2_14_section1_ok.txt"); OK14_2 = os.path.join(H.PROC_OUT, "P2_14_section2_ok.txt")
POS14 = os.path.join(H.PROC_OUT, "P2_14_positions.npz"); ENGINE14 = os.path.join(H.PROC_OUT, "P2_14_engine.npz"); TMAP14 = os.path.join(H.TABLES, "P2_14_transcript_map.tsv")
RES14 = os.path.join(H.TABLES, "P2_14_Q2ctx_results.tsv")
OK1 = os.path.join(H.PROC_OUT, "P2_15_section1_ok.txt"); OK2 = os.path.join(H.PROC_OUT, "P2_15_section2_ok.txt")
VARS_NPZ = os.path.join(H.PROC_OUT, "P2_15_variants.npz"); ENGINE_NPZ = os.path.join(H.PROC_OUT, "P2_15_engine.npz")
SPECTRA = os.path.join(H.TABLES, "P2_15_germline_spectra.tsv"); PARTS = os.path.join(H.PROC_OUT, "P2_15_parts"); CAPAB = os.path.join(H.TABLES, "P2_15_unit_capability.tsv")
VERDICT = os.path.join(H.PROC_OUT, "P2_15_controls_verdict.txt"); RESULTS = os.path.join(H.TABLES, "P2_15_germline_results.tsv")
BIO_VERDICT = os.path.join(H.PROC_OUT, "P2_15_biology_control_verdict.txt")
def CTRL_TABLE(k): return os.path.join(H.TABLES, f"P2_15_C2C3_run{k}.tsv")
def REAL_TABLE(k): return os.path.join(H.TABLES, f"P2_15_germline_run{k}.tsv")

CONFIGS_A2 = [dict(name="main", kind="main", mask="germA", mode="ctx", spec="germA", cancer=False),
              dict(name="strict", kind="strict", mask="germA", mode="ctx", spec="germA", cancer=False),
              dict(name="uniform_restricted", kind="main", mask="germA", mode="uniform", spec=None, cancer=False),
              dict(name="sensB", kind="main", mask="germB", mode="ctx", spec="germB", cancer=False)]
CFG_NAMES = [c["name"] for c in CONFIGS_A2]
THREE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M",
         "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
AAI = {a: i for i, a in enumerate(P.AA)}
GPAIRS = [(r, a) for r in "ACGT" for a in "ACGT" if r != a]; GCH = {p: i + 1 for i, p in enumerate(GPAIRS)}
TOK = {"gnomAD": 1, "TOPMed": 2, "ExAC": 4, "1000Genomes": 8, "ESP": 16}; SOM_TOK = {"NCI-TCGA", "NCI-TCGACosmic", "cosmiccurated"}; SOM_ID = ("COSV", "COSM", "TCGA", "NCI")
PV = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$"); GV = re.compile(r":g\.\d+([ACGT])>([ACGT])$")
POPC = np.array([bin(i).count("1") for i in range(32)])

def need_ok(path, what):
    if not os.path.exists(path) or not open(path).read().startswith("OK"): raise RuntimeError(f"{what} has not finished ({os.path.basename(path)} missing)")

def decode(uk):
    return uk >> 25, (uk >> 9) & 65535, (uk >> 4) & 31, uk & 15     # accession index, position, alt index, genomic change index

# ======================= Section 1: variants =======================
def section1(log):
    need_ok(OK14_2, "P2_14 section 2")
    log("Pre-specification: header of this script (fixed 2026-10-31). Post hoc, added after the primary analysis.")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__}")
    got = H.md5_file(VAR); log(f"md5 {os.path.basename(VAR)}: {got} ({'same as recorded' if got == VAR_MD5 else 'DIFFERENT, recorded ' + VAR_MD5})")
    if got != VAR_MD5: raise RuntimeError("variation file md5 differs: stop and decide")
    U = H.load_uniprot_main(); accs = sorted(U); aidx = {a: i for i, a in enumerate(accs)}
    if accs != [str(x) for x in C.load_npz(POS14)["accs"]]: raise RuntimeError("accession order differs from SqnPAN_14")
    keys = array("q"); flags = array("b"); T = collections.Counter(); t0 = time.time(); n = 0
    with gzip.open(VAR, "rt", encoding="utf-8", newline="") as f:
        for line in f:
            if line.startswith("Gene Name"): hdr = [h.strip() for h in line.rstrip("\n").split("\t")]; break
        else: raise RuntimeError("data header line not found")
        need = ["AC", "Variant AA Change", "Source DB ID", "Consequence Type", "Chromosome Coordinate", "Evidence"]
        miss = [k for k in need if k not in hdr]
        if miss: raise RuntimeError(f"columns missing: {miss}")
        iac, iaa, isid, icon, icoord, iev = (hdr.index(k) for k in need); next(f)
        log(f"columns found: {need} at {[hdr.index(k) for k in need]}")
        for line in f:
            n += 1
            if "missense variant" not in line: continue
            v = line.rstrip("\n").split("\t")
            if v[icon].strip() != "missense variant": continue
            T["missense variant rows"] += 1; ac = v[iac].strip(); i = aidx.get(ac)
            if i is None: T["rows on isoform accessions (left out)" if "-" in ac else "rows on accessions not in the reviewed set (left out)"] += 1; continue
            m = PV.match(v[iaa].strip())
            if not m: T["rows with an unparsable HGVS change (left out)"] += 1; continue
            ref = THREE.get(m.group(1)); alt = THREE.get(m.group(3)); pos = int(m.group(2))
            if ref is None or alt is None or ref == alt: T["rows with a non-standard or unchanged residue (left out)"] += 1; continue
            sq = U[ac]["seq"]
            if pos > len(sq) or pos > 65535 or sq[pos - 1] != ref: T["rows whose reference residue differs from UniProt 2026_03 (left out)"] += 1; continue
            gm = GV.search(v[icoord].strip()); g = GCH[(gm.group(1), gm.group(2))] if gm else 0
            fl = 0
            for tok in v[iev].split(","):
                tok = tok.strip(); fl |= TOK.get(tok, 0)
                if tok in SOM_TOK: fl |= 32
            if v[isid].strip().startswith(SOM_ID): fl |= 32
            keys.append((((i << 16) | pos) << 5 | AAI[alt]) << 4 | g); flags.append(fl); T["rows kept"] += 1
            if n % 5000000 == 0: log(f"  rows read {n} | {time.time() - t0:.0f} s")
    log(f"variation file: {n} data rows | {time.time() - t0:.0f} s"); [log(f"  {k}: {T[k]}") for k in sorted(T)]
    K = np.frombuffer(keys, dtype=np.int64); F = np.frombuffer(flags, dtype=np.int8).astype(np.int16)
    o = np.argsort(K, kind="stable"); Ks = K[o]; Fs = F[o]; uk, st = np.unique(Ks, return_index=True); uf = np.bitwise_or.reduceat(Fs, st)
    pv, inv = np.unique(uk >> 4, return_inverse=True); som_pv = np.zeros(len(pv), np.int16); np.maximum.at(som_pv, inv.ravel(), (uf & 32))
    som = som_pv[inv.ravel()] > 0; germ = uf & 31
    optA = (germ > 0) & ~som; optB = (POPC[germ] >= 2) & ~som
    log(f"distinct variants (accession, position, alt, genomic change) {len(uk)} | protein variants {len(pv)} | with somatic flag in any row of the protein variant {int(som.sum())}")
    log(f"  option A (germline evidence, no somatic flag) {int(optA.sum())} | option B (two or more germline resources) {int(optB.sum())} | without an SNV coordinate {int(((uk & 15) == 0).sum())}")
    for nm, bit in TOK.items(): log(f"  germline resource {nm}: {int(((uf & bit) > 0).sum())} variants")
    np.savez_compressed(VARS_NPZ, uk=uk, uf=uf, optA=optA, optB=optB, som=som)
    with open(OK1, "w") as o_: o_.write(f"OK\n{os.path.basename(VARS_NPZ)}\t{H.md5_file(VARS_NPZ)}\n")
    log(f"saved {VARS_NPZ}")

# ======================= Section 2: germline spectra and engine inputs =======================
def section2(log):
    need_ok(OK1, "section 1")
    U = H.load_uniprot_main(); accs = sorted(U); G, _, _ = C.load_gencode(); t0 = time.time()
    M = {}
    with open(TMAP14, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            vl = G[r["transcript"]]; vi = [x["ver"] for x in vl].index(r["version"]); M[r["accession"]] = (r["transcript"], vi)
    bad = sum(1 for a, (tid, vi) in M.items() if C.translate(G[tid][vi]["seq"][G[tid][vi]["cs"]:G[tid][vi]["ce"]]) != U[a]["seq"])
    log(f"SqnPAN_14 transcript map: {len(M)} proteins | translation differs from UniProt: {bad} (must be 0)")
    if bad: raise RuntimeError("transcript map does not reproduce")
    Vz = C.load_npz(VARS_NPZ); uk = Vz["uk"]; ai, ps, al, gi = decode(uk); optA = Vz["optA"]; optB = Vz["optB"]
    first = np.searchsorted(ai, np.arange(len(accs) + 1))
    opp = collections.Counter(); evA = np.zeros(96); evB = np.zeros(96); T = collections.Counter(); cache = {}
    GINV = {i + 1: p for i, p in enumerate(GPAIRS)}; AAL = P.AA
    for a_i, acc in enumerate(accs):
        mp_ = M.get(acc)
        if mp_ is None: continue
        d = G[mp_[0]][mp_[1]]; Tx = d["seq"]; cs = d["cs"]; sq = U[acc]["seq"]
        excl = set()
        for p in H.sequon_positions(sq)[0]: excl.update((p, p + 1, p + 2))
        five = {}
        for p in range(1, len(sq) + 1):
            if p in excl: continue
            i = cs + 3 * (p - 1)
            if i < 1 or i + 4 > len(Tx): T["positions without both flanks (left out)"] += 1; continue
            five[p] = Tx[i - 1:i + 4]
        opp.update(five.values())
        for j in range(first[a_i], first[a_i + 1]):
            p = int(ps[j])
            if p not in five: continue
            ck = (five[p], AAL[al[j]], int(gi[j])); ch = cache.get(ck)
            if ch is None:
                fv, alt_aa, g = ck; cod = fv[1:4]; cand = []
                for k, alt, aa2 in C.CHG["any"].get(cod, ()):
                    if aa2 != alt_aa: continue
                    if g:
                        r_, a_ = GINV[g]
                        if not ((cod[k] == r_ and alt == a_) or (cod[k] == C.comp(r_) and alt == C.comp(a_))): continue
                    c = C.chan_tx(fv[k:k + 3], alt)
                    if c >= 0: cand.append(c)
                cache[ck] = ch = cand
            if not ch: T["variants without a channel (codon or genomic change does not fit)"] += 1; continue
            w = 1.0 / len(ch); T["variants with a channel" + (" (split)" if len(ch) > 1 else "")] += 1
            for c in ch:
                if optA[j]: evA[c] += w
                if optB[j]: evB[c] += w
    oppc = np.zeros(96)
    for fv, nn in opp.items():
        for k, alt, _a in C.CHG["any"].get(fv[1:4], ()):
            c = C.chan_tx(fv[k:k + 3], alt)
            if c >= 0: oppc[c] += nn
    for k in sorted(T): log(f"  {k}: {T[k]}")
    log(f"spectrum inputs: {time.time() - t0:.0f} s | possible missense changes {int(oppc.sum())} | option A variants with a channel {evA.sum():.0f} | option B {evB.sum():.0f}")
    if (oppc <= 0).any(): raise RuntimeError("a channel has no opportunity")
    rA = evA / oppc; rB = evB / oppc; others = [k for k in range(96) if k not in C.CPG_CT]; ok = True
    with open(SPECTRA, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["channel", "opportunities", "variants_A", "rate_A", "variants_B", "rate_B"])
        for k in range(96): w.writerow([C.CHANNELS[k], int(oppc[k]), f"{evA[k]:.3f}", f"{rA[k]:.6g}", f"{evB[k]:.3f}", f"{rB[k]:.6g}"])
    for nm, r in (("A", rA), ("B", rB)):
        ratio = r[C.CPG_CT].mean() / r[others].mean(); log(f"sanity option {nm}: mean rate CpG C>T {r[C.CPG_CT].mean():.4g} / other channels {r[others].mean():.4g} = {ratio:.2f} (limit >= 2)")
        ok = ok and ratio >= 2.0
    if not ok: raise RuntimeError("germline spectrum sanity failed: stop and decide")
    E14 = C.load_npz(ENGINE14); P_posn = C.load_npz(POS14)["P_posn"]
    offs = np.zeros(len(accs) + 1, np.int64)
    for i, a in enumerate(accs): offs[i + 1] = offs[i] + len(U[a]["seq"])
    posmap = np.full(int(offs[-1]), -1, np.int64); posmap[offs[E14["P_prot"].astype(np.int64)] + P_posn.astype(np.int64) - 1] = np.arange(len(P_posn))
    pidx = posmap[offs[ai] + ps - 1]; resv = np.where(pidx >= 0, E14["P_res"][np.maximum(pidx, 0)], -1)
    st_alt = np.isin(al, [AAI["S"], AAI["T"]]); counted = (pidx >= 0) & ((resv == 0) | ((resv == 1) & ~st_alt)) & optA
    ev = np.flatnonzero(counted)
    E = {k: E14[k] for k in ("accs", "P_prot", "P_res", "P_role", "P_usable", "P_pb", "P_rb", "P_ctx", "OPP")}
    E.update(projects=np.array([], dtype="<U1"), ev_pos=pidx[ev].astype(np.int32), ev_cancer=np.zeros(len(ev), np.int8), ev_germA=optA[ev], ev_germB=optB[ev],
             R_germA=rA[:, None], R_germB=rB[:, None])
    np.savez_compressed(ENGINE_NPZ, **E)
    role = E["P_role"]; res = E["P_res"]; ctx = E["P_ctx"]; evp = E["ev_pos"]
    log("--- germline events at Q2 positions (option A; all / at positions with codon opportunities) ---")
    for t in (0, 1, 2, 3):
        for rc, rn in ((0, "N"), (1, "S/T")):
            mk = (role == t) & (res == rc); log(f"  {'background' if t == 0 else 'glycosite T' + str(t)} {rn}: events {int(mk[evp].sum())} / {int((mk & ctx)[evp].sum())}")
    log(f"engine inputs {ENGINE_NPZ} | md5 {H.md5_file(ENGINE_NPZ)} | events {len(ev)} (option B {int(optB[ev].sum())})")
    with open(OK2, "w") as o: o.write(f"OK\n{os.path.basename(ENGINE_NPZ)}\t{H.md5_file(ENGINE_NPZ)}\n")

# ======================= Section 3: controls C2, C3 =======================
def section3(log):
    need_ok(OK2, "section 2"); C.check_seed_blocks(log, BLOCK, SEEDS.values())
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__} | workers {WORKERS} | repeats {S8.N_REP} | bootstrap {S8.B}")
    C.run_controls(log, ENGINE_NPZ, CONFIGS_A2, True, CTRL_BASES, PARTS, CTRL_TABLE, WORKERS)
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
    eng = C.Engine(C.load_npz(ENGINE_NPZ), CONFIGS_A2, restricted=True)
    C.real_runs(log, eng, REAL_BASES, REAL_TABLE); cap = C.read_cap(CAPAB)
    runs, fin = C.assemble(log, REAL_TABLE, cap, CFG_NAMES, RESULTS)
    pk = ("main", "T1", "both", "pooled")
    log("=== DECISION UNIT (germline, option A, T1, N + S/T, germline expected rates); expectation stated before the result: Depleted ===")
    log("  " + C.unit_line(runs, fin, cap, pk))
    log("=== POOLED UNITS ===")
    for cfg in CFG_NAMES:
        for t in ("T1", "T2", "T3"):
            for g in ("N", "ST", "both"): log(f"  {cfg} {t} {g}: " + C.unit_line(runs, fin, cap, (cfg, t, g, "pooled")))
    passed = fin[1][pk][0] == "Depleted"
    p209 = C.read_run(C.RES_P209)[pk]["final_label"]
    p214 = C.read_run(RES14)[pk]["final_label"] if os.path.exists(RES14) else "(P2_14 results missing)"
    sfx = "" if passed else "_biology_control_FAILED"
    lines = ["PASSED" if passed else "FAILED", f"decision unit final label: {fin[1][pk][0]}" + (f" [{fin[1][pk][1]}]" if fin[1][pk][1] else ""),
             f"reading of the tumour primary (P2_09): {p209}{sfx}", f"reading of the A1 primary (P2_14): {p214}{sfx}", f"written {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    with open(BIO_VERDICT, "w") as o: o.write("\n".join(lines) + "\n")
    for x in lines: log("BIOLOGY CONTROL " + x)

SECTIONS = {"1": ("variants", section1), "2": ("germline_spectra", section2), "3": ("controls_C2_C3", section3), "4": ("real_result_decision", section4)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    C.run_chain("P2_15", SECTIONS, "1" if which == "all" else which)
# ---- END OF SCRIPT ----
