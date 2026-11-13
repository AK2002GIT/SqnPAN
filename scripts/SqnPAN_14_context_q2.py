# SqnPAN_14_context_q2.py - A1 context-aware Q2 (Amendment 4; added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-10-31). Post hoc: no label of SqnPAN_09 (Q2 primary unit Equivalent)
#   changes; if A1 disagrees, both are reported and the primary stays primary.
# Question: does the Q2 result (glycosite N and S/T against the same residue type elsewhere in the same protein) hold when every position's exposure is its expected
#   somatic mutation rate from its trinucleotide context, and not a plain count of positions?
# Section 1, structure check (stops the run if a gate fails):
#   - md5 of every input (GENCODE v36 pc_transcripts 61eaea20..., the SqnPAN_06 tables, UniProt main, P2_08_inputs.npz).
#   - Transcript per protein: a GENCODE v36 transcript whose CDS translation equals the UniProt 2026_03 sequence exactly (internal stop or non-ACGT: not used);
#     preference: listed in the UniProt Ensembl field without isoform tag or with tag -1, then lowest ENST number, then version order (not PAR_Y first).
#     Proteins without an exact match are left out of A1 and counted (positions and events reported).
#   - The SqnPAN_08 position list is rebuilt (with position numbers) and must equal P2_08_inputs.npz (P_prot, P_res, P_role).
#   - All 10,640 MAF files: MAF Transcript_ID matched to GENCODE without version; for every missense row the reference codon (Codons column) must equal the
#     GENCODE CDS codon at Protein_position. IDs with two versions in GENCODE (PAR_Y copies): a row matches if any version matches; the version used is counted.
#     Every mismatch is written to P2_14_codon_mismatches.tsv with its reason (this traces the 37 mismatches seen in a pilot sample). GATE: matched share >= 0.995 of
#     parsed missense rows on GENCODE transcripts. CONTEXT: GATE: share of SNV rows whose 11-letter CONTEXT centre base equals Reference_Allele >= 0.999.
#     Reported: trinucleotide from the transcript against the MAF CONTEXT (equal, opposite strand, flank differs = exon edge, centre differs).
#   - By-product for A3 (SqnPAN_16): every chosen-file missense row is mapped to (case, accession, position, ref, alt) with the SqnPAN_06 logic (copied) and gets its channel;
#     the key set is compared with the chosen-file rows of P2_06_mutations.tsv.gz (reported; SqnPAN_16 gates its own coverage).
# Section 2, spectra and weights:
#   - Spectrum per cancer type: coding SNVs (Missense_Mutation, Silent, Nonsense_Mutation; Variant_Type SNP) in the chosen file of every case, on genes that have a mapped
#     transcript; channel = pyrimidine-centred 96 classes from CONTEXT[4], Reference_Allele, Tumor_Seq_Allele2, CONTEXT[6]. Rate of channel k in cancer type c =
#     (count + 0.5) / number of CDS positions of the mapped transcripts with the same trinucleotide (either strand). Set "nohyper": spectra without hypermutated cases.
#   - Expected rate of a position = sum of the rates of the single-nucleotide changes of its codon that make a counted substitution (N: any missense; S/T: missense to
#     a residue that is not S or T); flanks from the transcript. Exposure of a stratum class = sum over its positions; pooled scope = sum over cancer types.
#   - GATE (spectrum sanity, a check of channel orientation): SKCM C>T share >= 0.60; LUAD and LUSC C>A share >= 0.20.
#   - Calibration control (judged): background positions with codon opportunities, split inside each protein x residue stratum at the median pooled expected rate
#     (above = "high"); the Q2 code compares high against low with the new exposure (main model, pooled, N + S/T, seed sub-base 32780000). PASSED if its 90% CI lies
#     inside [-0.32, 0.32] log2. The same split with position counts is reported (shows how much context matters). If FAILED, every A1 final label gets the suffix
#     _calibration_FAILED (results descriptive).
# Test (sections 3 and 4): the SqnPAN_08 / SqnPAN_09 code unchanged except the exposure (Engine in SqnPAN_14_context_helpers.py): Mantel-Haenszel rate ratio, strata protein x residue
#   (strict: also pLDDT and exposure bins), 2,000 bootstrap resamples over proteins, tested-unit rule (>= 20 informative proteins, expected glycosite count >= 10), labels
#   (Changed = BH < 0.05 and |log2 RR| >= 0.32; Equivalent = 90% CI inside +-0.14), BH families, final-label rule 1A of SqnPAN_09. Units: the 144 units of SqnPAN_08 / SqnPAN_09
#   (main pooled 9 + main per cancer type T1 99; strict, nohyper, union, norescue pooled 9 each) with the new exposure, plus 9 pooled units "uniform_restricted"
#   (position counts on the same restricted proteins, to separate the effect of the restriction from the effect of context). 153 units.
# Controls (section 3): C-S: the Engine in uniform mode on the full SqnPAN_08 positions must reproduce P2_08_C2C3_run1.tsv repeats 1 and 2 line for line (GATE).
#   C2 natural null and C3 planted change as in SqnPAN_08 (20 repeats; runs 1 and 2 must be identical; run 3 other seed): capability and Equivalent-allowed per unit.
# Seeds (detail fixed while writing, before any data: sub-bases inside the block reserved for A1, 32700000 to 32799999, checked in the code): C2/C3 run1 and run2
#   32700000, run3 32720000; real runs 1 and 2 32740000, run 3 32760000; calibration 32780000.
# Not blind: Fan et al. 2018 reported the position-3 excess; Q2 primary already read (SqnPAN_09). Reported beside the primary in every text.
# Usage (from the repository folder):  python scripts/SqnPAN_14_context_q2.py all [workers=6]
#   "all" runs sections 1 to 4 and stops at the first section that does not finish; "3" restarts from section 3 (C2/C3 repeats already saved are skipped).
# Logs: results\logs\P2_14_section<k>_<name>.txt. Outputs listed in each section.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, csv, gzip, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_06_helpers as P
import SqnPAN_08_selection_controls as S8
import SqnPAN_14_context_helpers as C

ARGS = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
WORKERS = int(ARGS.get("workers", "6"))
BLOCK = 32700000
SEEDS = dict(cc1=32700000, cc3=32720000, real1=32740000, real3=32760000, calib=32780000)
CTRL_BASES = {1: SEEDS["cc1"], 2: SEEDS["cc1"], 3: SEEDS["cc3"]}; REAL_BASES = {1: SEEDS["real1"], 2: SEEDS["real1"], 3: SEEDS["real3"]}
CODON_MIN = 0.995; CONTEXT_MIN = 0.999; CAL_LIM = 0.32
SANITY = [("TCGA-SKCM", "C>T", 0.60), ("TCGA-LUAD", "C>A", 0.20), ("TCGA-LUSC", "C>A", 0.20)]

OK1 = os.path.join(H.PROC_OUT, "P2_14_section1_ok.txt"); OK2 = os.path.join(H.PROC_OUT, "P2_14_section2_ok.txt")
POS_NPZ = os.path.join(H.PROC_OUT, "P2_14_positions.npz"); ENGINE_NPZ = os.path.join(H.PROC_OUT, "P2_14_engine.npz")
CASE_SPEC = os.path.join(H.PROC_OUT, "P2_14_case_spectra.tsv"); A3_TABLE = os.path.join(H.PROC_OUT, "P2_14_a3_mutation_channels.tsv.gz")
TMAP = os.path.join(H.TABLES, "P2_14_transcript_map.tsv"); MISM = os.path.join(H.TABLES, "P2_14_codon_mismatches.tsv")
SPECTRA = os.path.join(H.TABLES, "P2_14_spectra.tsv"); CALIB = os.path.join(H.TABLES, "P2_14_calibration.tsv")
CALIB_VERDICT = os.path.join(H.PROC_OUT, "P2_14_calibration_verdict.txt"); VERDICT = os.path.join(H.PROC_OUT, "P2_14_controls_verdict.txt")
PARTS = os.path.join(H.PROC_OUT, "P2_14_parts"); CAPAB = os.path.join(H.TABLES, "P2_14_unit_capability.tsv")
RESULTS = os.path.join(H.TABLES, "P2_14_Q2ctx_results.tsv")
def CTRL_TABLE(k): return os.path.join(H.TABLES, f"P2_14_C2C3_run{k}.tsv")
def REAL_TABLE(k): return os.path.join(H.TABLES, f"P2_14_Q2ctx_run{k}.tsv")

CONFIGS_A1 = [dict(name="main", kind="main", mask="main", mode="ctx", spec="main", cancer=True),
              dict(name="strict", kind="strict", mask="main", mode="ctx", spec="main", cancer=False),
              dict(name="nohyper", kind="main", mask="nohyper", mode="ctx", spec="nohyper", cancer=False),
              dict(name="union", kind="main", mask="union", mode="ctx", spec="main", cancer=False),
              dict(name="norescue", kind="main", mask="norescue", mode="ctx", spec="main", cancer=False),
              dict(name="uniform_restricted", kind="main", mask="main", mode="uniform", spec=None, cancer=False)]
CFG_NAMES = [c["name"] for c in CONFIGS_A1]
CONFIGS_P208 = [dict(name=n, kind="strict" if n == "strict" else "main", mask="main" if n == "strict" else n, mode="uniform", spec=None, cancer=(n == "main"))
                for n in S8.CONFIGS]
CLS = {"M": "Missense_Mutation", "S": "Silent", "X": "Nonsense_Mutation"}

def need_ok(path, what):
    if not os.path.exists(path) or not open(path).read().startswith("OK"): raise RuntimeError(f"{what} has not finished ({os.path.basename(path)} missing)")

def check_md5(log, pairs):
    for path, rec in pairs:
        got = H.md5_file(path); same = got == rec
        log(f"md5 {os.path.basename(path)}: {got} ({'same as recorded' if same else 'DIFFERENT, recorded ' + rec})")
        if not same: raise RuntimeError("an input md5 differs from the recorded value: stop and decide")

# ======================= Section 1: structure check =======================
def section1(log):
    import multiprocessing as mp
    log("Pre-specification: header of this script (fixed 2026-10-31). Post hoc, added after the primary analysis.")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__} | workers {WORKERS}")
    check_md5(log, list(S8.MD5.items()) + [(C.GENCODE, C.GENCODE_MD5), (C.NPZ8, C.NPZ8_MD5)])
    t0 = time.time(); U = H.load_uniprot_main(); G, ntx, nbad = C.load_gencode()
    log(f"GENCODE v36: {ntx} records, {len(G)} IDs without version, records without one CDS field {nbad}, IDs with more than one version {sum(1 for v in G.values() if len(v) > 1)} | {time.time() - t0:.0f} s")
    M, mc = C.build_transcript_map(U, G)
    for k in sorted(mc): log(f"  {k}: {mc[k]}")
    with open(TMAP, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["accession", "gene", "transcript", "version", "gene_id", "n_candidates", "listed_by_uniprot"])
        for acc in sorted(M):
            tid, vi, nc, li = M[acc]; d = G[tid][vi]; w.writerow([acc, U[acc]["gene"], tid, d["ver"], d["gene"], nc, int(li)])
    accs, P_prot, P_res, P_role, P_posn = C.rebuild_positions(U); D8 = C.load_npz(C.NPZ8)
    same = list(D8["accs"]) == accs and all(np.array_equal(a, D8[k]) for a, k in ((P_prot, "P_prot"), (P_res, "P_res"), (P_role, "P_role")))
    log(f"SqnPAN_08 position list rebuilt: {len(P_prot)} positions | equal to P2_08_inputs.npz (accessions, P_prot, P_res, P_role): {'yes' if same else 'NO'}")
    if not same: raise RuntimeError("rebuilt position list differs from SqnPAN_08: stop and decide")
    t0 = time.time(); OPP, P_ctx, oc = C.opportunities(accs, P_prot, P_res, P_posn, U, G, M)
    for k in sorted(oc): log(f"  opportunities: {k}: {oc[k]}")
    if oc["codon does not code the UniProt residue (must be 0)"]: raise RuntimeError("a mapped codon does not code the UniProt residue")
    tri_opp = C.cds_trinucleotides(G, M); log(f"opportunities and CDS trinucleotides: {time.time() - t0:.0f} s | CDS trinucleotides counted {int(tri_opp.sum())}")
    mapped_genes = {G[v[0]][v[1]]["gene"] for v in M.values()}
    # ---- MAF pass ----
    revmap = P.build_revmap(U)
    chosen = {}
    with open(S8.CASES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"): chosen[r["chosen_file"]] = (r["case_id"], r["project"], r["hypermutated"])
    with open(H.MANIFEST, newline="", encoding="utf-8") as f: man = sorted((r["id"], r["filename"]) for r in csv.DictReader(f, delimiter="\t"))
    jobs = [(fid, os.path.join(H.MAF_DIR, fid, fn), fid in chosen) for fid, fn in man]
    log(f"MAF files {len(jobs)} | chosen files {sum(1 for j in jobs if j[2])} | cases {len(chosen)}")
    T = collections.Counter(); fdates = collections.Counter(); spec = {v[0]: np.zeros(96, np.int64) for v in chosen.values()}; mism = []; a3 = {}
    t0 = time.time(); done = 0
    pool = mp.get_context("spawn").Pool(max(1, WORKERS + 2))
    try:
        for res in pool.imap_unordered(C.maf_worker, jobs, chunksize=4):
            if "error" in res: raise RuntimeError(f"file {res['fid']}: {res['error']}")
            done += 1; T.update(res["tally"]); fdates[res["fdate"]] += 1; info = chosen.get(res["fid"])
            for (c, vt, tid, ppos, aas, codons, cdspos, ctx, refal, alt2, sw, hugo, ae) in res["rows"]:
                cl = CLS[c]; T[f"coding rows: {cl}"] += 1; ctxu = ctx.upper(); ctx_ok = len(ctxu) == 11 and ctxu[5] == refal
                if info is not None and c == "M":
                    key, how = C.map_missense(sw, ppos, aas, hugo, ae, U, revmap); T[f"A3: chosen-file missense rows {how}"] += 1
                    if key is not None:
                        ch = C.channel(ctxu[4], refal, alt2, ctxu[6]) if (vt == "SNP" and ctx_ok) else (-2 if vt in ("DNP", "TNP", "ONP") else -1)
                        ks = f"{info[0]}|{key[0]}|{key[1]}|{key[2]}|{key[3]}"; old = a3.get(ks)
                        if old is None: a3[ks] = (ch, how)
                        elif old[0] != ch: a3[ks] = (-3, old[1])
                vl = G.get(tid.split(".")[0])
                if not vl: T[f"rows on transcripts absent from GENCODE v36 pc file: {cl}"] += 1; continue
                if c == "M":
                    try:
                        p = int(ppos.split("/")[0].split("-")[0]); refc = codons.split("/")[0].upper()
                        if len(refc) != 3 or p < 1: raise ValueError
                    except Exception:
                        T["codon test: missense rows unparsable"] += 1; p = None
                    if p is not None:
                        T["codon test: parsed missense rows on GENCODE transcripts"] += 1
                        if len(vl) > 1: T["codon test: rows on IDs with more than one version"] += 1
                        hits = [3 * p <= d["ce"] - d["cs"] and d["seq"][d["cs"] + 3 * p - 3:d["cs"] + 3 * p] == refc for d in vl]
                        if hits[0]: T["codon test: match (first version)"] += 1
                        elif any(hits): T["codon test: match only with another version"] += 1; T["codon test: version used " + vl[hits.index(True)]["ver"]] += 1
                        else:
                            d = vl[0]; L_ = d["ce"] - d["cs"]
                            reason = "protein position beyond the CDS" if 3 * p > L_ else ("CDS length not a multiple of 3" if L_ % 3 else "codon differs")
                            T["codon test: mismatch, " + reason] += 1
                            mism.append([res["fid"], tid, ";".join(x["ver"] for x in vl), p, refc, ";".join(x["seq"][x["cs"] + 3 * p - 3:x["cs"] + 3 * p] for x in vl), reason, hugo, aas])
                if vt == "SNP":
                    d = vl[0]
                    try: cp = int(cdspos.split("/")[0].split("-")[0])
                    except Exception: cp = None
                    if cp is not None and ctx_ok:
                        idx = d["cs"] + cp - 1
                        if 1 <= idx < len(d["seq"]) - 1:
                            tt = d["seq"][idx - 1:idx + 2]; tg = ctxu[4:7]
                            if tt == tg: T["context check: transcript trinucleotide equals CONTEXT (same strand)"] += 1
                            elif tt == C.revcomp(tg): T["context check: equals reverse complement of CONTEXT (opposite strand)"] += 1
                            elif tt[1] in (tg[1], C.comp(tg[1])): T["context check: centre equal, a flank differs (exon edge expected)"] += 1
                            else: T["context check: centre base differs"] += 1
                    if info is not None:
                        if d["gene"] not in mapped_genes: T["spectrum: chosen-file coding SNVs on genes without a mapped transcript (left out)"] += 1
                        elif not ctx_ok: T["spectrum: chosen-file coding SNVs with unusable CONTEXT (left out)"] += 1
                        else:
                            ch = C.channel(ctxu[4], refal, alt2, ctxu[6])
                            if ch < 0: T["spectrum: chosen-file coding SNVs with a non-ACGT base (left out)"] += 1
                            else: spec[info[0]][ch] += 1; T["spectrum: chosen-file coding SNVs counted"] += 1
            if done % 1000 == 0: log(f"  MAF files done {done}/{len(jobs)} | elapsed {time.time() - t0:.0f} s")
    finally:
        pool.close(); pool.join()
    log(f"MAF pass: {done} files | {time.time() - t0:.0f} s")
    for k in sorted(fdates): log(f"  #filedate {k}: {fdates[k]} files")
    for k in sorted(T): log(f"  {k}: {T[k]}")
    parsed = T["codon test: parsed missense rows on GENCODE transcripts"]; matched = T["codon test: match (first version)"] + T["codon test: match only with another version"]
    share = matched / parsed if parsed else 0.0
    snv = T["SNV rows"]; cshare = T["SNV rows whose CONTEXT centre base equals Reference_Allele"] / snv if snv else 0.0
    log(f"GATE codon test: matched {matched} of {parsed} = {share:.5f} (limit {CODON_MIN}) | GATE CONTEXT centre: {cshare:.5f} of {snv} SNV rows (limit {CONTEXT_MIN})")
    with open(MISM, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["file_id", "transcript", "gencode_versions", "protein_position", "maf_ref_codon", "gencode_codons", "reason", "Hugo_Symbol", "Amino_acids"]); w.writerows(mism)
    log(f"codon mismatches written: {len(mism)} rows | {MISM}")
    for r in mism[:40]: log("  mismatch: " + " | ".join(str(x) for x in r))
    if share < CODON_MIN or cshare < CONTEXT_MIN: raise RuntimeError("structure gate failed (codon test or CONTEXT): stop and decide")
    # ---- outputs ----
    projh = {v[0]: (v[1], v[2]) for v in chosen.values()}
    with open(CASE_SPEC, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["case_id", "project", "hypermutated", "n_snv"] + C.CHANNELS)
        for case in sorted(spec): w.writerow([case, projh[case][0], projh[case][1], int(spec[case].sum())] + [int(x) for x in spec[case]])
    nm = collections.Counter()
    with gzip.open(A3_TABLE, "wt", newline="", encoding="utf-8", compresslevel=6) as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["case_id", "accession", "position", "ref", "alt", "channel", "how"])
        for ks in sorted(a3):
            ch, how = a3[ks]; name = C.CHANNELS[ch] if ch >= 0 else {-1: "invalid", -2: "multi", -3: "ambiguous"}[ch]; nm[name if ch < 0 else "single channel"] += 1
            w.writerow(ks.split("|") + [name, how])
    log(f"A3 table {A3_TABLE}: {len(a3)} keys | channels {dict(nm)}")
    p06 = set()
    with gzip.open(S8.MUTS, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["in_chosen_file"] == "1": p06.add(f"{r['case_id']}|{r['accession']}|{r['position']}|{r['ref']}|{r['alt']}")
    ka = set(a3); log(f"A3 keys against SqnPAN_06 chosen-file rows: SqnPAN_06 {len(p06)} | this pass {len(ka)} | both {len(ka & p06)} | only SqnPAN_06 {len(p06 - ka)} | only this pass {len(ka - p06)}")
    np.savez_compressed(POS_NPZ, accs=np.array(accs), P_posn=P_posn, P_ctx=P_ctx, OPP=OPP, tri_opp=tri_opp)
    out_md5 = {os.path.basename(p): H.md5_file(p) for p in (POS_NPZ, CASE_SPEC, A3_TABLE, TMAP, MISM)}
    for k, v in out_md5.items(): log(f"output {k} md5 {v}")
    with open(OK1, "w") as o: o.write("OK\n" + "".join(f"{k}\t{v}\n" for k, v in out_md5.items()))

# ======================= Section 2: spectra, weights, calibration =======================
def section2(log):
    need_ok(OK1, "section 1"); check_md5(log, [(C.NPZ8, C.NPZ8_MD5)])
    pos = C.load_npz(POS_NPZ); D8 = C.load_npz(C.NPZ8); projects = [str(x) for x in D8["projects"]]; pj = {p: i for i, p in enumerate(projects)}
    cnt = {s: np.zeros((len(projects), 96)) for s in ("main", "nohyper")}; ncase = {s: np.zeros(len(projects), int) for s in ("main", "nohyper")}
    with open(CASE_SPEC, newline="", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t"); hdr = next(rd)
        if hdr[4:] != C.CHANNELS: raise RuntimeError("channel columns of the case table differ")
        for r in rd:
            i = pj[r[1]]; v = np.array([float(x) for x in r[4:]])
            cnt["main"][i] += v; ncase["main"][i] += 1
            if r[2] == "0": cnt["nohyper"][i] += v; ncase["nohyper"][i] += 1
    opp_ch = pos["tri_opp"][[C.TRI_IDX[t] for t in C.CH_TRI]]
    if (opp_ch <= 0).any(): raise RuntimeError("a trinucleotide has no CDS opportunity")
    R = {}; rows = []
    for s in ("main", "nohyper"):
        r = (cnt[s] + 0.5) / opp_ch[None, :]; R[s] = np.column_stack([r.sum(axis=0), r.T])
        for i, p in enumerate(projects):
            for k in range(96): rows.append([s, p, ncase[s][i], C.CHANNELS[k], int(cnt[s][i, k]), f"{r[i, k]:.6g}"])
    with open(SPECTRA, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["set", "project", "n_cases", "channel", "count", "rate"]); w.writerows(rows)
    log(f"spectra table {SPECTRA} | md5 {H.md5_file(SPECTRA)}")
    subs = sorted(set(C.CH_SUB)); bad = []
    log("substitution-class shares of the coding SNV counts (set main): " + " ".join(subs))
    for i, p in enumerate(projects):
        tot = cnt["main"][i].sum(); sh = {u: (cnt["main"][i][[k for k in range(96) if C.CH_SUB[k] == u]].sum() / tot if tot else 0.0) for u in subs}
        log(f"  {p} cases {ncase['main'][i]} SNVs {int(tot)}: " + " ".join(f"{sh[u]:.3f}" for u in subs))
        for proj, u, lim in SANITY:
            if p == proj: log(f"  sanity {p} {u} share {sh[u]:.3f} (limit >= {lim}): {'ok' if sh[u] >= lim else 'FAILED'}"); bad += [p] if sh[u] < lim else []
    if bad: raise RuntimeError("spectrum sanity failed (channel orientation?): stop and decide")
    E = dict(D8); E["P_ctx"] = pos["P_ctx"]; E["OPP"] = pos["OPP"]; E["R_main"] = R["main"]; E["R_nohyper"] = R["nohyper"]
    np.savez_compressed(ENGINE_NPZ, **E); log(f"engine inputs {ENGINE_NPZ} | md5 {H.md5_file(ENGINE_NPZ)}")
    ctx = E["P_ctx"]; role = E["P_role"]; res = E["P_res"]; evp = E["ev_pos"][E["ev_main"]]
    log("--- coverage (positions; main-model events): class: all / with codon opportunities ---")
    for t in (0, 1, 2, 3):
        for rc, rn in ((0, "N"), (1, "S/T")):
            mk = (role == t) & (res == rc); ev_all = int(mk[evp].sum()); ev_ctx = int((mk & ctx)[evp].sum())
            log(f"  {'background' if t == 0 else 'glycosite T' + str(t)} {rn}: positions {int(mk.sum())} / {int((mk & ctx).sum())} | events {ev_all} / {ev_ctx}")
    # ---- calibration control ----
    C.check_seed_blocks(log, BLOCK, SEEDS.values())
    O = E["OPP"].astype(np.int64); wp = np.where(O >= 0, R["main"][np.maximum(O, 0), 0], 0.0).sum(axis=1)
    strat = E["P_prot"].astype(np.int64) * 2 + res; idx = np.flatnonzero((role == 0) & ctx)
    order = np.lexsort((wp[idx], strat[idx])); si = strat[idx][order]; wi = wp[idx][order]
    uq, start, cn = np.unique(si, return_index=True, return_counts=True); med = 0.5 * (wi[start + (cn - 1) // 2] + wi[start + cn // 2])
    high = wp[idx] > med[np.searchsorted(uq, strat[idx])]
    role_cal = np.full(len(role), -1, np.int8); role_cal[idx] = np.where(high, 1, 0)
    log(f"calibration split: background positions with opportunities {len(idx)} | high {int(high.sum())} | low {int((~high).sum())}")
    cfgs = [dict(name="calib_ctx", kind="main", mask="main", mode="ctx", spec="main", cancer=True), dict(name="calib_uniform", kind="main", mask="main", mode="uniform", spec=None, cancer=False)]
    eng = C.Engine(E, cfgs, restricted=True, tiers=(1,), role=role_cal)
    recs = eng.evaluate(eng.real_labels(), np.random.default_rng(np.random.SeedSequence(SEEDS["calib"])), None, C.LEVEL1, extras=True)
    with open(CALIB, "w", newline="", encoding="utf-8") as o: o.write("\t".join(C.COLS_REAL) + "\n" + "".join(C.real_line(0, r) + "\n" for r in recs))
    log(f"calibration table {CALIB} | md5 {H.md5_file(CALIB)} (glyco = high-rate half, background = low-rate half)")
    for r in recs:
        if r["scope"] == "pooled" or r["group"] == "both":
            log(f"  {r['cfg']} {r['group']} {r['scope']}: log2RR {r['est']:.3f} (90% CI {r['lo90']:.3f} to {r['hi90']:.3f}) events high {r['ev_g']} low {r['ev_b']} tested {r['tested']}")
    pr = [r for r in recs if r["cfg"] == "calib_ctx" and r["group"] == "both" and r["scope"] == "pooled"][0]
    ok = (not np.isnan(pr["lo90"])) and pr["lo90"] >= -CAL_LIM and pr["hi90"] <= CAL_LIM
    with open(CALIB_VERDICT, "w") as o: o.write(("PASSED" if ok else "FAILED") + f"\nmain, pooled, N + S/T, context exposure: log2RR {pr['est']:.4f}, 90% CI {pr['lo90']:.4f} to {pr['hi90']:.4f}, limit +-{CAL_LIM}\nwritten {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    log(f"CALIBRATION {'PASSED' if ok else 'FAILED (every A1 final label gets the suffix _calibration_FAILED)'}")
    with open(OK2, "w") as o: o.write(f"OK\n{os.path.basename(ENGINE_NPZ)}\t{H.md5_file(ENGINE_NPZ)}\n")

# ======================= Section 3: controls C-S, C2, C3 =======================
def section3(log):
    need_ok(OK2, "section 2"); C.check_seed_blocks(log, BLOCK, SEEDS.values())
    check_md5(log, [(C.RUN1_P208, C.RUN1_P208_MD5), (C.NPZ8, C.NPZ8_MD5)])
    t0 = time.time(); eng = C.Engine(C.load_npz(C.NPZ8), CONFIGS_P208, restricted=False)
    mine = C.control_repeat(eng, S8.RUN_BASE[1], 1) + C.control_repeat(eng, S8.RUN_BASE[1], 2); del eng
    with open(C.RUN1_P208, encoding="utf-8") as f: ref = [ln.rstrip("\n") for ln in f if ln.split("\t", 1)[0] in ("1", "2")]
    eq = sum(1 for a, b in zip(mine, ref) if a == b); cs_ok = len(mine) == len(ref) and eq == len(ref)
    log(f"C-S: Engine (uniform, full SqnPAN_08 positions, SqnPAN_08 seeds) reproduced {eq} of {len(ref)} lines of P2_08_C2C3_run1.tsv repeats 1 and 2 | {time.time() - t0:.0f} s | {'PASSED' if cs_ok else 'FAILED'}")
    if not cs_ok:
        for a, b in list(zip(mine, ref))[:200]:
            if a != b: log("  first differing line, Engine: " + a); log("  first differing line, SqnPAN_08: " + b); break
        raise RuntimeError("control C-S failed: the Engine does not reproduce SqnPAN_08: stop and decide")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__} | workers {WORKERS} | repeats {S8.N_REP} | bootstrap {S8.B}")
    C.run_controls(log, ENGINE_NPZ, CONFIGS_A1, True, CTRL_BASES, PARTS, CTRL_TABLE, WORKERS)
    cap, c2ok, same12 = C.capability(log, CTRL_TABLE, CFG_NAMES, CAPAB)
    pk = ("main", "T1", "both", "pooled")
    lines = ["PASSED" if (cs_ok and same12) else "FAILED", f"C-S: {'PASSED' if cs_ok else 'FAILED'}", f"C2/C3 run1 and run2 identical: {'yes' if same12 else 'no'}",
             *[f"C2 run{k} {cfg}: {'PASSED' if c2ok[(k, cfg)] else 'FAILED'}" for k in (1, 3) for cfg in CFG_NAMES],
             f"primary unit (main, T1, both, pooled): {cap[pk][0]}; Equivalent allowed: {'yes' if cap[pk][1] else 'no'}", f"written {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    with open(VERDICT, "w") as o: o.write("\n".join(lines) + "\n")
    for x in lines: log("VERDICT " + x)

# ======================= Section 4: real result =======================
def section4(log):
    if not os.path.exists(VERDICT) or open(VERDICT).readline().strip() != "PASSED": raise RuntimeError("controls verdict is not PASSED: no real result is computed")
    calib = open(CALIB_VERDICT).readline().strip(); suffix = "" if calib == "PASSED" else "_calibration_FAILED"
    log(f"controls verdict PASSED | calibration {calib} | label suffix '{suffix}'")
    C.check_seed_blocks(log, BLOCK, SEEDS.values()); check_md5(log, [(C.RES_P209, C.RES_P209_MD5)])
    eng = C.Engine(C.load_npz(ENGINE_NPZ), CONFIGS_A1, restricted=True)
    C.real_runs(log, eng, REAL_BASES, REAL_TABLE)
    p209 = C.read_run(C.RES_P209); cap = C.read_cap(CAPAB)
    def extra(key):
        r = p209.get(key); return [calib, r["log2rr"] if r else "", r["final_label"] if r else ""]
    runs, fin = C.assemble(log, REAL_TABLE, cap, CFG_NAMES, RESULTS, suffix, ("calibration", "p209_log2rr", "p209_final_label"), extra)
    pk = ("main", "T1", "both", "pooled")
    log("=== PRIMARY UNIT with context-aware exposure (main, T1, N + S/T, pooled); added after the primary analysis ===")
    log("  A1:    " + C.unit_line(runs, fin, cap, pk))
    log(f"  SqnPAN_09 (primary, position counts, all proteins): log2RR {C.ff(p209[pk]['log2rr'])} (90% CI {C.ff(p209[pk]['lo90'])} to {C.ff(p209[pk]['hi90'])}) FINAL {p209[pk]['final_label']}")
    kr = ("uniform_restricted", "T1", "both", "pooled"); log("  uniform_restricted (position counts, same restricted proteins): " + C.unit_line(runs, fin, cap, kr))
    log("=== POOLED UNITS, all set-ups ===")
    for cfg in CFG_NAMES:
        for t in ("T1", "T2", "T3"):
            for g in ("N", "ST", "both"):
                k = (cfg, t, g, "pooled"); r9 = p209.get(k)
                log(f"  {cfg} {t} {g}: " + C.unit_line(runs, fin, cap, k) + (f" | SqnPAN_09 {C.ff(r9['log2rr'])} {r9['final_label']}" if r9 else ""))
    cnt = collections.Counter(fin[1][k][0] for k in cap if k[3] != "pooled")
    log(f"=== per cancer type (main, T1; 99 units): final labels {dict(cnt)} ===")
    for k in sorted(cap):
        if k[3] != "pooled" and fin[1][k][0] not in (C.NOT_TESTED, C.NOTSENS, "Inconclusive" + suffix, C.NOTSENS + suffix): log(f"  {k[3]} {k[2]}: " + C.unit_line(runs, fin, cap, k))

SECTIONS = {"1": ("structure_check", section1), "2": ("spectra_weights_calibration", section2), "3": ("controls_CS_C2_C3", section3), "4": ("real_result", section4)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    C.run_chain("P2_14", SECTIONS, "1" if which == "all" else which)
# ---- END OF SCRIPT ----
