# SqnPAN_23_dndscv_prep.py - dNdScv addition, section 1: input files (added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Section 1 of the dNdScv addition. The pre-specification is in the header of SqnPAN_23_dndscv.R (fixed 2026-11-09,
# D1 A, D2 A, D3 B). This script only builds inputs; no mutation-at-site count or result is computed here.
#   data\processed\P2_23_mutations.tsv.gz  sampleID chr pos ref mut: the SNVs of the chosen file of every case (rule 4A, P2_06_cases.tsv),
#                                          one row per case and change (a duplicate row in a file is counted once); indels and
#                                          DNP/TNP/ONP rows are counted and left out
#   data\processed\P2_23_cases.tsv         case_id project hypermutated n_snv
#   data\processed\P2_23_positions.tsv.gz  accession gene position residue role: the SqnPAN_08 position list (role 1, 2, 3 = sequon N or S/T of
#                                          that tier; 0 = N or S/T in no sequon), rebuilt by P2_14_context_helpers.rebuild_positions
#   data\processed\P2_23_uniprot.tsv       accession gene sequence (UniProt reviewed 2026_03; proteins of the position list)
#   data\processed\P2_23_bailey.txt        the 299 Bailey driver genes after the SqnPAN_20 alias step
#   data\processed\P2_23_section1_ok.txt   gate file: md5 of the five outputs (checked by SqnPAN_23_dndscv.R)
# Run (from the repository folder):  python scripts/SqnPAN_23_dndscv_prep.py
# Log: results\logs\P2_23_section1_inputs.txt (last line starts with "---- END OF SECTION 1 ----")
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, csv, gzip, glob, time, collections
import multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
import SqnPAN_05_helpers as H
import SqnPAN_08_selection_controls as S8
import SqnPAN_14_context_helpers as C
import SqnPAN_20_hotspots as P20

OUT_MUT = os.path.join(H.PROC, "P2_23_mutations.tsv.gz"); OUT_CASES = os.path.join(H.PROC, "P2_23_cases.tsv")
OUT_POS = os.path.join(H.PROC, "P2_23_positions.tsv.gz"); OUT_UNI = os.path.join(H.PROC, "P2_23_uniprot.tsv")
OUT_BAIL = os.path.join(H.PROC, "P2_23_bailey.txt"); GATE = os.path.join(H.PROC, "P2_23_section1_ok.txt")
NEED = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type"]
N_POS_P208 = 1965194; WORKERS = 8

def read_chosen(job):
    cid, path = job
    t = collections.Counter(); snv = set(); hdr = None; ix = None
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        for line in f:
            if line[0] == "#": continue
            if hdr is None:
                hdr = line.rstrip("\n").split("\t"); miss = [k for k in NEED if k not in hdr]
                if miss: return dict(cid=cid, error="missing columns: " + ", ".join(miss))
                ix = [hdr.index(k) for k in NEED]; continue
            v = line.rstrip("\n").split("\t"); vt = v[ix[4]]; t["rows with Variant_Type " + vt] += 1
            if vt != "SNP": continue
            ref = v[ix[2]]; alt = v[ix[3]]
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                t["SNP rows left out (alleles not one of A, C, G, T, or equal)"] += 1; continue
            key = (v[ix[0]], int(v[ix[1]]), ref, alt)
            if key in snv: t["SNP rows repeated in the same file (counted once)"] += 1; continue
            snv.add(key)
    return dict(cid=cid, tally=t, snv=sorted(snv))

def section1(log):
    t0 = time.time()
    log("SqnPAN_23 section 1: inputs for dNdScv (added after the primary analysis). Python", sys.version.split()[0])
    for p in (H.UNIPROT_MAIN, S8.SITES, S8.CASES):
        m = H.md5_file(p); exp = S8.MD5.get(p)
        log(f"md5 {os.path.basename(p)} {m}" + ("" if exp is None else (" (equal to the SqnPAN_08 value)" if m == exp else f" DIFFERS from the SqnPAN_08 value {exp}")))
        if exp is not None and m != exp: raise RuntimeError("an input md5 differs from SqnPAN_08: " + p)
    cases = []
    with open(S8.CASES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"): cases.append((r["case_id"], r["project"], r["chosen_file"], r["hypermutated"]))
    hv = collections.Counter(c[3] for c in cases)
    log("cases", len(cases), "; projects", len({c[1] for c in cases}), "; hypermutated values", dict(hv))
    if set(hv) - {"0", "1"}: raise RuntimeError("unexpected values in the hypermutated column")
    if len({c[0] for c in cases}) != len(cases): raise RuntimeError("a case_id occurs twice in P2_06_cases.tsv")
    jobs, missing = [], []
    for cid, proj, fid, hyp in cases:
        g = glob.glob(os.path.join(H.MAF_DIR, fid, "*.maf.gz"))
        if len(g) == 1: jobs.append((cid, g[0]))
        else: missing.append((cid, fid, len(g)))
    log("chosen MAF files found", len(jobs), "; not found or not unique", len(missing))
    if missing: raise RuntimeError(f"chosen MAF file not found for {len(missing)} cases, first {missing[:3]}")
    tally = collections.Counter(); per = {}; done = 0; errs = []
    with mp.Pool(WORKERS) as pool:
        for res in pool.imap_unordered(read_chosen, jobs, chunksize=8):
            done += 1
            if "error" in res: errs.append((res["cid"], res["error"])); continue
            tally.update(res["tally"]); per[res["cid"]] = res["snv"]
            if done % 1000 == 0: log(f"  files read {done} of {len(jobs)} ({round(time.time() - t0)} s)")
    if errs: raise RuntimeError(f"{len(errs)} files with errors, first {errs[:3]}")
    log("rows of the chosen files:")
    for k in sorted(tally): log(f"  {k}: {tally[k]}")
    chrc = collections.Counter(); nsnv = {}
    with gzip.open(OUT_MUT, "wt", encoding="utf-8", newline="") as f:
        f.write("sampleID\tchr\tpos\tref\tmut\n")
        for cid, proj, fid, hyp in cases:
            rows = per[cid]; nsnv[cid] = len(rows)
            for ch, pos, ref, alt in rows:
                chrc[ch] += 1; f.write(f"{cid}\t{ch}\t{pos}\t{ref}\t{alt}\n")
    log("SNV rows written", sum(nsnv.values()), "; cases with no SNV", sum(1 for v in nsnv.values() if v == 0))
    log("chromosome names:", " ".join(f"{k}:{v}" for k, v in sorted(chrc.items())))
    byp = collections.Counter()
    with open(OUT_CASES, "w", encoding="utf-8", newline="") as f:
        f.write("case_id\tproject\thypermutated\tn_snv\n")
        for cid, proj, fid, hyp in cases:
            byp[proj] += nsnv[cid]; f.write(f"{cid}\t{proj}\t{hyp}\t{nsnv[cid]}\n")
    log("SNVs per project (largest first):", "; ".join(f"{k} {v}" for k, v in byp.most_common()))
    U = H.load_uniprot_main(); accs, P_prot, P_res, P_role, P_posn = C.rebuild_positions(U)
    log("UniProt entries", len(U), "; positions", len(P_posn), "(expected 1,965,194 as in SqnPAN_08)")
    if len(P_posn) != N_POS_P208: raise RuntimeError("the position list differs from SqnPAN_08")
    cnt = collections.Counter(); used = set()
    with gzip.open(OUT_POS, "wt", encoding="utf-8", newline="") as f:
        f.write("accession\tgene\tposition\tresidue\trole\n")
        for i in range(len(P_posn)):
            a = accs[int(P_prot[i])]; p = int(P_posn[i]); r = U[a]["seq"][p - 1]; ro = int(P_role[i])
            if (r == "N") != (int(P_res[i]) == 0) or r not in "NST": raise RuntimeError(f"residue check failed at {a} {p}")
            cnt[(ro, "N" if r == "N" else "S/T")] += 1; used.add(a)
            f.write(f"{a}\t{U[a]['gene']}\t{p}\t{r}\t{ro}\n")
    for k in sorted(cnt): log(f"  positions role {k[0]} {k[1]}: {cnt[k]}")
    with open(OUT_UNI, "w", encoding="utf-8", newline="") as f:
        f.write("accession\tgene\tseq\n")
        for a in sorted(used): f.write(f"{a}\t{U[a]['gene']}\t{U[a]['seq']}\n")
    log("proteins written", len(used), "; without a primary gene name", sum(1 for a in used if not U[a]["gene"]))
    drv = P20.bailey_genes(log, U)
    with open(OUT_BAIL, "w", encoding="utf-8", newline="") as f:
        for g in sorted(drv): f.write(g + "\n")
    lines = []
    for p in (OUT_MUT, OUT_CASES, OUT_POS, OUT_UNI, OUT_BAIL):
        m = H.md5_file(p); lines.append(f"{os.path.basename(p)}\t{m}"); log(f"md5 {os.path.basename(p)} {m}")
    with open(GATE, "w", encoding="utf-8", newline="") as f: f.write("\n".join(lines) + "\n")
    log("gate file written:", GATE)

if __name__ == "__main__":
    H.run_section("P2_23", 1, "inputs", section1)
