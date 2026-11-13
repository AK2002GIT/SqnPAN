# SqnPAN_04_check_alphafold_archive.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Structure check of the AlphaFold human archive v6 against UniProt reviewed human 2026_03. No mutation or glycosite result of interest is computed.
# AlphaFold v6 is synced to UniProt 2025_03 (Bertoni et al. 2026, Nucleic Acids Res 54(D1):D358-D362, doi 10.1093/nar/gkaf1226).
# Rules (fixed 2026-10-22): 1A a model is used only if its full sequence equals the UniProt 2026_03 sequence, others are excluded and counted;
# 2A most central fragment (see P2_04_alphafold_helpers.choose_fragment); 3A the tar stays packed, an index of member positions is written.
# Usage (from any folder):  python scripts/SqnPAN_04_check_alphafold_archive.py [all|1|2|3|4]
# Logs: results\logs\P2_04_section<k>_<name>.txt (append mode). Each section reads only saved files, so it can be run on its own.
import sys, os, csv, gzip, re, hashlib, tarfile, time, collections, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import SqnPAN_04_alphafold_helpers as H

P2 = H.PAPER2
LOGS = os.path.join(P2, "results", "logs"); TABLES = os.path.join(P2, "results", "tables"); PROC = os.path.join(P2, "data", "processed")
UNIPROT = os.path.join(P2, "data", "raw", "Uniprot", "stream.gz")
RECORDED_TAR_MD5 = "5bc6bf3812832ceb6fe8f7c9d64fd6be"    # computed after download, 2026-10-22
RECORDED_TAR_SIZE = 5177506304
RECORDED_UNIPROT_MD5 = "af9af22f1e6eea80e96751c93c9da2b3"  # recorded at download
MODELS = os.path.join(PROC, "P2_04_alphafold_models.tsv"); SEQS = os.path.join(PROC, "P2_04_alphafold_sequences.tsv")
USE = os.path.join(PROC, "P2_04_alphafold_use.tsv")
NLINK_RE = re.compile(r'CARBOHYD\s+(\d+);\s*/note="N-linked')  # rough count; evidence codes not used here

class Log:
    def __init__(self, k, name):
        for d in (LOGS, TABLES, PROC): os.makedirs(d, exist_ok=True)
        self.path = os.path.join(LOGS, f"P2_04_section{k}_{name}.txt"); self.k = k
        self.f = open(self.path, "a", encoding="utf-8")
        self("=" * 70); self(f"P2_04 section {k} ({name}) started {time.strftime('%Y-%m-%d %H:%M:%S')}")
    def __call__(self, *a):
        s = " ".join(str(x) for x in a); print(s, flush=True); self.f.write(s + "\n"); self.f.flush()
    def close(self):
        self(f"---- END OF SECTION {self.k} ---- {time.strftime('%Y-%m-%d %H:%M:%S')}"); self.f.close()

def md5_file(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""): h.update(b)
    return h.hexdigest()

def run(k, name, fn):
    log = Log(k, name); t0 = time.time()
    try:
        fn(log)
    except Exception:
        log("ERROR:"); log(traceback.format_exc())
    finally:
        log(f"seconds: {round(time.time() - t0)}"); log.close()

# ---------------- Section 1: archive file and index ----------------
def section1(log):
    size = os.path.getsize(H.TAR_PATH)
    log("tar:", H.TAR_PATH); log("size:", size, "| recorded:", RECORDED_TAR_SIZE, "| same:", size == RECORDED_TAR_SIZE)
    m = md5_file(H.TAR_PATH); log("md5:", m, "| recorded:", RECORDED_TAR_MD5, "| same:", m == RECORDED_TAR_MD5)
    rows, bad, kinds, vers, dates = [], [], collections.Counter(), collections.Counter(), collections.Counter()
    with tarfile.open(H.TAR_PATH, "r:") as tf:
        for mem in tf:
            if not mem.isfile(): continue
            p = H.parse_name(mem.name)
            if p is None: bad.append(mem.name); continue
            acc, frag, ver, kind = p
            kinds[kind] += 1; vers[ver] += 1; dates[time.strftime("%Y-%m-%d", time.gmtime(mem.mtime))] += 1
            rows.append((acc, frag, kind, mem.offset_data, mem.size, mem.name))
    log("member files:", len(rows) + len(bad), "| names not fitting AF-<acc>-F<n>-model_v<k>.(cif|pdb).gz:", len(bad), bad[:5])
    log("by kind:", dict(kinds), "| by version:", dict(vers), "| file dates (UTC):", dict(dates))
    keys = collections.Counter((r[0], r[1], r[2]) for r in rows)
    log("duplicate (accession, fragment, kind):", sum(1 for v in keys.values() if v > 1))
    cif = {(a, f) for (a, f, k) in keys if k == "cif"}; pdb = {(a, f) for (a, f, k) in keys if k == "pdb"}
    log("models with cif:", len(cif), "| with pdb:", len(pdb), "| cif without pdb:", len(cif - pdb), "| pdb without cif:", len(pdb - cif))
    accs = {a for a, f in cif}; frag = collections.Counter(a for a, f in cif)
    log("accessions:", len(accs), "| isoform accessions (with -n):", sum(1 for a in accs if "-" in a),
        "| with several fragments:", sum(1 for v in frag.values() if v > 1), "| extra fragments:", sum(v - 1 for v in frag.values()),
        "| most fragments:", max(frag.items(), key=lambda x: x[1]))
    gaps = [a for a in accs if sorted(f for b, f in cif if b == a) != list(range(1, frag[a] + 1))]
    log("accessions whose fragment numbers are not 1..n:", len(gaps), gaps[:5])
    rows.sort(key=lambda r: r[3])
    with open(H.INDEX_PATH, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["accession", "fragment", "kind", "offset_data", "size", "member_name"]); w.writerows(rows)
    log("index written:", H.INDEX_PATH, "| rows:", len(rows), "| md5:", md5_file(H.INDEX_PATH))

# ---------------- Section 2: read every model through the index ----------------
def section2(log):
    idx = H.load_index(); log("index rows:", len(idx))
    keys = sorted(idx, key=lambda k: idx[k][0])
    c = collections.Counter(); models = []; problems = []; pdb_ca = {}
    with H.open_archive() as fh:
        for n, key in enumerate(keys, 1):
            acc, frag, kind = key
            try:
                text = H.read_member(fh, *idx[key])
            except Exception as e:
                c[f"{kind}_unreadable"] += 1; problems.append((acc, frag, kind, "unreadable " + repr(e))); continue
            c[f"{kind}_read_ok"] += 1
            if kind == "pdb":
                pdb_ca[(acc, frag)] = H.pdb_ca_count(text); continue
            d = H.parse_cif(text); ca = d["ca"]; L = len(d["sequence"])
            checks = {
                "accession_in_file_equals_name": d["accession"] == acc,
                "ca_count_equals_seq_length": len(ca) == L,
                "seq_length_equals_align_span": L == d["align_end"] - d["align_begin"] + 1,
                "label_seq_is_1_to_n": [r[0] for r in ca] == list(range(1, L + 1)),
                "uniprot_numbers_are_begin_to_end": [r[2] for r in ca] == list(range(d["align_begin"], d["align_end"] + 1)),
                "fragment_start_is_1_plus_200_per_step": d["align_begin"] == 1 + 200 * (frag - 1),
                "fragment_length_at_most_1400": L <= 1400,
            }
            for k2, ok in checks.items():
                c[k2 + ("_yes" if ok else "_NO")] += 1
                if not ok: problems.append((acc, frag, "cif", k2))
            models.append((acc, frag, d["gene"], d["align_begin"], d["align_end"], L, len(ca), d["seq_version_date"], d["seq_checksum"], d["sequence"]))
            if n % 5000 == 0: log(f"  read {n} of {len(keys)}")
    for (acc, frag, *_r) in models:
        ok = pdb_ca.get((acc, frag)) == _r[4]
        c["pdb_ca_count_equals_cif_" + ("yes" if ok else "NO")] += 1
        if not ok: problems.append((acc, frag, "pdb", "ca count differs from cif"))
    for k2, v in sorted(c.items()): log(f"{k2}: {v}")
    # stitch fragments into one sequence per accession
    by = collections.defaultdict(list)
    for m in models: by[m[0]].append(m)
    seqs, stitch_bad, lens_last = [], [], collections.Counter()
    for acc, ms in by.items():
        ms.sort(key=lambda m: m[1]); full = {}
        ok = True
        for m in ms:
            for i, aa in enumerate(m[9]):
                p = m[3] + i
                if full.setdefault(p, aa) != aa: ok = False
        n = ms[-1][4]
        if sorted(full) != list(range(1, n + 1)) or not ok: stitch_bad.append(acc)
        if len(ms) > 1:
            if any(m[5] != 1400 for m in ms[:-1]): stitch_bad.append(acc)
            lens_last[ms[-1][5]] += 1
        seqs.append((acc, ms[0][2], len(ms), n, "".join(full[p] for p in range(1, n + 1)) if acc not in stitch_bad else ""))
    log("accessions:", len(seqs), "| fragments that disagree in an overlap, leave a gap, or are not 1400 long before the last:", len(stitch_bad), stitch_bad[:5])
    log("last-fragment lengths in split proteins: min", min(lens_last), "max", max(lens_last), "| split protein lengths: min",
        min(s[3] for s in seqs if s[2] > 1), "max", max(s[3] for s in seqs if s[2] > 1))
    log("residues covered by k fragments (split proteins):", sorted(collections.Counter(
        sum(1 for m in by[a] if m[3] <= p <= m[4]) for a, _g, nf, n, _s in seqs if nf > 1 for p in range(1, n + 1)).items()))
    log("sequence version dates in the models: earliest", min(m[7] for m in models), "latest", max(m[7] for m in models))
    log("problems:", len(problems), problems[:20])
    with open(MODELS, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "fragment", "gene", "align_begin", "align_end", "seq_length", "ca_count", "seq_version_date", "seq_checksum"])
        w.writerows(m[:9] for m in sorted(models))
    with open(SEQS, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["accession", "gene", "n_fragments", "length", "sequence"]); w.writerows(sorted(seqs))
    log("written:", MODELS, "md5", md5_file(MODELS)); log("written:", SEQS, "md5", md5_file(SEQS))

# ---------------- Section 3: comparison with UniProt reviewed 2026_03 (decision 1A) ----------------
def section3(log):
    m = md5_file(UNIPROT); log("UniProt file md5:", m, "| recorded:", RECORDED_UNIPROT_MD5, "| same:", m == RECORDED_UNIPROT_MD5)
    up = {}
    with gzip.open(UNIPROT, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            up[r["Entry"]] = (r["Gene Names (primary)"], r["Sequence"], NLINK_RE.findall(r["Glycosylation"] or ""))
    af = {}
    with open(SEQS, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"): af[r["accession"]] = (int(r["n_fragments"]), r["sequence"])
    log("UniProt reviewed entries:", len(up), "| archive accessions:", len(af))
    both = set(up) & set(af); nomodel = set(up) - set(af); notrev = set(af) - set(up)
    same = {a for a in both if up[a][1] == af[a][1]}; diff = both - same
    difflen = {a for a in diff if len(up[a][1]) != len(af[a][1])}
    log("in both:", len(both), "| identical full sequence (used):", len(same), "| sequence differs (excluded):", len(diff),
        "(length differs", len(difflen), "| same length, residues differ", len(diff - difflen), ")")
    log("reviewed entries without a model:", len(nomodel), "| archive accessions not in the reviewed set (not used):", len(notrev))
    ns = lambda S: sum(len(up[a][2]) for a in S)
    log("annotated N-linked CARBOHYD sites (rough count): all reviewed", ns(up), "| in used proteins", ns(same),
        "| in excluded proteins", ns(diff), "| in proteins without a model", ns(nomodel))
    log("expected from the first pass (2026-10-22): both 20,290; identical 20,249; differ 41 (35 length, 6 same length); no model 141; not reviewed 260; sites 16,244 / 61 / 109")
    with open(USE, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "gene", "status", "uniprot_length", "alphafold_length", "n_fragments", "n_nlinked_sites"])
        for a in sorted(up):
            st = "use" if a in same else ("excluded_sequence_differs" if a in diff else "no_model")
            w.writerow([a, up[a][0], st, len(up[a][1]), len(af[a][1]) if a in af else "", af[a][0] if a in af else "", len(up[a][2])])
    log("written:", USE, "md5", md5_file(USE))
    t1 = os.path.join(TABLES, "P2_04_excluded_sequence_mismatch.tsv")
    with open(t1, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["accession", "gene", "uniprot_length", "alphafold_length", "kind", "n_nlinked_sites"])
        for a in sorted(diff, key=lambda a: up[a][0]):
            w.writerow([a, up[a][0], len(up[a][1]), len(af[a][1]), "length differs" if a in difflen else "same length, residues differ", len(up[a][2])])
    t2 = os.path.join(TABLES, "P2_04_reviewed_without_model.tsv")
    with open(t2, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["accession", "gene", "uniprot_length", "n_nlinked_sites"])
        for a in sorted(nomodel, key=lambda a: up[a][0]): w.writerow([a, up[a][0], len(up[a][1]), len(up[a][2])])
    log("written:", t1); log("written:", t2)
    log("excluded proteins (gene, UniProt length, AlphaFold length, N-linked sites):")
    for a in sorted(diff, key=lambda a: up[a][0]): log(f"  {a}\t{up[a][0]}\t{len(up[a][1])}\t{len(af[a][1])}\t{len(up[a][2])}")
    lens = collections.Counter("<=30 aa" if len(up[a][1]) <= 30 else ">30 aa" for a in nomodel)
    log("reviewed without model by length:", dict(lens), "| largest site counts:",
        sorted(((len(up[a][2]), up[a][0]) for a in nomodel), reverse=True)[:10])

# ---------------- Section 4: fragment rule check (decision 2A) ----------------
def section4(log):
    ranges = collections.defaultdict(list)
    with open(MODELS, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"): ranges[r["accession"]].append((int(r["fragment"]), int(r["align_begin"]), int(r["align_end"])))
    use = {}
    with open(USE, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["status"] == "use": use[r["accession"]] = (r["gene"], int(r["uniprot_length"]))
    split = [a for a in use if len(ranges[a]) > 1]
    log("used proteins:", len(use), "| of them split into fragments:", len(split))
    nres = notcov = 0; mind = collections.Counter(); chosen_count = collections.Counter()
    for a in split:
        n = use[a][1]
        for p in range(1, n + 1):
            f = H.choose_fragment(ranges[a], p, n); nres += 1
            if f is None: notcov += 1; continue
            _, b, e = next(x for x in ranges[a] if x[0] == f)
            d = min(p - b if b > 1 else 10 ** 9, e - p if e < n else 10 ** 9)
            mind["<200" if d < 200 else ("200-599" if d < 600 else ">=600")] += 1
    log("residues in split proteins:", nres, "| not covered by any chosen fragment:", notcov)
    log("distance of each residue to the nearest artificial cut in its chosen fragment:", dict(mind))
    up = {}
    with gzip.open(UNIPROT, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["Entry"] in split: up[r["Entry"]] = [int(x) for x in NLINK_RE.findall(r["Glycosylation"] or "")]
    sites = [(a, p) for a in split for p in up.get(a, [])]
    cover = collections.Counter(sum(1 for _f, b, e in ranges[a] if b <= p <= e) for a, p in sites)
    log("annotated N-linked sites in split used proteins:", len(sites), "| in proteins:", sum(1 for a in split if up.get(a)),
        "| sites by number of covering fragments:", sorted(cover.items()))
    log("tie-break check: position 701 of a 3000-residue protein with fragments F1 1-1400, F2 201-1600 ->",
        H.choose_fragment([(1, 1, 1400), (2, 201, 1600), (3, 401, 1800), (4, 601, 2000), (5, 801, 2200), (6, 1001, 2400), (7, 1201, 2600), (8, 1401, 2800), (9, 1601, 3000)], 701, 3000),
        "(F1 distance to its cut end 699 vs F2 500 ... F4 100 -> F1 expected)")
    with H.open_archive() as fh:
        idx = H.load_index(); d = H.parse_cif(H.read_model(fh, idx, "A0A024R1R8", 1))
        log("read-back test A0A024R1R8 F1 via index:", d["gene"], d["align_begin"], d["align_end"], len(d["sequence"]), "residues")

SECTIONS = {"1": ("index", section1), "2": ("read_models", section2), "3": ("uniprot_compare", section3), "4": ("fragment_rule", section4)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for k in (["1", "2", "3", "4"] if which == "all" else [which]):
        run(k, *SECTIONS[k])
# ---- END OF SCRIPT ----
