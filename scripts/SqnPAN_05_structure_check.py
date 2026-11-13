# SqnPAN_05_structure_check.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Structure check for SqnPAN (the analysis plan, section 12): inputs, md5, shapes, and counts for analysis plan rules 2A to 5A and the decision 10 table.
# No mutation-at-site result is computed: the MAF section counts rows, file choice and identifier agreement only (no counts per gene, site or tier).
# There is no random step in this script (nothing to seed); a second run must give the same logs apart from the time lines.
# Sections (each reads only saved files, so each can be run alone):
#   1 inputs           md5 and size guards, MAF file sizes against the manifest, case-project table, TCGA-CDR and Bailey table shapes
#   2 uniprot_tables   the second UniProt table against stream.gz (decision 10A: sequence identity), column and feature-text census
#   3 sites_tiers      N-linked records, evidence codes, sequon status, atlas rule 3A, tiers T1/T2/T3 (rule 2A), strata sizes, decision 10 table
#   4 maf_files        one pass over all MAF files: file/case/barcode checks, file choice per case (rule 4A), transcript rule counts (rule 5A)
# Usage:  python scripts/SqnPAN_05_structure_check.py [all|1|2|3|4]
# Logs: results\logs\P2_05_section<k>_<name>.txt (append mode). Section 4 also writes data\processed\P2_05_maf_file_summary.tsv and
# results\tables\P2_05_case_file_selection.tsv. Test options (quick checks only): limit=N (first N MAF files), P2_LOGDIR / P2_TABLEDIR / P2_PROCDIR.
import sys, os, re, csv, gzip, time, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import SqnPAN_05_helpers as H

# ---- recorded values (guards; a difference is printed, never hidden) ----
REC_MD5 = {
    "UniProt stream.gz": (H.UNIPROT_MAIN, "af9af22f1e6eea80e96751c93c9da2b3"),
    "UniProt features table": (H.UNIPROT_FEAT, "bf26668a884779609e7c1fe46bb6aec4"),
    "N-GlycositeAtlas HumanAll.xlsx": (H.ATLAS, "23d97d828400e59f1ca0837b5d59cf8a"),
    "TCGA-CDR table": (H.CDR, "a4591b2dcee39591f59e5e25a6ce75fa"),
    "Bailey 2018 Table S1 file": (H.BAILEY, "144a9ccbbf038236decc1226df24620b"),
    "case-project table": (H.CASE_PROJECT, "b6bd436e630584d4738d9bdacf4d5b7a"),
    "AlphaFold index": (os.path.join(H.PROC, "P2_04_alphafold_index.tsv"), "b674fa1d5c67da8bbfa5145f30ebadf1"),
    "AlphaFold use table": (H.USE_TSV, "bc009789438ac60310524d4ef3e6f901"),
}
REC_SIZE = {"UniProt features table": 9497308, "TCGA-CDR table": 2945129, "Bailey 2018 Table S1 file": 247939}
REC_TAR = (os.path.join(H.RAW, "alphafold", "UP000005640_9606_HUMAN_v6.tar"), 5177506304)
MAIN_COLS = ['Entry', 'Entry Name', 'Gene Names (primary)', 'Length', 'Sequence', 'Glycosylation', 'Ensembl', 'AlphaFoldDB']
FEAT_COLS = ['Entry', 'Length', 'Sequence', 'Topological domain', 'Transmembrane', 'Intramembrane', 'Signal peptide', 'Subcellular location [CC]',
             'Domain [FT]', 'Disulfide bond', 'Active site', 'Binding site']
ATLAS_COLS = ['Database', 'Accession (UniProtKB)', 'Gene names', 'Protein names', 'Glycosylation location', 'Glycosite-containing peptides',
              'Glycosite_20AAs', 'Identification resources', 'Years', 'References']
REC_MAF = {"files": 10640, "rows": 2570542, "empty_files": 91, "columns": 140, "cases": 10190, "projects": 33, "cases_multi": 385, "max_files_per_case": 6}
LIMIT = None; SELECTED = None
for a in sys.argv[2:]:
    if a.startswith("limit="): LIMIT = int(a[6:])

def same(a, b): return "same" if a == b else f"DIFFERENT (recorded {b})"
def pct(a, b): return f"{100.0 * a / b:.1f}%" if b else "n/a"

# ======================= Section 1: inputs =======================
def section1(log):
    log("Recorded values are from the first download and the analysis plan. A line saying DIFFERENT is a finding, not an error in the script.")
    for label, (path, md5) in REC_MD5.items():
        if not os.path.exists(path): log(f"{label}: FILE NOT FOUND {path}"); continue
        m = H.md5_file(path); sz = os.path.getsize(path)
        extra = f" | size {sz} {same(sz, REC_SIZE[label])}" if label in REC_SIZE else f" | size {sz}"
        log(f"{label}: md5 {m} {same(m, md5)}{extra}")
    p, sz = REC_TAR
    log(f"AlphaFold tar: size {os.path.getsize(p) if os.path.exists(p) else 'NOT FOUND'} {same(os.path.getsize(p) if os.path.exists(p) else None, sz)} (md5 of the tar was checked by SqnPAN_04_check_alphafold_archive.py and is not repeated)")
    fx = [x for x in os.listdir(os.path.join(H.PAPER2, "foldx5_Windows")) if x.lower().startswith("foldx") and x.lower().endswith(".exe")] if os.path.isdir(os.path.join(H.PAPER2, "foldx5_Windows")) else []
    for x in fx:
        fp = os.path.join(H.PAPER2, "foldx5_Windows", x); log(f"FoldX file {x}: size {os.path.getsize(fp)} | md5 {H.md5_file(fp)} (first md5 record, none was recorded before)")
    if not fx: log("FoldX exe: NOT FOUND in the foldx5_Windows folder")

    log("--- MAF files against the manifest (size only; md5 of every file was checked in SqnPAN_01/SqnPAN_02) ---")
    man = list(csv.DictReader(open(H.MANIFEST, newline=""), delimiter="\t"))
    c = collections.Counter(); folders = set(os.listdir(H.MAF_DIR))
    for r in man:
        fp = os.path.join(H.MAF_DIR, r["id"], r["filename"])
        if not os.path.exists(fp): c["missing"] += 1
        elif os.path.getsize(fp) != int(r["size"]): c["size_differs"] += 1
        else: c["ok"] += 1
    log(f"manifest rows: {len(man)} {same(len(man), REC_MAF['files'])} | states: {dict(collections.Counter(r['state'] for r in man))} | {dict(c)}")
    log(f"folders in maf_tcga: {len(folders)} | not in manifest: {len(folders - {r['id'] for r in man})} | in manifest without folder: {len({r['id'] for r in man} - folders)}")

    log("--- case-project table ---")
    cp = list(csv.DictReader(open(H.CASE_PROJECT, newline=""), delimiter="\t"))
    log("columns:", list(cp[0].keys()))
    cases = collections.defaultdict(list)
    for r in cp: cases[r["cases.0.case_id"]].append(r["file_id"])
    mids = {r["id"] for r in man}
    log(f"rows: {len(cp)} | unique file_id: {len({r['file_id'] for r in cp})} | files in manifest and in table: {len(mids & {r['file_id'] for r in cp})} of {len(mids)}")
    log(f"cases: {len(cases)} {same(len(cases), REC_MAF['cases'])} | projects: {len({r['cases.0.project.project_id'] for r in cp})} {same(len({r['cases.0.project.project_id'] for r in cp}), REC_MAF['projects'])}")
    fpc = collections.Counter(len(v) for v in cases.values())
    log(f"files per case: {dict(sorted(fpc.items()))} | cases with more than one file: {sum(v for k, v in fpc.items() if k > 1)} {same(sum(v for k, v in fpc.items() if k > 1), REC_MAF['cases_multi'])} | max {max(fpc)} {same(max(fpc), REC_MAF['max_files_per_case'])}")
    log("cases per project:", dict(sorted(collections.Counter(r["cases.0.project.project_id"] for r in {x["cases.0.case_id"]: x for x in cp}.values()).items())))
    log(f"submitter_id format TCGA-XX-XXXX for all rows: {all(re.match(r'^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$', r['cases.0.submitter_id']) for r in cp)}")

    log("--- TCGA-CDR table ---")
    names = H.xlsx_sheet_names(H.CDR); log("sheets:", names)
    rows = H.read_xlsx(H.CDR, "TCGA-CDR"); hdr = rows[0]; body = [r + [""] * (len(hdr) - len(r)) for r in rows[1:] if any(x != "" for x in r)]
    log(f"header row = sheet row 1 | columns ({len(hdr)}): {hdr}")
    ix = {h: i for i, h in enumerate(hdr)}
    log(f"patients: {len(body)} | unique bcr_patient_barcode: {len({r[ix['bcr_patient_barcode']] for r in body})} | cancer types ({len({r[ix['type']] for r in body})}):",
        dict(sorted(collections.Counter(r[ix['type']] for r in body).items())))
    na = {"", "#N/A", "[Not Available]", "[Not Applicable]", "[Unknown]", "[Discrepancy]", "[Completed]"}
    for col in ["age_at_initial_pathologic_diagnosis", "gender", "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time", "Redaction"]:
        vals = collections.Counter(r[ix[col]] for r in body)
        log(f"  {col}: missing-type values {sum(v for k, v in vals.items() if k in na)} of {len(body)} | distinct values {len(vals)} | most common {vals.most_common(4)}")
    cdr_ids = {r[ix['bcr_patient_barcode']] for r in body}; maf_ids = {r["cases.0.submitter_id"] for r in {x["cases.0.case_id"]: x for x in cp}.values()}
    log(f"MAF cases in TCGA-CDR: {len(maf_ids & cdr_ids)} of {len(maf_ids)} | CDR patients without MAF case: {len(cdr_ids - maf_ids)}")
    t = {r[ix['bcr_patient_barcode']]: r[ix['type']] for r in body}; bad = collections.Counter()
    for r in {x["cases.0.case_id"]: x for x in cp}.values():
        s = r["cases.0.submitter_id"]
        if s in t and t[s] != r["cases.0.project.project_id"].replace("TCGA-", ""): bad[(r["cases.0.project.project_id"], t[s])] += 1
    log(f"cases whose CDR type differs from the project name: {sum(bad.values())} {dict(bad)}")
    log("TCGA-CDR_Notes sheet: endpoint advice is general (PFI preferred over OS, DSS accurate only for CESC, PAAD, UVM, DFI absent for SKCM, THYM, UVM) and refers to Table 3 of Liu et al. 2018 for the per-type choice: the [confirm] item in analysis plan section 9 is NOT settled by the Notes sheet.")

    log("--- Bailey et al. 2018 Table S1 ---")
    rows = H.read_xlsx(H.BAILEY, "Table S1"); hi = next(i for i, r in enumerate(rows) if r and r[0] == "Gene"); hdr = rows[hi]
    body = [r + [""] * (len(hdr) - len(r)) for r in rows[hi + 1:] if r and r[0] != ""]
    log(f"header at sheet row {hi + 1}: {[h for h in hdr if h]}")
    log(f"data rows with a gene: {len(body)} | rows after the header without a gene: {len(rows) - hi - 1 - len(body)} | unique genes: {len({r[0] for r in body})}")
    for j, h in enumerate(hdr):
        if h: log(f"  column {j + 1} {h!r}: non-empty {sum(1 for r in body if r[j] not in ('', '#N/A', 'NA'))} of {len(body)}")
    di = hdr.index("Decision"); ci = hdr.index("Cancer")
    log("Decision values:", dict(collections.Counter(r[di] for r in body)), "| unique genes per Decision:", {k: len({r[0] for r in body if r[di] == k}) for k in {r[di] for r in body}})
    log("Cancer values:", dict(sorted(collections.Counter(r[ci] for r in body).items())))
    U = H.load_uniprot_main(); genes = set(); alias = set()
    for u in U.values(): genes.add(u["gene"].split(" ")[0] if u["gene"] else "")
    bg = {r[0] for r in body}
    log(f"Bailey genes equal to a UniProt primary gene name: {len(bg & genes)} of {len(bg)} | not found: {sorted(bg - genes)[:40]}")

# ======================= Section 2: UniProt tables =======================
def section2(log):
    for label, path, cols in (("stream.gz", H.UNIPROT_MAIN, MAIN_COLS), ("features table", H.UNIPROT_FEAT, FEAT_COLS)):
        with gzip.open(path, "rt", encoding="utf-8") as f: hdr = f.readline().rstrip("\n").split("\t")
        log(f"{label}: columns {hdr} | as expected: {hdr == cols}")
    U = H.load_uniprot_main(); F = H.load_uniprot_features()
    log(f"entries: stream.gz {len(U)} | features table {len(F)} | same entry set: {set(U) == set(F)} | only in one: {len(set(U) ^ set(F))}")
    seq_same = sum(1 for k in U if k in F and U[k]["seq"] == F[k]["Sequence"]); len_same = sum(1 for k in U if k in F and U[k]["length"] == int(F[k]["Length"]))
    log(f"sequence identical between the two tables: {seq_same} of {len(U)} | Length equal: {len_same} | Length equals len(Sequence) in stream.gz: {sum(1 for u in U.values() if u['length'] == len(u['seq']))} | in features table: {sum(1 for r in F.values() if int(r['Length']) == len(r['Sequence']))}")
    log("Decision 10A check (same release): " + ("PASSED, sequences identical in all entries" if seq_same == len(U) == len(F) else "NOT PASSED, see counts above; stop and decide"))
    alpha = collections.Counter(ch for u in U.values() for ch in u["seq"])
    log("residue letters in all sequences:", dict(sorted(alpha.items(), key=lambda x: -x[1])))
    keys = {"Topological domain": "TOPO_DOM", "Transmembrane": "TRANSMEM", "Intramembrane": "INTRAMEM", "Signal peptide": "SIGNAL", "Domain [FT]": "DOMAIN",
            "Disulfide bond": "DISULFID", "Active site": "ACT_SITE", "Binding site": "BINDING"}
    for col, key in keys.items():
        nonempty = 0; nfeat = 0; keyc = collections.Counter(); kinds = collections.Counter(); out_of_range = 0; unp = 0; notes = collections.Counter(); raw_keys = collections.Counter()
        for acc, r in F.items():
            t = r[col]
            if not t.strip(): continue
            nonempty += 1; L = int(r["Length"])
            raw_keys.update(re.findall(r"(?:^|(?<=; ))([A-Z_]{3,10})(?=\s+[<>?\d])", t))
            for ft in H.parse_features(t):
                nfeat += 1; keyc[ft["key"]] += 1
                kinds["uncertain(<,>,?)" if ft["uncertain"] else ("single" if ft["start"] == ft["end"] else "range")] += 1
                if ft["start"] is None or ft["end"] is None: unp += 1
                elif ft["start"] > L or ft["end"] > L or ft["start"] > ft["end"]: out_of_range += 1
                if key in ("TOPO_DOM", "TRANSMEM", "INTRAMEM"): notes[re.sub(r"\d+", "#", ft["note"])] += 1
        log(f"{col}: entries with text {nonempty} | features parsed {nfeat} | keys {dict(keyc)} (expected {key}) | raw key tokens {dict(raw_keys)} | position forms {dict(kinds)} | without a number {unp} | outside 1..Length or reversed {out_of_range}")
        if notes: log(f"   notes (top 14): {notes.most_common(14)}")
    n = sum(1 for r in F.values() if r["Subcellular location [CC]"].strip()); n2 = sum(1 for r in F.values() if r["Subcellular location [CC]"].startswith("SUBCELLULAR LOCATION:"))
    log(f"Subcellular location [CC]: entries with text {n} | starting with 'SUBCELLULAR LOCATION:' {n2} (free text, not parsed in this script)")
    log("entries with a signal peptide:", sum(1 for r in F.values() if r["Signal peptide"].strip()), "| with topological domain:", sum(1 for r in F.values() if r["Topological domain"].strip()),
        "| with transmembrane:", sum(1 for r in F.values() if r["Transmembrane"].strip()))

# ======================= Section 3: sites, tiers, strata, decision 10 table =======================
def section3(log):
    U = H.load_uniprot_main(); F = H.load_uniprot_features()
    use = {r["accession"]: r for r in csv.DictReader(open(H.USE_TSV, newline=""), delimiter="\t")}
    T1CODES = {"ECO:0000269", "ECO:0007744"}
    note_cls = collections.Counter(); kinds = collections.Counter(); rough = 0; ev_rough = collections.Counter(); ev_asn = collections.Counter()
    sites = {}; dup = 0
    for acc, u in U.items():
        for r in H.parse_glyco(u["glyco"]):
            ncls = H.note_class(r["note"]); note_cls[ncls] += 1
            kind, pos = H.pos_kind(r["pos_raw"]); kinds[(kind, "N-linked" if r["note"].startswith("N-linked") else "not N-linked")] += 1
            if not r["note"].startswith("N-linked") or kind != "single": continue
            combo = "+".join(c.replace("ECO:", "") for c in r["codes"]) or "no evidence code"
            rough += 1; ev_rough[combo] += 1
            if ncls == "N-linked asparagine": ev_asn[combo] += 1
            if (acc, pos) in sites: dup += 1; sites[(acc, pos)]["codes"] |= set(r["codes"]); continue
            sites[(acc, pos)] = {"codes": set(r["codes"]), "ncls": ncls, "status": H.sequon_status(u["seq"], pos)}
    log(f"CARBOHYD records by note class: {dict(note_cls)}")
    log(f"records: N-linked/not by position form: {dict(kinds)}")
    log(f"N-linked records at a single canonical position (the count used in SqnPAN_04 and the analysis plan): {rough} {same(rough, 16244)} | unique (entry, position): {len(sites)} | repeated positions merged: {dup}")
    log("evidence-code combinations, all N-linked single-position records:", dict(ev_rough.most_common()))
    log("evidence-code combinations, N-linked asparagine records only:", dict(ev_asn.most_common()))
    log("analysis plan rule 2A counts for comparison: 12718 ECO:0000255 only, 2384 ECO:0000269 only, 638 ECO:0000269 + ECO:0007744, 131 ECO:0000305 only, 91 without code; ~3200 with 269 or 7744")
    n_t1 = sum(1 for s in sites.values() if s["codes"] & T1CODES); log(f"records with ECO:0000269 or ECO:0007744 (any other codes allowed): {n_t1} (analysis plan: about 3200)")
    st = collections.Counter(); bytier = collections.Counter()
    for k, s in sites.items():
        tier = "T1-code" if s["codes"] & T1CODES else "other-code"; st[(tier, s["status"])] += 1; bytier[(s["ncls"], s["status"])] += 1
    log("sequon status of UniProt N-linked sites by evidence class:", dict(sorted(st.items())))
    log("sequon status by note class:", dict(sorted(bytier.items())))

    # ---- atlas (rule 3A) ----
    rows = H.read_xlsx(H.ATLAS, "Glycosite and glycoprotein DB"); hi = next(i for i, r in enumerate(rows) if r and r[0] == "Database"); hdr = rows[hi]
    log(f"atlas header at sheet row {hi + 1}, as expected: {hdr[:10] == ATLAS_COLS}")
    body = [r + [""] * (len(hdr) - len(r)) for r in rows[hi + 1:] if r and r[0] != ""]; ix = {h: i for i, h in enumerate(hdr)}
    dbv = collections.Counter(r[ix["Database"]] for r in body)
    log(f"atlas data rows: {len(body)} | Database values: {dict(dbv)} | unique (accession, location) all: {len({(r[ix['Accession (UniProtKB)']], r[ix['Glycosylation location']]) for r in body})}")
    rev = [r for r in body if r[ix["Database"]] == "Reviewed"]
    uniq = {}; badloc = 0; notin = set(); win = collections.Counter()
    for r in rev:
        acc = r[ix["Accession (UniProtKB)"]]; loc = r[ix["Glycosylation location"]]
        try: pos = int(float(loc))
        except ValueError: badloc += 1; continue
        if acc not in U: notin.add(acc); continue
        uniq.setdefault((acc, pos), []).append(r[ix["Glycosite_20AAs"]])
    log(f"Reviewed rows: {len(rev)} | location not a number: {badloc} | accessions not in UniProt 2026_03 reviewed: {len(notin)} | unique (accession, position) in UniProt 2026_03: {len(uniq)} | proteins: {len({a for a, p in uniq})}")
    ast = collections.Counter(H.sequon_status(U[a]["seq"], p) for (a, p) in uniq)
    nn = ast["residue_not_N"] + ast["position_outside_sequence"]
    log(f"atlas sites by sequon status in 2026_03 sequence: {dict(ast)}")
    log(f"analysis plan rule 3A: 12961 reviewed sites, 429 with no N at the stated position | now: {len(uniq)} sites {same(len(uniq), 12961)}, {nn} with no N (residue not N or position beyond the sequence) {same(nn, 429)}")
    # window check: is the 20-residue window centred on the stated position (a test of the position convention)
    for (a, p), ws in uniq.items():
        w = ws[0]; s = U[a]["seq"]
        if len(w) == 41 and p - 21 >= 0 and p + 20 <= len(s): win["41 letters; centre = seq[pos-1]" if w[20].upper() == s[p - 1] and w.upper() == s[p - 21:p + 20] else "41 letters; does not match"] += 1
        else: win[f"length {len(w)} or near a terminus"] += 1
    log("atlas Glycosite_20AAs window against the 2026_03 sequence:", dict(win))
    atlas_valid = {k for k in uniq if H.sequon_status(U[k[0]]["seq"], k[1]) == "N-X-S/T"}
    atlas_dropped = len(uniq) - len(atlas_valid)

    # ---- tiers (rule 2A) ----
    uni_ann = set(sites)
    T1u = {k for k, s in sites.items() if s["status"] == "N-X-S/T" and s["codes"] & T1CODES}
    T2u = {k for k, s in sites.items() if s["status"] == "N-X-S/T" and not (s["codes"] & T1CODES)}
    atypical = {k for k, s in sites.items() if s["status"] != "N-X-S/T"}
    allseq = set(); xc = 0; nps = 0
    for acc, u in U.items():
        a, b, c = H.sequon_positions(u["seq"]); allseq |= {(acc, i) for i in a}; xc += len(b); nps += len(c)
    T1 = T1u | atlas_valid; T2 = T2u - T1; T3 = allseq - T1 - T2
    log(f"UniProt-annotated N-linked sites that are N-X-S/T: T1 by code {len(T1u)}, other codes {len(T2u)} | annotated but not N-X-S/T (counted, not analysed): {len(atypical)}")
    log(f"atlas sites passing rule 3A (N-X-S/T at the stated position): {len(atlas_valid)} | dropped: {atlas_dropped} of {len(uniq)}")
    log(f"atlas sites among UniProt T1-code sites: {len(atlas_valid & T1u)} | atlas sites at UniProt other-code sites (T2 to T1): {len(atlas_valid & T2u)} | atlas sites with no UniProt N-linked annotation: {len(atlas_valid - uni_ann)} | atlas sites at UniProt-annotated non-sequon positions: {len(atlas_valid & atypical)}")
    log(f"all N-X-S/T sequons in the 2026_03 sequences: {len(allseq)} | N-X-C (not counted): {xc} | N-P-S/T (not a sequon): {nps}")
    log(f"TIERS: T1 {len(T1)} sites in {len({a for a, p in T1})} proteins | T2 {len(T2)} sites in {len({a for a, p in T2})} proteins | T3 {len(T3)} sites in {len({a for a, p in T3})} proteins | sum {len(T1) + len(T2) + len(T3)} equals all sequons: {len(T1) + len(T2) + len(T3) == len(allseq)} | tiers disjoint: {not (T1 & T2 or T1 & T3 or T2 & T3)}")
    # the +3 proline flag and the position-3 S/T
    pro3 = collections.Counter()
    tier_of = {}
    for k in T1: tier_of[k] = "T1"
    for k in T2: tier_of[k] = "T2"
    for k in T3: tier_of[k] = "T3"
    for (a, p), t in tier_of.items():
        s = U[a]["seq"]
        pro3[(t, "P at +3 (N-X-S/T-P)" if p + 3 <= len(s) and s[p + 2] == "P" else "no P at +3")] += 1
    log("sequons by tier and proline at position +3 (flag only):", dict(sorted(pro3.items())))

    # ---- strata and denominators for Q2, split by model availability ----
    prot_sites = collections.defaultdict(lambda: collections.defaultdict(set))   # acc -> tier -> N positions
    for (a, p), t in tier_of.items(): prot_sites[a][t].add(p)
    agg = collections.Counter()
    for acc, u in U.items():
        s = u["seq"]; ps = prot_sites.get(acc, {})
        seqonN = {p for t in ps for p in ps[t]}; sequonST = {p + 2 for p in seqonN}
        nN = s.count("N"); nS = s.count("S"); nT = s.count("T")
        bgN = nN - len(seqonN); bgST = (nS + nT) - len({q for q in sequonST if q <= len(s)})
        agg["proteins"] += 1; agg["N residues"] += nN; agg["S/T residues"] += nS + nT
        for t in ("T1", "T2", "T3"):
            n = len(ps.get(t, ())); agg[f"{t} sequon N"] += n
            if n: agg[f"proteins with {t}"] += 1
            if n and bgN > 0: agg[f"{t}: proteins with sequon N and a non-sequon N (informative N stratum)"] += 1
            if n and bgST > 0: agg[f"{t}: proteins with sequon S/T and a non-sequon S/T (informative S/T stratum)"] += 1
        agg["non-sequon N"] += bgN; agg["non-sequon S/T"] += bgST
    log("Q2 denominators (position classes of analysis plan section 6, all 20431 proteins):", dict(agg))
    ST_roles = collections.Counter()
    for (a, p) in tier_of:
        s = U[a]["seq"]
        ST_roles["position-3 residue is S" if s[p + 1] == "S" else "position-3 residue is T"] += 1
        if s[p] in "ST": ST_roles["X (position 2) is itself S/T"] += 1
        if s[p] == "N": ST_roles["X (position 2) is N"] += 1
    log("sequon composition:", dict(ST_roles))
    sub = collections.Counter()
    for (a, p), t in tier_of.items():
        stt = use.get(a, {}).get("status", "not_in_use_table"); frag = use.get(a, {}).get("n_fragments", "0")
        sub[(t, stt, "split into fragments" if frag not in ("0", "1", "") else "single model")] += 1
    log("sequons by tier, AlphaFold status (SqnPAN_04 use table) and fragment model:", dict(sorted(sub.items())))
    roughs = collections.Counter()
    for (a, p) in sites: roughs[use.get(a, {}).get("status", "not_in_use_table")] += 1
    log("UniProt N-linked single-position sites by AlphaFold status:", dict(roughs), "| SqnPAN_04 reported 16074 used, 61 excluded, 109 without model (counted with the rough regex)")

    # ---- decision 10 table: topology around every N ----
    topo = {}; notes = collections.Counter()
    for acc, r in F.items():
        topo[acc] = {"topo": [(f["start"], f["end"], f["note"]) for f in H.parse_features(r["Topological domain"]) if f["start"] is not None and f["end"] is not None],
                     "tm": [(f["start"], f["end"]) for f in H.parse_features(r["Transmembrane"]) if f["start"] is not None and f["end"] is not None],
                     "im": [(f["start"], f["end"]) for f in H.parse_features(r["Intramembrane"]) if f["start"] is not None and f["end"] is not None],
                     "sig": [(f["start"], f["end"]) for f in H.parse_features(r["Signal peptide"]) if f["start"] is not None and f["end"] is not None]}
    def classify(acc, p):
        d = topo[acc]; cyto = noncyto = False
        for a, b, n in d["topo"]:
            if a <= p <= b:
                if n.lower().startswith("cytoplasmic"): cyto = True
                else: noncyto = True
        if noncyto: cat = "TOPO_DOM non-cytoplasmic"
        elif cyto: cat = "TOPO_DOM cytoplasmic"
        elif any(a <= p <= b for a, b in d["tm"]): cat = "TRANSMEM"
        elif any(a <= p <= b for a, b in d["im"]): cat = "INTRAMEM"
        elif any(a <= p <= b for a, b in d["sig"]): cat = "SIGNAL"
        elif d["topo"] or d["tm"]: cat = "no region at this position (protein has topology)"
        else: cat = "protein without topology annotation"
        return cat, cyto, noncyto, bool(d["sig"])
    hasT1 = {a for a, p in T1}
    cnt = collections.defaultdict(collections.Counter)
    for acc, u in U.items():
        s = u["seq"]
        for m in re.finditer("N", s):
            p = m.start() + 1; cat, cyto, noncyto, sig = classify(acc, p)
            grp = ["all N"] + ([tier_of[(acc, p)] + " sequon N"] if (acc, p) in tier_of else [])
            readA = (acc in hasT1 and noncyto) or (sig and not cyto)      # reading A of the brief sentence: (T1 AND non-cytoplasmic) OR (signal AND not cytoplasmic)
            readB = acc in hasT1 and (noncyto or (sig and not cyto))      # reading B: T1 AND (non-cytoplasmic OR (signal AND not cytoplasmic))
            for g in grp:
                cnt[g][cat] += 1; cnt[g]["plausible, reading A"] += readA; cnt[g]["plausible, reading B"] += readB; cnt[g]["_total"] += 1
    for g, c in cnt.items(): log(f"decision 10 table, {g}: {dict(c)}")
    log("Analysis plan section 5 says: 'plausible if the protein carries at least one T1 site AND the new N lies in a non-cytoplasmic topological domain, or the protein has a signal peptide and no cytoplasmic annotation at that position'.")
    log("The sentence can be read two ways (A: the T1 condition covers only the first clause; B: it covers both). Both are counted above. The reading is fixed in the pre-specification of SqnPAN_06; the plan text is not changed here.")
    log("Definition used for 'non-cytoplasmic' in this count (from the data, not from memory): any TOPO_DOM whose note does not start with 'Cytoplasmic'. The notes seen are listed in section 2.")

# ======================= Section 4: MAF files =======================
BC = re.compile(r"^(TCGA)-(\w{2})-(\w{4})-(\d{2})([A-Z])-(\d{2})([A-Z])-")
def section4(log):
    U = H.load_uniprot_main()
    revmap = collections.defaultdict(list); forms = collections.Counter()
    for acc, u in U.items():
        for m in re.finditer(r"(ENST\d+)(?:\.\d+)?\s*(?:\[([^\]]*)\])?", u["ensembl"]):
            br = m.group(2)
            if br is None: tag = None; forms["no bracket"] += 1
            else:
                mm = re.match(r"^(.+?)(?:-(\d+))?$", br); tag = mm.group(2)
                forms["bracket without isoform number" if tag is None else ("bracket -1" if tag == "1" else "bracket other isoform")] += 1
                if not br.startswith(acc): forms["bracket accession differs from entry"] += 1
            revmap[m.group(1)].append((acc, tag in (None, "1")))
    log(f"UniProt Ensembl field: transcripts {len(revmap)} | forms {dict(forms)} | transcripts listed under more than one entry: {sum(1 for v in revmap.values() if len({a for a, o in v}) > 1)}")
    man = list(csv.DictReader(open(H.MANIFEST, newline=""), delimiter="\t"))
    if LIMIT: man = man[:LIMIT]; log(f"TEST RUN: only the first {LIMIT} manifest files are read")
    cp = {r["file_id"]: r for r in csv.DictReader(open(H.CASE_PROJECT, newline=""), delimiter="\t")}
    need = ["Variant_Classification", "SWISSPROT", "Protein_position", "Amino_acids", "CANONICAL", "all_effects", "Tumor_Sample_Barcode", "case_id", "HGVSp_Short"]
    hdr_variants = collections.Counter(); vc_total = collections.Counter(); meta = []; t0 = time.time()
    reasons = collections.defaultdict(collections.Counter)   # file id -> counter of pass / fail / rescue categories (by canonical flag)
    for n, r in enumerate(man, 1):
        fid = r["id"]; d = os.path.join(H.MAF_DIR, fid); path = os.path.join(d, r["filename"]); c = reasons[fid]
        hdr = None; bcs = collections.Counter(); case_mismatch = 0; rows = 0
        with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
            for line in f:
                if line[0] == "#": continue
                if hdr is None:
                    hdr = line.rstrip("\n").split("\t"); hdr_variants[tuple(hdr)] += 1; ix = {k: hdr.index(k) for k in need}; continue
                v = line.rstrip("\n").split("\t"); rows += 1
                vc = v[ix["Variant_Classification"]]; vc_total[vc] += 1; bcs[v[ix["Tumor_Sample_Barcode"]]] += 1
                if v[ix["case_id"]] != cp[fid]["cases.0.case_id"]: case_mismatch += 1
                if vc != "Missense_Mutation": continue
                canon = "canonical" if v[ix["CANONICAL"]] == "YES" else "non-canonical"; c["missense"] += 1; c[canon + " missense"] += 1
                acc = v[ix["SWISSPROT"]].split(".")[0].split("&")[0]
                if "&" in v[ix["SWISSPROT"]]: c["SWISSPROT has several accessions (&)"] += 1
                reason = None
                if not acc: reason = "no SWISSPROT accession"
                elif acc not in U: reason = "accession not in UniProt reviewed 2026_03"
                else:
                    try:
                        pos_s, tot_s = v[ix["Protein_position"]].split("/"); pos = int(pos_s); tot = int(tot_s); ref, alt = v[ix["Amino_acids"]].split("/")
                    except Exception: pos = None; reason = "position or amino acids unparsable"
                    if reason is None:
                        s = U[acc]["seq"]; lenok = tot == len(s); resok = 1 <= pos <= len(s) and len(ref) == 1 and s[pos - 1] == ref
                        reason = None if (lenok and resok) else ("length and residue both differ" if not lenok and not resok else ("length differs" if not lenok else "reference residue differs"))
                if reason is None: c["pass"] += 1; c[canon + " pass"] += 1; continue
                c["fail"] += 1; c[canon + " fail"] += 1; c["fail: " + reason] += 1
                # rescue through all_effects (rule 5A): transcript listed by UniProt for an entry without isoform tag or with -1, reference residue matches
                ok_t = set(); other_t = set(); n_eff = 0
                for e in v[ix["all_effects"]].split(";"):
                    fl = e.split(",")
                    if len(fl) < 4 or "missense_variant" not in fl[1]: continue
                    m = re.match(r"^p\.([A-Z])(\d+)([A-Z])$", fl[2])
                    if not m: continue
                    n_eff += 1; rp = int(m.group(2))
                    for acc2, okflag in revmap.get(fl[3], ()):
                        s2 = U[acc2]["seq"]
                        if 1 <= rp <= len(s2) and s2[rp - 1] == m.group(1): (ok_t if okflag else other_t).add((acc2, rp, m.group(1), m.group(3)))
                if ok_t:
                    c["rescue: tag none or -1"] += 1; c[canon + " rescue: tag none or -1"] += 1
                    c["rescue: one target" if len(ok_t) == 1 else "rescue: several different targets"] += 1
                    if len({(a, p) for a, p, _r, _a in ok_t}) > 1: c["rescue: targets differ in entry or position"] += 1
                    # literal reading of rule 5A: the transcript must be listed for the row's own entry (its SWISSPROT accession)
                    if not acc: c["rescue: row has no SWISSPROT accession (target found through the transcript only)"] += 1
                    elif any(a == acc for a, _p, _r, _a in ok_t): c["rescue: target in the row's own entry"] += 1
                    else: c["rescue: target only in a different entry than the row's SWISSPROT"] += 1
                elif other_t: c["rescue: only through other isoform tags (not used)"] += 1
                elif n_eff == 0: c["not rescuable: no missense effect with a parsable protein change"] += 1
                else: c["not rescuable: no listed transcript with matching residue"] += 1
        c["rows"] = rows
        first = bcs.most_common(1)[0][0] if bcs else ""
        m = BC.match(first) if first else None
        meta.append({"file_id": fid, "case_id": cp[fid]["cases.0.case_id"], "project": cp[fid]["cases.0.project.project_id"], "submitter": cp[fid]["cases.0.submitter_id"],
                     "barcode": first, "n_barcodes": len(bcs), "type": m.group(4) if m else "", "vial": m.group(5) if m else "", "analyte": m.group(7) if m else "",
                     "rows": rows, "missense": c["missense"], "case_mismatch": case_mismatch, "barcode_ok": bool(m) if first else None})
        if n % 1000 == 0 or n == len(man): log(f"  read {n} of {len(man)} files, {round(time.time() - t0)} s")
    # ---- file level checks ----
    log(f"files read: {len(meta)} | distinct header layouts: {len(hdr_variants)} | columns in the first layout: {len(next(iter(hdr_variants)))} {same(len(next(iter(hdr_variants))), REC_MAF['columns'])}")
    tot_rows = sum(m["rows"] for m in meta); empty = sum(1 for m in meta if m["rows"] == 0)
    log(f"data rows: {tot_rows} {same(tot_rows, REC_MAF['rows']) if not LIMIT else '(test run)'} | files without data rows: {empty} {same(empty, REC_MAF['empty_files']) if not LIMIT else '(test run)'}")
    log("Variant_Classification (all rows):", dict(vc_total.most_common()))
    ne = [m for m in meta if m["rows"]]
    log(f"files with more than one tumour barcode: {sum(1 for m in ne if m['n_barcodes'] > 1)} | barcode not in TCGA format: {sum(1 for m in ne if not m['barcode_ok'])} | rows whose case_id differs from the case-project table: {sum(m['case_mismatch'] for m in ne)} | submitter in barcode differs from table: {sum(1 for m in ne if m['barcode'][:12] != m['submitter'])}")
    log("sample type (barcode characters 14-15), all non-empty files:", dict(sorted(collections.Counter(m["type"] for m in ne).items())))
    log("vial letters:", dict(sorted(collections.Counter(m["vial"] for m in ne).items())), "| analyte letters:", dict(sorted(collections.Counter(m["analyte"] for m in ne).items())))
    log("LAML sample types [rule 4A says 03 for LAML, from memory]:", dict(collections.Counter(m["type"] for m in ne if m["project"] == "TCGA-LAML")))
    log("SKCM sample types [rule 4A says most are metastatic, from memory]:", dict(collections.Counter(m["type"] for m in ne if m["project"] == "TCGA-SKCM")))
    # ---- file choice (rule 4A) ----
    bycase = collections.defaultdict(list)
    for m in meta: bycase[m["case_id"]].append(m)
    chosen = {}; step = collections.Counter(); sel_rows = []; gaps = collections.Counter()
    for case, fl in sorted(bycase.items()):
        proj = fl[0]["project"]; prim = "03" if proj == "TCGA-LAML" else "01"
        cand = [m for m in fl if m["rows"] > 0]
        if not cand: cand = fl; gaps["case has only files without data rows"] += 1
        n0 = len(cand); decided = "only one non-empty file" if len(cand) == 1 else None
        if len(cand) == 1: s1 = cand
        else:
            s1 = [m for m in cand if m["type"] == prim]
            if not s1: gaps["several files, none of the primary type (rule does not say; all kept)"] += 1; s1 = cand
        if len(s1) == 1 and decided is None: decided = "sample type"
        s2 = [m for m in s1 if m["vial"] == "A"] or s1
        if len(s2) == 1 and decided is None: decided = "vial A"
        s3 = [m for m in s2 if m["analyte"] == "D"] or s2
        if len(s3) == 1 and decided is None: decided = "analyte D over W"
        if len({m["analyte"] for m in s2} - {"D", "W"}) and len(s2) > 1: gaps["tied files with an analyte other than D or W (rule silent)"] += 1
        pick = sorted(s3, key=lambda m: m["file_id"])[0]
        if decided is None: decided = "lowest file id"; gaps["decided by lowest file id: same barcode" if len({m["barcode"] for m in s3}) == 1 else "decided by lowest file id: different tumour samples"] += 1
        step[decided] += 1; chosen[case] = pick
        sel_rows.append([case, proj, len(fl), n0, pick["file_id"], pick["barcode"], pick["type"], pick["vial"], pick["analyte"], decided, pick["missense"]])
    log(f"cases: {len(bycase)} | chosen file per case: {len(chosen)} | decided by: {dict(step)}")
    log("rule 4A gaps and ties:", dict(gaps))
    log("sample type of the chosen files:", dict(sorted(collections.Counter(m["type"] for m in chosen.values()).items())), "| by project where type is not 01:",
        dict(sorted(collections.Counter((m["project"], m["type"]) for m in chosen.values() if m["type"] not in ("01",)).items())))
    log(f"chosen files without data rows: {sum(1 for m in chosen.values() if m['rows'] == 0)} | cases with more than one non-empty primary file (re-sequencing or different samples): {sum(1 for fl in bycase.values() if len([m for m in fl if m['rows'] and m['type'] == ('03' if fl[0]['project'] == 'TCGA-LAML' else '01')]) > 1)}")
    mc = [m["missense"] for m in chosen.values()]
    log(f"missense rows in chosen files: {sum(mc)} | cases with more than 1000 missense rows (hypermutated by the fixed threshold): {sum(1 for x in mc if x > 1000)} | by project:",
        dict(sorted(collections.Counter(m["project"] for m in chosen.values() if m["missense"] > 1000).items())))
    os.makedirs(H.TABLES, exist_ok=True); os.makedirs(H.PROC_OUT, exist_ok=True)
    with open(os.path.join(H.TABLES, "P2_05_case_file_selection.tsv"), "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["case_id", "project", "n_files", "n_candidates", "chosen_file", "barcode", "sample_type", "vial", "analyte", "decided_by", "missense_rows"]); w.writerows(sel_rows)
    cols = ["file_id", "case_id", "project", "submitter", "barcode", "n_barcodes", "type", "vial", "analyte", "rows", "missense", "case_mismatch", "barcode_ok"]
    with open(os.path.join(H.PROC_OUT, "P2_05_maf_file_summary.tsv"), "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(cols + ["pass", "fail", "rescue_ok"])
        for m in meta: c = reasons[m["file_id"]]; w.writerow([m[k] for k in cols] + [c["pass"], c["fail"], c["rescue: tag none or -1"]])
    log("tables written: results\\tables\\P2_05_case_file_selection.tsv, data\\processed\\P2_05_maf_file_summary.tsv")
    # ---- transcript rule (rule 5A), all files and chosen files ----
    def summ(label, fids):
        t = collections.Counter()
        for fid in fids: t.update(reasons[fid])
        log(f"--- rule 5A counts, {label} ---")
        log(f"missense rows {t['missense']} | canonical {t['canonical missense']} | non-canonical {t['non-canonical missense']} | rows with several SWISSPROT accessions (&): {t['SWISSPROT has several accessions (&)']}")
        for k in ("canonical", "non-canonical"):
            log(f"  {k}: pass {t[k + ' pass']} of {t[k + ' missense']} = {pct(t[k + ' pass'], t[k + ' missense'])} | fail {t[k + ' fail']} | failed rows rescuable (tag none or -1): {t[k + ' rescue: tag none or -1']}")
        log(f"  all: pass {t['pass']} = {pct(t['pass'], t['missense'])} | fail {t['fail']}")
        log("  failure reasons:", {k[6:]: v for k, v in sorted(t.items()) if k.startswith("fail: ")})
        f = t["fail"]
        log("  rescuable rows by entry (literal reading of rule 5A = only 'row's own entry'):", {k[8:]: v for k, v in sorted(t.items()) if k.startswith("rescue: ") and ("entry" in k or "SWISSPROT" in k) and "differ in entry or position" not in k})
        log(f"  failed rows: rescuable with tag none or -1 {t['rescue: tag none or -1']} ({pct(t['rescue: tag none or -1'], f)}; one target {t['rescue: one target']}, several different targets {t['rescue: several different targets']}, of which differing in entry or position {t['rescue: targets differ in entry or position']}) | only other isoform tags {t['rescue: only through other isoform tags (not used)']} | not rescuable {t['not rescuable: no missense effect with a parsable protein change'] + t['not rescuable: no listed transcript with matching residue']} ({ {k: v for k, v in t.items() if k.startswith('not rescuable')} })")
    summ("all files", [m["file_id"] for m in meta])
    summ("chosen files only (one per case)", [m["file_id"] for m in chosen.values()])
    log("Pilot sample of 400 files (seed 20800000) for comparison: pass 93.5% canonical (32496/34765), 92.2% non-canonical (11856/12859); of 3272 failed rows, 1359 (41.5%) rescuable with tag none or -1, 848 only other isoform tags, 1065 not rescuable")

SECTIONS = {"1": ("inputs", section1), "2": ("uniprot_tables", section2), "3": ("sites_tiers", section3), "4": ("maf_files", section4)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for k in (SECTIONS if which == "all" else [which]):
        name, fn = SECTIONS[k]; H.run_section("P2_05", k, name, fn)
# ---- END OF SCRIPT ----
