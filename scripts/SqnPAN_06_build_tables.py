# SqnPAN_06_build_tables.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-10-24; choices 1A and 2A; analysis plan with Amendments 1 and 2)
# Question: turn the TCGA mutations into three tables for the later scripts: sites (every N-X-S/T sequon with its tier), proteins, and mutations
#   (every missense substitution with its loss/gain class, its role in the Q2 position classes and, for gains, the Amendment 2 plausibility state).
# Order: section 1 control C1 (classifier), section 2 site and protein tables, section 3 mutation table. Section 3 refuses to run unless C1 wrote PASSED.
# C1 passes if: every hand-built case gives the expected lost/gained sequons in BOTH the set-difference classifier and the rule-by-rule version; the two versions
#   agree on 200,000 random windows for seed 31000000 (run twice, identical) and seed 31100000 (once); each of loss, gain and both occurs often enough in the
#   random windows to be tested (>= 5000, >= 5000, >= 50). Published examples (KIT N486D loss, NOTCH1 A465T gain) are checked on the UniProt sequence only; if the
#   reference residue or the sequon context is not as published, that is printed as a finding and does not decide the verdict.
# Section 2 stops if the tier counts differ from the recorded ones (T1 13247, T2 8087, T3 38218, 59552 sequons).
# Section 3 rules: 4A file choice (chosen file) plus the other primary-tumour files of the case flagged as union rows (choice 1A); 5A transcript rule;
#   Amendment 1 (rows without a SWISSPROT accession rescued through all_effects with the gene-symbol check; rows with an accession are rescued only onto their own
#   entry, the literal 5A reading; rows with several different targets are dropped and counted); one row per case, protein, position, substitution;
#   hypermutated = more than 1,000 missense rows in the chosen file (flag only). Rescued rows are flagged (only_rescued = every source row was rescued).
# Printed in section 3 (choice 2A): pipeline counts, class counts pooled over all tiers for the chosen-file mutations, the three plausibility states of gains, the
#   gene-symbol check. NOT printed: counts by tier, gene, cancer type or site (the raw material of Q2; first read in SqnPAN_09 after C2 and C3).
# Seeds: base 31000000 (C1 cross-check) and 31100000, at least 100,000 apart (checked in the code).
# Usage: python scripts/SqnPAN_06_build_tables.py [all|1|2|3] [limit=N]     (limit = first N cases, quick test runs only)
# Logs: results\logs\P2_06_section<k>_<name>.txt. Outputs in data\processed: P2_06_C1_verdict.txt, P2_06_sites.tsv, P2_06_proteins.tsv,
#   P2_06_cases.tsv, P2_06_mutations.tsv.gz.
import sys, os, re, csv, gzip, time, random, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import SqnPAN_05_helpers as H
import SqnPAN_06_helpers as P

LIMIT = None
for a in sys.argv[2:]:
    if a.startswith("limit="): LIMIT = int(a[6:])
SEEDS = (31000000, 31100000)
VERDICT = os.path.join(H.PROC_OUT, "P2_06_C1_verdict.txt")
SITES = os.path.join(H.PROC_OUT, "P2_06_sites.tsv"); PROTS = os.path.join(H.PROC_OUT, "P2_06_proteins.tsv")
CASES = os.path.join(H.PROC_OUT, "P2_06_cases.tsv"); MUTS = os.path.join(H.PROC_OUT, "P2_06_mutations.tsv.gz")
REC_TIERS = {"T1": 13247, "T2": 8087, "T3": 38218, "all_sequons": 59552}
HYPER = 1000
PUBLISHED = {("P10721", 486, "N", "D"), ("P46531", 465, "A", "T")}   # KIT N486D, NOTCH1 A465T (Fan et al. 2018), for the presence check only

def tsv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(rows)

# ======================= Section 1: control C1 =======================
HAND = [  # (name, sequence, position, ref, alt, expected lost sequon N positions, expected gained sequon N positions)
    ("loss: N changed to D (N-G-T)", "AANGTAA", 3, "N", "D", [3], []), ("loss: N changed to Q", "AANGTAA", 3, "N", "Q", [3], []),
    ("loss: X changed to P", "AANGTAA", 4, "G", "P", [3], []), ("loss: T changed to A", "AANGTAA", 5, "T", "A", [3], []),
    ("loss: T changed to P (non-S/T)", "AANGTAA", 5, "T", "P", [3], []), ("no loss: T changed to S", "AANGTAA", 5, "T", "S", [], []),
    ("no change: X changed to A", "AANGTAA", 4, "G", "A", [], []), ("no change: X changed to N", "AANGTAA", 4, "G", "N", [], []),
    ("no change: residue outside the sequon", "AANGTAA", 1, "A", "G", [], []), ("no loss: S changed to T", "AANGSAA", 5, "S", "T", [], []),
    ("N-X-C is not a sequon: N changed to D", "AANGCAA", 3, "N", "D", [], []), ("gain: C changed to S in N-G-C", "AANGCAA", 5, "C", "S", [], [3]),
    ("gain: C changed to T in N-G-C", "AANGCAA", 5, "C", "T", [], [3]), ("none: C changed to A in N-G-C", "AANGCAA", 5, "C", "A", [], []),
    ("N-P-T is not a sequon: N changed to D", "AANPTAA", 3, "N", "D", [], []), ("gain: P at position 2 changed to A", "AANPTAA", 4, "P", "A", [], [3]),
    ("gain: P at position 2 changed to G", "AANPTAA", 4, "P", "G", [], [3]), ("none: T changed to A in N-P-T", "AANPTAA", 5, "T", "A", [], []),
    ("none: P changed to A, no S/T at +2", "AANPAAA", 4, "P", "A", [], []), ("gain: new N (A to N before G-T)", "AAAGTAA", 3, "A", "N", [], [3]),
    ("none: new N with P at X", "AAAPTAA", 3, "A", "N", [], []), ("none: new N with no S/T at +2", "AAAGAAA", 3, "A", "N", [], []),
    ("gain: new T at +2", "AANGAAA", 5, "A", "T", [], [3]), ("gain: new S at +2", "AANGAAA", 5, "A", "S", [], [3]),
    ("none: A changed to G at +2", "AANGAAA", 5, "A", "G", [], []), ("none: new T at +2 with P at X", "AANPAAA", 5, "A", "T", [], []),
    ("loss at the N-terminus", "NGT", 1, "N", "D", [1], []), ("loss at the C-terminus (T changed to A)", "NGT", 3, "T", "A", [1], []),
    ("gain at the very start (A to T at residue 3)", "NGA", 3, "A", "T", [], [1]), ("none: new N too close to the C-terminus", "AAAG", 3, "A", "N", [], []),
    ("none: two-residue protein", "NG", 2, "G", "T", [], []), ("none: new N at the last residue", "AAAA", 4, "A", "N", [], []),
    ("both: one sequon lost and another gained (T to N)", "NATGS", 3, "T", "N", [1], [3]),
    ("overlapping sequons: N1 changed, only sequon 1 lost", "NNST", 1, "N", "D", [1], []), ("overlapping sequons: S3 changed, sequon 1 lost, sequon 2 stays", "NNST", 3, "S", "A", [1], []),
    ("overlapping sequons: N2 changed to D, only sequon 2 lost", "NNST", 2, "N", "D", [2], []), ("overlapping sequons: T4 changed to P, sequon 2 lost", "NNST", 4, "T", "P", [2], []),
]

def random_check(seed, n_draws=200000):
    rng = random.Random(seed); pool = "NNNNPPSSSTTTAAGGCDEKLQ"; res = collections.Counter(); examples = []
    for _ in range(n_draws):
        n = rng.randint(6, 12); seq = "".join(rng.choice(pool) for _ in range(n)); pos = rng.randint(1, n); ref = seq[pos - 1]
        alt = rng.choice([a for a in (P.AA if rng.random() < 0.5 else "NSTPA") if a != ref]); a = P.classify_change(seq, pos, ref, alt); b = P.classify_literal(seq, pos, ref, alt)
        if a != b:
            res["disagreements"] += 1
            if len(examples) < 5: examples.append((seq, pos, ref, alt, a, b))
        res["both" if a[0] and a[1] else ("loss" if a[0] else ("gain" if a[1] else "none"))] += 1
    return res, examples

def section1(log):
    assert abs(SEEDS[0] - SEEDS[1]) >= 100000, "seed bases closer than 100,000"
    log("--- hand-built cases (every rule of analysis plan section 5, protein ends, overlapping and double events) ---")
    bad = 0
    for name, seq, pos, ref, alt, el, eg in HAND:
        a = P.classify_change(seq, pos, ref, alt); b = P.classify_literal(seq, pos, ref, alt); ok = a == b == (el, eg)
        if not ok: bad += 1
        log(f"  {'ok  ' if ok else 'FAIL'} {name}: {seq} pos {pos} {ref}->{alt} | expected lost {el} gained {eg} | set-difference {a} | rule-by-rule {b}")
    log(f"hand-built cases: {len(HAND)} | wrong: {bad}")
    log("--- random windows, set-difference classifier against rule-by-rule version ---")
    runs = {}
    for label, seed in (("seed 31000000 run 1", SEEDS[0]), ("seed 31000000 run 2", SEEDS[0]), ("seed 31100000", SEEDS[1])):
        res, ex = random_check(seed); runs[label] = res; log(f"  {label}: {dict(res)}")
        for e in ex: log("    disagreement example:", e)
    same_runs = runs["seed 31000000 run 1"] == runs["seed 31000000 run 2"]
    cover = all(r["loss"] >= 5000 and r["gain"] >= 5000 and r["both"] >= 50 for r in runs.values())
    nodis = all(r["disagreements"] == 0 for r in runs.values())
    log(f"same seed twice identical: {same_runs} | no disagreements in any run: {nodis} | coverage (loss >= 5000, gain >= 5000, both >= 50) in every run: {cover}")
    log("--- published examples (UniProt 2026_03 sequence only; findings, not part of the verdict) ---")
    U = H.load_uniprot_main()
    for label, acc, pos, ref, alt, want in (("KIT N486D (Fan et al. 2018, loss)", "P10721", 486, "N", "D", "loss"), ("NOTCH1 A465T (Fan et al. 2018, gain)", "P46531", 465, "A", "T", "gain")):
        s = U[acc]["seq"]; found = s[pos - 1] if len(s) >= pos else "(beyond sequence)"
        if found != ref: log(f"  {label}: UniProt {acc} residue {pos} is {found}, not {ref}: finding (sequence differs from the one used in the paper)"); continue
        lost, gained = P.classify_change(s, pos, ref, alt); log(f"  {label}: UniProt {acc}, window {s[max(pos - 4, 0):pos + 3]}, lost {lost}, gained {gained} | classifier says {'loss' if lost else ('gain' if gained else 'none')}, paper says {want}")
    verdict = "PASSED" if (bad == 0 and same_runs and cover and nodis) else "FAILED"
    log(f"C1 VERDICT: {verdict}")
    os.makedirs(H.PROC_OUT, exist_ok=True)
    with open(VERDICT, "w") as o: o.write(f"{verdict} {time.strftime('%Y-%m-%d %H:%M:%S')} hand cases {len(HAND)}, wrong {bad}, random disagreements {sum(r['disagreements'] for r in runs.values())}\n")
    log(f"verdict file: {VERDICT}")

# ======================= Section 2: site and protein tables =======================
def section2(log):
    U = H.load_uniprot_main(); F = H.load_uniprot_features(); topo = P.load_topology(F)
    use = {r["accession"]: r for r in csv.DictReader(open(H.USE_TSV, newline=""), delimiter="\t")}
    T = P.build_tiers(U); tier_of = T["tier_of"]; c = T["counts"]
    log("counts:", c)
    for k, v in REC_TIERS.items(): log(f"  {k}: {c[k]} (recorded in SqnPAN_05: {v}) {'same' if c[k] == v else 'DIFFERENT'}")
    log(f"  atlas reviewed sites {c['atlas_reviewed_sites']} (SqnPAN_05: 12959) | atlas sites with N-X-S/T {c['atlas_valid']} (SqnPAN_05: 12523) | UniProt N-linked sites not N-X-S/T {c['uniprot_annotated_not_sequon']} (SqnPAN_05: 120)")
    if any(c[k] != v for k, v in REC_TIERS.items()) or not c["tiers_disjoint_and_complete"]: raise RuntimeError("tier counts differ from the recorded ones: stop and decide")
    notes_case = sum(1 for d in topo.values() for _a, _b, n in d["topo"] if n.lower().startswith("cytoplasmic") != n.startswith("Cytoplasmic"))
    log(f"TOPO_DOM notes that start with 'cytoplasmic' only in lower case or other case (would change the literal Amendment 2 definition): {notes_case}")
    hasT1 = {a for (a, p), t in tier_of.items() if t == "T1"}
    rows = []; sig_n = 0
    for (acc, p), t in sorted(tier_of.items()):
        u = U[acc]; s = u["seq"]; ur = use.get(acc, {}); nf = ur.get("n_fragments", "0")
        codes = "+".join(x.replace("ECO:", "") for x in sorted(T["codes"].get((acc, p), ())))
        rows.append([acc, u["gene"], p, s[p - 1:p + 2], t, codes, "yes" if (acc, p) in T["annotated"] else "no", "yes" if (acc, p) in T["atlas"] else "no",
                     "yes" if p + 3 <= len(s) and s[p + 2] == "P" else "no", len(s), ur.get("status", "not_in_use_table"), nf, "yes" if nf not in ("0", "1", "") else "no",
                     P.topology_category(topo[acc], p), "yes" if topo[acc]["sig"] else "no", "yes" if acc in hasT1 else "no"])
    os.makedirs(H.PROC_OUT, exist_ok=True)
    tsv(SITES, ["accession", "gene", "N_position", "motif", "tier", "uniprot_eco_codes", "uniprot_annotated", "in_atlas", "proline_at_plus3", "protein_length", "alphafold_status",
                "n_fragments", "fragment_model", "topology_at_N", "protein_has_signal_peptide", "protein_has_T1"], rows)
    seqN = collections.defaultdict(set); per = collections.defaultdict(collections.Counter)
    for (acc, p), t in tier_of.items(): seqN[acc].add(p); per[acc][t] += 1
    prow = []
    for acc, u in U.items():
        s = u["seq"]; nN = s.count("N"); nS = s.count("S"); nT = s.count("T"); k = len(seqN.get(acc, ())); ur = use.get(acc, {}); nf = ur.get("n_fragments", "0")
        prow.append([acc, u["name"], u["gene"], len(s), ur.get("status", "not_in_use_table"), nf, "yes" if nf not in ("0", "1", "") else "no", "yes" if topo[acc]["sig"] else "no",
                     "yes" if (topo[acc]["topo"] or topo[acc]["tm"]) else "no", per[acc]["T1"], per[acc]["T2"], per[acc]["T3"], nN, nS, nT, nN - k, nS + nT - k])
    tsv(PROTS, ["accession", "entry_name", "gene", "length", "alphafold_status", "n_fragments", "fragment_model", "has_signal_peptide", "has_topology", "n_T1_sequons", "n_T2_sequons",
                "n_T3_sequons", "n_N", "n_S", "n_T", "n_N_in_no_sequon", "n_ST_in_no_sequon"], prow)
    log(f"site table: {len(rows)} rows -> {SITES} | md5 {H.md5_file(SITES)}")
    log(f"protein table: {len(prow)} rows (UniProt entries {len(U)}) -> {PROTS} | md5 {H.md5_file(PROTS)}")
    log(f"proteins with at least one T1 sequon: {len(hasT1)} | check: sum of n_N_in_no_sequon + sequon N = all N: {sum(r[15] for r in prow) + len(tier_of) == sum(r[12] for r in prow)}")

# ======================= Section 3: mutation table =======================
def section3(log):
    if not os.path.exists(VERDICT) or not open(VERDICT).read().startswith("PASSED"): raise RuntimeError("control C1 has not PASSED (verdict file missing or FAILED): section 3 does not run")
    log("C1 verdict file:", open(VERDICT).read().strip())
    if LIMIT: log(f"TEST RUN: only the first {LIMIT} cases")
    U = H.load_uniprot_main(); topo = P.load_topology(H.load_uniprot_features()); revmap = P.build_revmap(U)
    tier_of = {(r["accession"], int(r["N_position"])): r["tier"] for r in csv.DictReader(open(SITES, newline=""), delimiter="\t")}
    hasT1 = {r["accession"] for r in csv.DictReader(open(SITES, newline=""), delimiter="\t") if r["tier"] == "T1"}
    gene_odd = sum(1 for u in U.values() if (" " in u["gene"]) or (";" in u["gene"]))
    log(f"site table sequons: {len(tier_of)} | proteins with T1: {len(hasT1)} | UniProt entries whose primary gene field holds more than one word: {gene_odd}")
    sel = {r["case_id"]: r for r in csv.DictReader(open(P.SELECTION, newline=""), delimiter="\t")}
    fsum = collections.defaultdict(list)
    for r in csv.DictReader(open(P.FILESUM, newline=""), delimiter="\t"): fsum[r["case_id"]].append(r)
    fname = {r["id"]: r["filename"] for r in csv.DictReader(open(H.MANIFEST, newline=""), delimiter="\t")}
    subm = {r["cases.0.case_id"]: r["cases.0.submitter_id"] for r in csv.DictReader(open(H.CASE_PROJECT, newline=""), delimiter="\t")}
    cases = sorted(sel)
    if LIMIT: cases = cases[:LIMIT]
    CT = {"chosen": collections.Counter(), "union": collections.Counter()}
    tally = collections.Counter(); guard = collections.Counter(); sym = collections.Counter(); seen = {}
    tally["chosen-file mutations with more than one source row in the chosen file (same case, protein, position, substitution)"] += 0
    best_rank = {"plausible": 0, "unknown": 1, "implausible": 2}
    cf = open(CASES, "w", newline="", encoding="utf-8"); cw = csv.writer(cf, delimiter="\t", lineterminator="\n")
    cw.writerow(["case_id", "submitter", "project", "chosen_file", "missense_rows_chosen", "hypermutated", "n_files_used", "n_mutations", "n_mutations_in_chosen"])
    mf = gzip.open(MUTS, "wt", newline="", encoding="utf-8", compresslevel=6)
    mw = csv.writer(mf, delimiter="\t", lineterminator="\n")
    mw.writerow(["case_id", "submitter", "project", "accession", "gene", "position", "ref", "alt", "n_rows_all_files", "n_rows_chosen_file", "in_chosen_file", "only_rescued", "canonical_any", "hypermutated", "res_role",
                 "lost_sequon_N", "lost_tiers", "gained_sequon_N", "gain_states", "gain_state_best", "class"])
    t0 = time.time(); need = ["Variant_Classification", "Hugo_Symbol", "SWISSPROT", "Protein_position", "Amino_acids", "CANONICAL", "all_effects"]
    for ci, case in enumerate(cases, 1):
        s = sel[case]; proj = s["project"]; prim = "03" if proj == "TCGA-LAML" else "01"; chosen = s["chosen_file"]
        extra = [r["file_id"] for r in fsum[case] if int(r["rows"]) > 0 and r["type"] == prim and r["file_id"] != chosen]
        hyper = int(s["missense_rows"]) > HYPER; local = {}
        for fid, scope in [(chosen, "chosen")] + [(x, "union") for x in sorted(extra)]:
            c = CT[scope]; c["files read"] += 1; hdr = None
            with gzip.open(os.path.join(H.MAF_DIR, fid, fname[fid]), "rt", encoding="utf-8", newline="") as f:
                for line in f:
                    if line[0] == "#": continue
                    if hdr is None: hdr = line.rstrip("\n").split("\t"); ix = {k: hdr.index(k) for k in need}; continue
                    if "Missense_Mutation" not in line: continue
                    v = line.rstrip("\n").split("\t")
                    if v[ix["Variant_Classification"]] != "Missense_Mutation": continue
                    c["missense rows"] += 1; canon = v[ix["CANONICAL"]] == "YES"; sw = v[ix["SWISSPROT"]]; acc = sw.split(".")[0].split("&")[0]; hugo = v[ix["Hugo_Symbol"]]; reason = None
                    if not acc: reason = "no SWISSPROT accession"
                    elif acc not in U: reason = "accession not in UniProt reviewed 2026_03"
                    else:
                        try: pos = int(v[ix["Protein_position"]].split("/")[0]); tot = int(v[ix["Protein_position"]].split("/")[1]); ref, alt = v[ix["Amino_acids"]].split("/")
                        except Exception: reason = "position or amino acids unparsable"
                        if reason is None:
                            sq = U[acc]["seq"]; lenok = tot == len(sq); resok = 1 <= pos <= len(sq) and len(ref) == 1 and sq[pos - 1] == ref
                            reason = None if (lenok and resok) else "length or reference residue differs"
                    how = None
                    if reason is None:
                        if len(alt) != 1 or alt not in P.AA or alt == ref: c["pass, but alt not a single different standard residue (row skipped)"] += 1; continue
                        c["pass"] += 1; c["pass: Hugo_Symbol equals the UniProt primary gene name" if hugo == U[acc]["gene"] else "pass: Hugo_Symbol differs from the UniProt primary gene name"] += 1
                        key = (acc, pos, ref, alt); how = "direct"
                    else:
                        c["fail"] += 1; c["fail: " + reason] += 1; targets = set(); other = False
                        for e in v[ix["all_effects"]].split(";"):
                            fl = e.split(",")
                            if len(fl) < 4 or "missense_variant" not in fl[1]: continue
                            m = re.match(r"^p\.([A-Z])(\d+)([A-Z])$", fl[2])
                            if not m: continue
                            rp = int(m.group(2))
                            for acc2, okflag in revmap.get(fl[3], ()):
                                s2 = U[acc2]["seq"]
                                if 1 <= rp <= len(s2) and s2[rp - 1] == m.group(1):
                                    if okflag: targets.add((acc2, rp, m.group(1), m.group(3)))
                                    else: other = True
                        if acc:
                            own = {t for t in targets if t[0] == acc}
                            if len(own) == 1: how = "rescued"; key = next(iter(own)); c["rescued: row with accession, target in its own entry"] += 1
                            elif len(own) > 1: c["dropped: row with accession, several different targets in its own entry"] += 1
                            elif targets: c["dropped: row with accession, target only in a different entry"] += 1
                            else: c["dropped: row with accession, " + ("target only through other isoform tags" if other else "no listed transcript with matching residue")] += 1
                        else:
                            if not targets: c["dropped: row without accession, " + ("target only through other isoform tags" if other else "no listed transcript with matching residue")] += 1
                            elif len(targets) > 1:
                                c["dropped: row without accession, several different targets"] += 1
                                if sum(1 for t in targets if U[t[0]]["gene"] == hugo) == 1: c["   of these, exactly one target has the matching gene symbol (information only, row stays dropped)"] += 1
                            else:
                                t = next(iter(targets))
                                if U[t[0]]["gene"] == hugo: how = "rescued"; key = t; c["rescued: row without accession, one target, Hugo_Symbol equals its primary gene name"] += 1
                                else: c["dropped: row without accession, one target but Hugo_Symbol differs from its primary gene name (gene-symbol check)"] += 1
                    if how is None: continue
                    if key[3] == key[2] or key[3] not in P.AA: c["skipped: alt not a different standard residue after rescue"] += 1; continue
                    e = local.setdefault(key, [0, 0, 0, 0])      # rows in all files, rows in the chosen file, direct rows, canonical any
                    e[0] += 1; e[1] += (scope == "chosen"); e[2] += (how == "direct"); e[3] |= canon
        n_in = 0
        for (acc, pos, ref, alt), e in sorted(local.items()):
            sq = U[acc]["seq"]; lost, gained = P.classify_change(sq, pos, ref, alt)
            if (lost, gained) != P.classify_literal(sq, pos, ref, alt): guard["classifier disagreements"] += 1
            lt = [tier_of.get((acc, p), "?") for p in lost]
            if "?" in lt: guard["lost sequon not in site table"] += 1
            gs = [P.plausibility(topo[acc], p, acc in hasT1) for p in gained]; gbest = min(gs, key=lambda x: best_rank[x]) if gs else ""
            if ref == "N": role = "seqN_" + tier_of[(acc, pos)] if (acc, pos) in tier_of else "nonseqN"
            elif ref in "ST": role = "seqST_" + tier_of[(acc, pos - 2)] if (acc, pos - 2) in tier_of else "nonseqST"
            else: role = "other"
            if role.startswith("seqN_") and pos not in lost: guard["N of a sequon changed but not lost"] += 1
            if role == "nonseqN" and pos in lost: guard["non-sequon N lost a sequon"] += 1
            if role.startswith("seqST_") and alt not in "ST" and (pos - 2) not in lost: guard["S/T of a sequon changed to non-S/T but not lost"] += 1
            cls = "loss_and_gain" if lost and gained else ("loss" if lost else ("gain" if gained else "none"))
            mw.writerow([case, subm[case], proj, acc, U[acc]["gene"], pos, ref, alt, e[0], e[1], int(e[1] > 0), int(e[2] == 0), e[3], int(hyper), role, ";".join(map(str, lost)), ";".join(lt),
                         ";".join(map(str, gained)), ";".join(gs), gbest, cls])
            if e[1] > 0:
                n_in += 1; tally["chosen-file mutations (unique per case, protein, position, substitution)"] += 1; tally["class: " + cls] += 1
                if gbest: tally["gain state (best state per mutation): " + gbest] += 1
            else: tally["union-only mutations (not in the chosen file)"] += 1
            if e[2] == 0: tally["mutations where every source row was rescued"] += 1
            if e[1] > 1: tally["chosen-file mutations with more than one source row in the chosen file (same case, protein, position, substitution)"] += 1
            if e[0] > e[1] and e[1] > 0: tally["chosen-file mutations also seen in another primary-tumour file of the case"] += 1
            if (acc, pos, ref, alt) in PUBLISHED: seen[(acc, pos, ref, alt)] = cls
            if hyper: tally["mutations in hypermutated cases"] += 1
        cw.writerow([case, subm[case], proj, chosen, s["missense_rows"], int(hyper), 1 + len(extra), len(local), n_in]); tally["cases"] += 1; tally["hypermutated cases"] += int(hyper)
        if ci % 1000 == 0 or ci == len(cases): log(f"  {ci} of {len(cases)} cases, {round(time.time() - t0)} s")
    mf.close(); cf.close()
    for scope in ("chosen", "union"):
        log(f"--- source rows, {scope} files ---")
        for k, v in sorted(CT[scope].items()): log(f"  {k}: {v}")
    c = CT["chosen"]
    if not LIMIT: log(f"SqnPAN_05 (chosen files): missense 1517858, pass 1417588, fail 100270 | here: missense {c['missense rows']}, pass {c['pass']}, fail {c['fail']} | SqnPAN_05 rescuable with tag none or -1: own entry 16807, no accession 31400, other entry 442, several different targets 324 (counted before the own-entry restriction)")
    log("--- mutation table ---")
    for k, v in sorted(tally.items()): log(f"  {k}: {v}")
    log("--- guards (all must be 0) ---", dict(guard) if guard else "all 0")
    log("--- gene-symbol check: rescued rows without accession (match / mismatch) are in the source-row counts above; baseline = passing rows ---")
    for scope in ("chosen", "union"):
        a = CT[scope]["pass: Hugo_Symbol equals the UniProt primary gene name"]; b = CT[scope]["pass: Hugo_Symbol differs from the UniProt primary gene name"]
        log(f"  {scope}: passing rows with Hugo_Symbol equal to the UniProt primary gene name {a}, different {b} ({100.0 * b / max(a + b, 1):.2f}% different)")
    for key, name in (((("P10721", 486, "N", "D")), "KIT N486D"), ((("P46531", 465, "A", "T")), "NOTCH1 A465T")):
        log(f"  published example {name}: {'present in the mutation table, class ' + str(seen[key]) if key in seen and seen[key] else 'not present in the table'}")
    for p in (MUTS, CASES): log(f"written {p} | size {os.path.getsize(p)} | md5 {H.md5_file(p)}")

SECTIONS = {"1": ("controls_C1", section1), "2": ("site_tables", section2), "3": ("mutation_table", section3)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for k in (SECTIONS if which == "all" else [which]):
        name, fn = SECTIONS[k]; H.run_section("P2_06", k, name, fn)
# ---- END OF SCRIPT ----
