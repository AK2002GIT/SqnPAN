# SqnPAN_13_fan_crosscheck.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-10-31; decisions: D1 A = four parts, D2 = no supplementary lists (printed tables only),
#   D3 A = three-level "reproduced" rule). Written before any TCGA value was compared with Fan et al.; no result of the cross-check has been read.
# LABELLED CROSS-CHECK, NOT BLIND: Fan Y, Hu Y, Yan C, Goldman R, Pan Y, Mazumder R, Dingerdissen HM. Loss and gain of N-linked glycosylation sequons due to
#   single-nucleotide variation in cancer. Sci Rep 2018;8:4322. doi 10.1038/s41598-018-22345-2, PMID 29531238, PMCID PMC5847511. Their data (COSMIC v73, IntOGen 2014-12,
#   ICGC v0.10a, TCGA 2015-01-27, ClinVar 2015-02; UniProt 2015_01 plus HPRD, dbPTM, neXtProt, CDD) differ from ours (TCGA masked MAFs from the GDC, UniProt 2026_03,
#   N-GlycositeAtlas). Their cancer types are Disease Ontology terms (overlapping: colon, rectum, large intestine), ours are 33 TCGA projects: counts of cancer types are
#   NOT like for like. Everything here is descriptive; no new test, no seed (no random step).
# Fan tables 4 (13 high-confidence LOG variants) and 5 (GOG variants) were taken on 2026-10-31 from the Europe PMC full-text XML (the PMC page
#   and the shell were blocked), NOT from the PDF, so a transcription slip is possible. The printed Table 5 has 39 rows while their text says 41: the two extra rows are
#   unresolved. Guards below: the reference residue of every variant is checked against UniProt 2026_03, and the gene symbol is compared.
#   Fan numbers quoted in Part 3 and Part 4 come from their text and Table 3 (read via PubMed Central full text by the same date).
# PARTS (D1 A)
#   Part 1  Each Fan variant (13 LOG, 39 GOG, plus CREBBP N1978D named in their text) looked up in our chosen-file mutation table (the Q1 set).
#           Label (D3 A, fixed now): "Reproduced" = the exact variant (accession, position, reference, substitution) is seen in >= 1 case; "Same position, other substitution"
#           = no exact match but >= 1 case with another substitution at that residue; "Not seen"; "Not checkable (reason)" = accession not in UniProt 2026_03 or reference
#           residue differs. Reported beside the label: cases, cases without hypermutated cases, cases in any primary-tumour file, number of TCGA projects (>= 3 = recurrent in
#           TCGA), sequence class from our classifier (loss / gain), lost sequon and its tier, plausibility state of observed gains.
#   Part 2  Position of the lost sequon hit (1 = N, 2 = X to P, 3 = S/T) by tier, main set and without hypermutated cases (same counts as the Q1 position table), with the
#           Q2 observed and expected counts of the same tier (P2_09_Q2_results.tsv, main set, pooled units) for N and S/T. Fan's statement: in somatic LOG variants the third
#           position is hit more often than the first. "Reproduced (raw)" = position 3 count > position 1 count. "Expected from background rates" = the ratio of observed S/T to
#           N events in the eligible Q2 strata is within the Q2 threshold of 1.25-fold of the ratio of their expected counts (|log2 of the ratio of ratios| < 0.32) AND neither the
#           N nor the S/T unit has a Q2 final label Enriched or Depleted. No new test.
#   Part 3  Counts side by side (Fan Table 3 and text against ours: variants, sequons, proteins for LOG and GOG; high-confidence sequons abolished; variants in >= 3 cancer
#           types) and the overlap of the variants seen in >= 3 TCGA projects (our categories loss_T12 = lost sequon of tier T1 or T2, loss_any, gain_pu = gain plausible or
#           unknown, gain_any) with Fan tables 4 and 5.
#   Part 4  Genes named in the Fan text (USH2A, CEACAM5, CEACAM7, KIT, THSD7A, NOTCH1, MUC2, TTN, LRP1B, TP53, LCN8, CREBBP): our counts of distinct loss and gain variants and
#           sequons and the rank of the gene among all proteins, beside the Fan statements.
# CONTROLS (section 1, judged before section 2 runs; section 2 refuses to run unless C-F2 PASSED)
#   C-F0  at least 80% of the 52 table variants are checkable (accession in UniProt 2026_03 and reference residue equal). If not: Part 1 labels carry _table_check_FAILED.
#   C-F1  our classifier (classify_change) reproduces Fan's category (LOG variant = sequon lost, GOG variant = sequon gained) for >= 80% of the checkable variants, and the
#         independent classify_literal agrees with classify_change for every checkable variant. If not: Part 1 labels carry _table_check_FAILED.
#   C-F2  self-consistency: the chosen-file class counts equal those of the SqnPAN_06 run (loss 11,735; gain 18,994; loss_and_gain 187; none 1,429,876; gain states plausible 3,654,
#         unknown 13,630, implausible 1,897); the pooled position counts rebuilt here equal every pooled line of P2_09_Q1_loss_position_in_sequon.tsv (main and nohyper, T1 to T3);
#         every line of P2_09_Q1_recurrent_variants.tsv (n_cases, n_cases_without_hypermutated, n_cancer_types) is rebuilt exactly and the number of variants seen in >= 3 types
#         equals the saved table for each of its five categories. Any difference = FAILED, nothing later is computed.
#   No effect estimate of the real data is printed in section 1 beyond equality checks of the saved tables.
# From-memory items: the Fan numbers typed in this file (marked FAN_ constants) are from the sources above; the gene CREBBP is resolved to its accession through the UniProt
#   gene field, not typed. Any change after seeing results is a dated amendment and is called post hoc.
# Order: section 1 inputs and controls; section 2 Parts 1 to 4. Usage (from the repository folder):  python scripts/SqnPAN_13_fan_crosscheck.py all
# Logs: results\logs\P2_13_section<k>_<name>.txt. Outputs: data\processed\P2_13_controls_verdict.txt; results\tables\P2_13_fan_variants.tsv, P2_13_fan_positions.tsv,
#   P2_13_fan_counts.tsv, P2_13_fan_recurrent_overlap.tsv, P2_13_fan_genes.tsv.
import os, sys, csv, gzip, collections, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import SqnPAN_05_helpers as H
import SqnPAN_06_helpers as H6

REAL_TABLES = os.path.join(H.PAPER2, "results", "tables")
MUTS = os.path.join(H.PROC, "P2_06_mutations.tsv.gz"); SITES = os.path.join(H.PROC, "P2_06_sites.tsv")
Q1_POS = os.path.join(REAL_TABLES, "P2_09_Q1_loss_position_in_sequon.tsv"); Q1_VARS = os.path.join(REAL_TABLES, "P2_09_Q1_recurrent_variants.tsv")
Q2_RES = os.path.join(REAL_TABLES, "P2_09_Q2_results.tsv")
VERDICT = os.path.join(H.PROC_OUT, "P2_13_controls_verdict.txt")
OUT = {k: os.path.join(H.TABLES, "P2_13_" + k + ".tsv") for k in ("fan_variants", "fan_positions", "fan_counts", "fan_recurrent_overlap", "fan_genes")}
MD5 = {H.UNIPROT_MAIN: "af9af22f1e6eea80e96751c93c9da2b3", MUTS: "3be28f9265312957c13ab2f0faa2b2b2", SITES: "34e0a114e0e7a6955418d837e9a5a245",
       Q1_POS: "a917402b84e6500e50ec8928f29a7085", Q1_VARS: "98dddb7348e0ccf0e158fa4b449f8326", Q2_RES: "2e25c9b14f3cc0984341648af07f8698",
       os.path.join(HERE, "SqnPAN_06_helpers.py"): "e8c4aeacb53a52e1a81255e5673cdf05"}
EXPECT_CLASS = {"loss": 11735, "gain": 18994, "loss_and_gain": 187, "none": 1429876}
EXPECT_GAIN = {"plausible": 3654, "unknown": 13630, "implausible": 1897}
RELNAME = {0: "1 (N)", 1: "2 (X to P)", 2: "3 (S/T)"}
VARCATS = ("loss_T1", "loss_T2", "loss_T3", "gain_plausible", "gain_unknown")

# Fan et al. 2018 Table 4 (high-confidence LOG variants in >= 3 cancer types) and Table 5 (somatic-only GOG variants in >= 3 cancers): accession|gene|position|ref|alt|cancer types
FAN_LOG = """
O15244|SLC22A2|74|T|M|pancreatic cancer;brain cancer;colon cancer
O15438|ABCC3|1009|S|F|skin cancer;stomach cancer;liver cancer
O75473|LGR5|77|N|S|colon cancer;lung cancer;large intestine cancer
P10163|PRB4|110|S|P|brain cancer;malignant glioma;head and neck cancer
P10721|KIT|486|N|D|bone cancer;rectum cancer;colon cancer;large intestine cancer
P28472|GABRB3|107|T|M|pharynx cancer;stomach cancer;head and neck cancer
Q96S37|SLC22A12|104|T|M|brain cancer;stomach cancer;colon cancer
Q9NTG1|PKDREJ|297|S|L|skin cancer;pharynx cancer;head and neck cancer
Q9NTG1|PKDREJ|925|S|L|stomach cancer;colon cancer;large intestine cancer
Q9NUN5|LMBRD1|349|S|C|brain cancer;lung cancer;malignant glioma
Q9P121|NTM|46|T|M|pharynx cancer;esophageal cancer;head and neck cancer
Q9UKZ4|TENM1|1759|S|L|cervical cancer;pharynx cancer;head and neck cancer
Q9UPZ6|THSD7A|236|T|M|pharynx cancer;rectum cancer;esophageal cancer;head and neck cancer
"""
FAN_GOG = """
O14522|PTPRT|281|A|T|brain cancer;stomach cancer;malignant glioma
O43526|KCNQ2|785|D|N|pancreatic cancer;breast cancer;prostate cancer
O43699|SIGLEC6|251|A|T|brain cancer;colon cancer;malignant glioma
O60469|DSCAM|213|A|T|uterine cancer;breast cancer;pancreatic cancer;esophageal cancer
O94973|AP2A2|115|A|T|breast cancer;colon cancer;large intestine cancer
P04264|KRT1|248|D|N|skin cancer;ovarian cancer;large intestine cancer
P07357|C8A|169|D|N|skin cancer;uterine cancer;melanoma
P10721|KIT|566|N|S|rectum cancer;kidney cancer;large intestine cancer
P21802|FGFR2|659|K|N|uterine cancer;breast cancer;lung cancer
P21817|RYR1|2861|D|N|oral cavity cancer;pancreatic cancer;uterine cancer;head and neck cancer
P46531|NOTCH1|465|A|T|colon cancer;malignant glioma;esophageal cancer;head and neck cancer;pharynx cancer;brain cancer
P54852|EMP3|42|D|N|oral cavity cancer;pharynx cancer;head and neck cancer
Q02817|MUC2|1750|T|N|breast cancer;liver cancer;thyroid cancer;ovarian cancer;prostate cancer;brain cancer;rectum cancer
Q13002|GRIK2|528|D|N|oral cavity cancer;breast cancer;urinary bladder cancer
Q13349|ITGAD|1070|D|N|uterine cancer;rectum cancer;colon cancer;large intestine cancer
Q4ZHG4|FNDC1|253|K|N|uterine cancer;pancreatic cancer;rectum cancer
Q685J3|MUC17|2784|T|N|ovarian cancer;head and neck cancer;pharynx cancer
Q6P1J6|PLB1|645|D|N|uterine cancer;pharynx cancer;head and neck cancer
Q6UWW8|CES3|161|D|N|skin cancer;colon cancer;large intestine cancer
Q6UX06|OLFM4|372|R|S|pancreatic cancer;breast cancer;esophageal cancer
Q7Z304|MAMDC2|319|D|N|skin cancer;breast cancer;urinary bladder cancer
Q7Z5H5|VN1R4|265|L|S|cervical cancer;pancreatic cancer;head and neck cancer
Q86TH1|ADAMTSL2|44|D|N|stomach cancer;colon cancer;large intestine cancer
Q8IZD2|KMT2E|902|D|N|uterine cancer;colon cancer;large intestine cancer
Q8N158|GPC2|200|D|N|breast cancer;esophageal cancer;large intestine cancer
Q8N8F6|YIPF7|141|D|N|cervical cancer;pharynx cancer;head and neck cancer
Q8NFZ4|NLGN2|542|A|T|oral cavity cancer;pharynx cancer;head and neck cancer
Q8NGZ4|OR2G3|159|H|N|uterine cancer;colon cancer;lung cancer
Q8TC71|SPATA18|404|K|N|uterine cancer;brain cancer;pancreatic cancer;lung cancer
Q8TDM6|DLG5|1799|D|N|uterine cancer;pharynx cancer;head and neck cancer
Q92556|ELMO1|55|D|N|skin cancer;brain cancer;malignant glioma
Q92736|RYR2|898|I|T|pharynx cancer;lung cancer;head and neck cancer
Q96PZ7|CSMD1|3053|D|N|skin cancer;stomach cancer;colon cancer
Q9BXX0|EMILIN2|759|K|N|uterine cancer;lung cancer;urinary bladder cancer
Q9H2B2|SYT4|89|K|N|uterine cancer;breast cancer;colon cancer
Q9H9P2|CHODL|186|P|S|skin cancer;breast cancer;lung cancer
Q9UKJ8|ADAM21|278|D|N|skin cancer;colon cancer;large intestine cancer
Q9UQP3|TNN|1091|D|N|pharynx cancer;lung cancer;head and neck cancer
Q9Y5F1|PCDHB12|556|D|N|uterine cancer;skin cancer;stomach cancer
"""
# Fan numbers from their Table 3 (somatic-only) and text (Part 3), and statements about named genes (Part 4)
FAN_COUNTS = [("LOG variants, somatic-only (all sequons)", 8895, "loss_any", "variants"), ("LOG sequons, somatic-only", 8061, "loss_any", "sequons"),
              ("LOG proteins, somatic-only", 5204, "loss_any", "proteins"), ("GOG variants, somatic-only (string search)", 12939, "gain_any", "variants"),
              ("GOG sequons, somatic-only", 12815, "gain_any", "sequons"), ("GOG proteins, somatic-only", 7356, "gain_any", "proteins")]
FAN_GENES = [("USH2A", "most somatic-only LOGs: 27 LOG variants abolishing 25 high-confidence sequons"), ("CEACAM7", "5 LOG variants at 5 positions, 265 residues"),
             ("CEACAM5", "12 LOG variants at 11 positions (secreted/membrane subset)"), ("KIT", "N486D (LOG) in 4 cancers; N566S (GOG) in 3"),
             ("THSD7A", "T236M (LOG) in 4 cancers"), ("NOTCH1", "A465T (GOG) in 6 cancers"), ("MUC2", "T1750N (GOG) in 7 cancers"),
             ("TTN", "most somatic-only GOGs: 92 GOG variants affecting 90 positions"), ("LRP1B", "24 GOG variants at 23 positions"),
             ("TP53", "15 GOG variants, highest GOG density (one per 27 residues)"), ("LCN8", "3 GOG variants (liver and skin cancer datasets)"),
             ("CREBBP", "N1978D (LOG, ovarian cancer; HGMD-listed)")]

def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f: return list(csv.DictReader(f, delimiter="\t"))
def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(header)
        for r in rows: w.writerow(r)

def check_md5(log):
    bad = 0
    for p, m in MD5.items():
        got = H.md5_file(p); ok = (got == m); bad += (not ok)
        log(f"  md5 {os.path.basename(p)}: {'same' if ok else 'DIFFERENT (' + got + ')'}")
    if bad: raise RuntimeError("an input md5 differs from the recorded value: stop and decide")

def fan_rows(U):
    out = []
    for setname, text in (("LOG", FAN_LOG), ("GOG", FAN_GOG)):
        for line in text.strip().splitlines():
            a, g, p, r, al, c = line.strip().split("|")
            out.append(dict(set=setname, acc=a, gene=g, pos=int(p), ref=r, alt=al, cancers=c, n_fan=len(c.split(";")), text_named=False))
    if (sum(x["set"] == "LOG" for x in out), sum(x["set"] == "GOG" for x in out)) != (13, 39): raise RuntimeError("Fan tables do not hold 13 LOG and 39 GOG rows")
    accs = [a for a, u in U.items() if u["gene"] == "CREBBP"]
    out.append(dict(set="LOG_text", acc=accs[0] if len(accs) == 1 else "", gene="CREBBP", pos=1978, ref="N", alt="D", cancers="ovarian cancer (Fan text, HGMD)", n_fan=1, text_named=True))
    return out

def fan_status(U, r):
    u = U.get(r["acc"])
    if u is None: return dict(status="accession_not_in_UniProt_2026_03", found="", gene_uniprot="", lost=[], gained=[], literal_same=None)
    seq = u["seq"]; p = r["pos"]; found = seq[p - 1] if 1 <= p <= len(seq) else ""
    if found != r["ref"]: return dict(status="reference_residue_differs", found=found, gene_uniprot=u["gene"], lost=[], gained=[], literal_same=None)
    lost, gained = H6.classify_change(seq, p, r["ref"], r["alt"]); l2, g2 = H6.classify_literal(seq, p, r["ref"], r["alt"])
    return dict(status="ok", found=found, gene_uniprot=u["gene"], lost=lost, gained=gained, literal_same=(lost == l2 and gained == g2))

def scan(log, fan):
    """One pass over the chosen-file mutation rows (Q1 set) plus the Fan look-ups in all rows. Returns dict of structures."""
    fan_keys = {(r["acc"], r["pos"], r["ref"], r["alt"]) for r in fan}; fan_pos = {(r["acc"], r["pos"]) for r in fan}
    posc = collections.Counter(); vd = {}; cls = collections.Counter(); gst = collections.Counter(); guard = collections.Counter()
    lostseq = collections.defaultdict(set); gainseq = collections.defaultdict(set); prots = collections.defaultdict(set)
    fan_proj = collections.defaultdict(collections.Counter); fan_ch = collections.Counter(); fan_nh = collections.Counter()
    fan_all = collections.defaultdict(set); fan_gs = collections.defaultdict(set)
    pos_alts = collections.defaultdict(collections.Counter); n_chosen = 0
    with gzip.open(MUTS, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            acc = r["accession"]; pos = int(r["position"]); pk = (acc, pos); chosen = (r["in_chosen_file"] == "1")
            if pk in fan_pos:
                key = (acc, pos, r["ref"], r["alt"])
                if key in fan_keys: fan_all[key].add(r["case_id"])
                if chosen: pos_alts[pk][(r["ref"], r["alt"])] += 1
            if not chosen: continue
            n_chosen += 1; cls[r["class"]] += 1
            proj = r["project"]; hyp = (r["hypermutated"] == "1"); variant = (acc, pos, r["ref"], r["alt"])
            if variant in fan_keys:
                fan_proj[variant][proj] += 1; fan_ch[variant] += 1; fan_nh[variant] += (0 if hyp else 1)
                if r["gain_state_best"]: fan_gs[variant].add(r["gain_state_best"])
            sets = ["main"] if hyp else ["main", "nohyper"]
            cats = []; lt = r["lost_tiers"]; gs = r["gain_state_best"]
            if lt:
                ns = r["lost_sequon_N"].split(";"); ts = lt.split(";")
                if len(ns) != len(ts): guard["lost_sequon_N and lost_tiers have different numbers of entries"] += 1; continue
                tset = sorted(set(ts)); cats += ["loss_" + t for t in tset] + ["loss_any"]
                if "T1" in tset or "T2" in tset: cats.append("loss_T12")
                for n_, t_ in zip(ns, ts):
                    rel = pos - int(n_)
                    if rel not in (0, 1, 2): guard["lost sequon position outside 0 to 2 from the N"] += 1; continue
                    for s in sets: posc[(s, t_, rel)] += 1
                    lostseq["loss_any"].add((acc, int(n_)))
                    if t_ in ("T1", "T2"): lostseq["loss_T12"].add((acc, int(n_)))
                    if t_ == "T1": lostseq["loss_T1"].add((acc, int(n_)))
            if gs:
                cats += ["gain_" + gs, "gain_any"]; gst[gs] += 1
                if gs in ("plausible", "unknown"): cats.append("gain_pu")
                gn = r["gained_sequon_N"].split(";"); gstates = r["gain_states"].split(";")
                if len(gn) != len(gstates): gstates = [gs] * len(gn)
                for n_, s_ in zip(gn, gstates):
                    if n_ == "": continue
                    gainseq["gain_any"].add((acc, int(n_)))
                    if s_ in ("plausible", "unknown"): gainseq["gain_pu"].add((acc, int(n_)))
            for c in cats:
                e = vd.get((c, variant))
                if e is None: e = vd[(c, variant)] = [acc, 0, 0, collections.Counter(), r["gene"]]
                e[1] += 1; e[2] += (0 if hyp else 1); e[3][proj] += 1; prots[c].add(acc)
    log(f"chosen-file mutation rows read: {n_chosen} | classes {dict(cls)} | gain states {dict(gst)} | guards (must be 0): {dict(guard) if guard else 0}")
    if guard: raise RuntimeError("scan guard failed")
    return dict(posc=posc, vd=vd, cls=cls, gst=gst, lostseq=lostseq, gainseq=gainseq, prots=prots, fan_proj=fan_proj, fan_ch=fan_ch, fan_nh=fan_nh,
                fan_all=fan_all, fan_gs=fan_gs, pos_alts=pos_alts, n_chosen=n_chosen)

def section1(log):
    check_md5(log)
    U = H.load_uniprot_main(); fan = fan_rows(U); stat = [fan_status(U, r) for r in fan]
    tab = [(r, s) for r, s in zip(fan, stat) if not r["text_named"]]
    n_all = len(tab); n_ok = sum(s["status"] == "ok" for r, s in tab)
    lit_bad = sum(1 for r, s in tab if s["status"] == "ok" and not s["literal_same"])
    agree = 0; notes = []
    for r, s in tab:
        if s["status"] != "ok":
            notes.append(f"  NOT CHECKABLE {r['set']} {r['gene']} {r['ref']}{r['pos']}{r['alt']} ({r['acc']}): {s['status']} (UniProt residue '{s['found']}')"); continue
        good = bool(s["lost"]) if r["set"] == "LOG" else bool(s["gained"]); agree += good
        if not good: notes.append(f"  CLASS DIFFERS {r['set']} {r['gene']} {r['ref']}{r['pos']}{r['alt']} ({r['acc']}): our classifier lost {s['lost']} gained {s['gained']}")
        if s["gene_uniprot"] != r["gene"]: notes.append(f"  gene symbol differs {r['acc']}: Fan {r['gene']}, UniProt {s['gene_uniprot']}")
    log(f"Fan table variants: {n_all} | checkable (accession in UniProt 2026_03 and reference residue equal): {n_ok} | classifier agrees with Fan's category: {agree} of {n_ok} | classify_literal differs from classify_change: {lit_bad}")
    for x in notes: log(x)
    c0 = (n_ok >= 0.8 * n_all); c1 = (n_ok > 0 and agree >= 0.8 * n_ok and lit_bad == 0)
    log(f"C-F0 (>= 80% of the table variants checkable): {'PASSED' if c0 else 'FAILED'}")
    log(f"C-F1 (classifier agrees for >= 80% of checkable variants, literal version identical): {'PASSED' if c1 else 'FAILED'}")
    sc = scan(log, fan); diffs = []
    if dict(sc["cls"]) != EXPECT_CLASS: diffs.append(f"class counts {dict(sc['cls'])} differ from {EXPECT_CLASS}")
    if dict(sc["gst"]) != EXPECT_GAIN: diffs.append(f"gain-state counts {dict(sc['gst'])} differ from {EXPECT_GAIN}")
    relnum = {v: k for k, v in RELNAME.items()}; n_pos = 0
    for r in read_tsv(Q1_POS):
        if r["scope"] != "pooled": continue
        n_pos += 1
        if sc["posc"][(r["set"], r["tier"], relnum[r["sequon_position"]])] != int(r["n_sequon_losses"]): diffs.append("position count differs: " + str(r))
    if n_pos != 2 * 3 * 3: diffs.append(f"expected 18 pooled position lines, found {n_pos}")
    saved = read_tsv(Q1_VARS); n_var = 0
    for r in saved:
        n_var += 1; e = sc["vd"].get((r["category"], (r["accession"], int(r["position"]), r["ref"], r["alt"])))
        if e is None or (e[1], e[2], len(e[3])) != (int(r["n_cases"]), int(r["n_cases_without_hypermutated"]), int(r["n_cancer_types"])): diffs.append("variant line differs: " + r["category"] + " " + r["accession"] + " " + r["position"])
    for c in VARCATS:
        mine = sum(1 for (cc, v), e in sc["vd"].items() if cc == c and len(e[3]) >= 3); theirs = sum(1 for r in saved if r["category"] == c and int(r["n_cancer_types"]) >= 3)
        if mine != theirs: diffs.append(f"{c}: variants in >= 3 types {mine} here, {theirs} in the saved table")
    log(f"C-F2: position lines compared {n_pos}, variant lines compared {n_var}, differences {len(diffs)}")
    for d in diffs[:20]: log("  " + d)
    c2 = (not diffs)
    log(f"C-F2 (rebuilt tables equal the saved SqnPAN_09 Q1 tables): {'PASSED' if c2 else 'FAILED'}")
    with open(VERDICT, "w", encoding="utf-8") as o:
        o.write(f"P2_13 controls verdict {time.strftime('%Y-%m-%d %H:%M:%S')}\nC-F0: {'PASSED' if c0 else 'FAILED'}\nC-F1: {'PASSED' if c1 else 'FAILED'}\nC-F2: {'PASSED' if c2 else 'FAILED'}\n")
    log(f"verdict file {VERDICT} written")
    if not c2: raise RuntimeError("C-F2 FAILED: section 2 will not run")

def require_verdict():
    if not os.path.exists(VERDICT): raise RuntimeError("no verdict file: run section 1 first")
    with open(VERDICT, encoding="utf-8") as vf: txt = vf.read()
    if "C-F2: PASSED" not in txt: raise RuntimeError("C-F2 did not pass: stop")
    return ("C-F0: PASSED" in txt and "C-F1: PASSED" in txt)

def part1(log, fan, stat, sc, tier_of, sfx):
    rows = []; summ = collections.defaultdict(collections.Counter)
    for r, s in zip(fan, stat):
        key = (r["acc"], r["pos"], r["ref"], r["alt"]); pk = (r["acc"], r["pos"])
        n_ch = sc["fan_ch"][key]; n_nh = sc["fan_nh"][key]; n_any = len(sc["fan_all"][key]); projs = sc["fan_proj"][key]; npj = len(projs)
        alts = sorted(((al, n) for (rf, al), n in sc["pos_alts"][pk].items() if rf == r["ref"] and al != r["alt"]), key=lambda x: (-x[1], x[0]))
        if s["status"] != "ok": label = "Not checkable (" + s["status"] + ")"
        elif n_ch >= 1: label = "Reproduced" + sfx
        elif alts: label = "Same position, other substitution" + sfx
        else: label = "Not seen" + sfx
        sclass = ("loss_and_gain" if s["lost"] and s["gained"] else "loss" if s["lost"] else "gain" if s["gained"] else "none") if s["status"] == "ok" else ""
        agrees = ("" if s["status"] != "ok" else "yes" if ((r["set"] != "GOG" and bool(s["lost"])) or (r["set"] == "GOG" and bool(s["gained"]))) else "no")
        tiers = ";".join(tier_of.get((r["acc"], n), "?") for n in s["lost"]); recur = "yes" if npj >= 3 else "no"
        rows.append([r["set"], r["gene"], s["gene_uniprot"], r["acc"], r["pos"], r["ref"], r["alt"], r["n_fan"], r["cancers"], s["status"], s["found"], sclass, agrees,
                     ";".join(map(str, s["lost"])), tiers, ";".join(map(str, s["gained"])), ";".join(sorted(sc["fan_gs"][key])), n_ch, n_nh, n_any, npj,
                     ";".join(f"{p}:{n}" for p, n in sorted(projs.items(), key=lambda x: (-x[1], x[0]))), ";".join(f"{al}:{n}" for al, n in alts), label, recur])
        summ[r["set"]][label.split(" (")[0] if label.startswith("Not checkable") else label] += 1
        if recur == "yes": summ[r["set"]]["recurrent in TCGA (>= 3 projects)"] += 1
        log(f"  {r['set']:8s} {r['gene']:9s} {r['ref']}{r['pos']}{r['alt']} ({r['acc']}): {label} | exact cases {n_ch} (no hypermutated {n_nh}, any primary file {n_any}), TCGA projects {npj}"
            f" | class {sclass or '-'} tier {tiers or '-'} | other substitutions here: {';'.join(f'{al}:{n}' for al, n in alts) or '-'}")
    hdr = ["fan_set", "gene_fan", "gene_uniprot", "accession", "position", "ref", "alt", "n_fan_cancer_types", "fan_cancer_types", "status", "residue_in_uniprot", "sequence_class",
           "class_agrees_with_fan", "lost_sequon_N", "lost_sequon_tier", "gained_sequon_N", "gain_state_observed", "cases_exact_chosen_file", "cases_exact_without_hypermutated",
           "cases_exact_any_primary_file", "n_tcga_projects_exact", "tcga_projects_exact", "other_substitutions_at_residue", "label", "recurrent_in_tcga"]
    write_tsv(OUT["fan_variants"], hdr, rows)
    log(f"table {OUT['fan_variants']} | md5 {H.md5_file(OUT['fan_variants'])}")
    log("=== Part 1 summary (label counts) ===")
    for k in ("LOG", "GOG", "LOG_text"): log(f"  {k}: {dict(summ[k])}")
    lt = collections.Counter(x[14].split(";")[0] if x[14] else "none" for x in rows if x[0] == "LOG" and x[9] == "ok")
    log(f"  tier of the sequon lost by Fan's 13 LOG variants in UniProt 2026_03 (first sequon): {dict(lt)}")

def part2(log, sc, sfx):
    q2 = {(r["tier"], r["group"]): r for r in read_tsv(Q2_RES) if r["config"] == "main" and r["family"] in ("primary", "pooled_other")}
    if len(q2) != 9: raise RuntimeError(f"expected 9 pooled main units in the Q2 table, found {len(q2)}")
    rows = []; log("=== Part 2: position of the lost sequon hit (loss mutations; chosen files) and Q2 observed/expected (main set, pooled) ===")
    for t in ("T1", "T2", "T3"):
        pm = [sc["posc"][("main", t, k)] for k in (0, 1, 2)]; pn = [sc["posc"][("nohyper", t, k)] for k in (0, 1, 2)]
        qn = q2[(t, "N")]; qs = q2[(t, "ST")]
        on, en, os_, es = float(qn["n_events_glyco"]), float(qn["expected_glyco"]), float(qs["n_events_glyco"]), float(qs["expected_glyco"])
        import math
        ratio_obs = os_ / on if on > 0 else float("nan"); ratio_exp = es / en if en > 0 else float("nan")
        lor = math.log2(ratio_obs / ratio_exp) if on > 0 and en > 0 and es > 0 and os_ > 0 else float("nan")
        changed = any(q["final_label"] in ("Enriched", "Depleted") for q in (qn, qs))
        expected = "yes" if (lor == lor and abs(lor) < 0.32 and not changed) else "no"
        raw = "yes" if pm[2] > pm[0] else "no"
        rows.append([t, *pm, round(pm[2] / pm[0], 3) if pm[0] else "", *pn, round(pn[2] / pn[0], 3) if pn[0] else "", raw, on, en, qn["final_label"], os_, es, qs["final_label"],
                     round(ratio_obs, 3), round(ratio_exp, 3), round(lor, 3), expected + sfx])
        log(f"  {t}: positions N / X to P / S/T = {pm[0]} / {pm[1]} / {pm[2]} (ratio S/T to N {pm[2] / pm[0]:.2f}); without hypermutated {pn[0]} / {pn[1]} / {pn[2]}")
        log(f"      Q2 eligible strata: N observed {on:.0f} expected {en:.1f} ({qn['final_label']}); S/T observed {os_:.0f} expected {es:.1f} ({qs['final_label']}); "
            f"observed ratio {ratio_obs:.2f}, expected ratio {ratio_exp:.2f}, log2 ratio of ratios {lor:.3f}")
        log(f"      raw position-3 excess (Fan statement) reproduced: {raw}; expected from background rates (rule of the pre-specification): {expected}")
    hdr = ["tier", "main_pos1_N", "main_pos2_XtoP", "main_pos3_ST", "main_ratio_3_to_1", "nohyper_pos1_N", "nohyper_pos2_XtoP", "nohyper_pos3_ST", "nohyper_ratio_3_to_1", "raw_position3_excess_reproduced",
           "q2_N_observed", "q2_N_expected", "q2_N_final_label", "q2_ST_observed", "q2_ST_expected", "q2_ST_final_label", "observed_ratio_ST_to_N", "expected_ratio_ST_to_N", "log2_ratio_of_ratios", "expected_from_background"]
    write_tsv(OUT["fan_positions"], hdr, rows); log(f"table {OUT['fan_positions']} | md5 {H.md5_file(OUT['fan_positions'])}")
    log("  Q2 counts are events (case x position) in strata with glycosite and background positions; Q1 position counts are mutation rows in all proteins, so the two ratios are not identical.")

def part3(log, fan, sc):
    vd = sc["vd"]; nvar = collections.Counter(c for (c, v) in vd)
    rec = {c: {v for (cc, v), e in vd.items() if cc == c and len(e[3]) >= 3} for c in ("loss_T12", "loss_any", "gain_pu", "gain_any", "loss_T1")}
    ours = {("loss_any", "variants"): nvar["loss_any"], ("loss_any", "sequons"): len(sc["lostseq"]["loss_any"]), ("loss_any", "proteins"): len(sc["prots"]["loss_any"]),
            ("gain_any", "variants"): nvar["gain_any"], ("gain_any", "sequons"): len(sc["gainseq"]["gain_any"]), ("gain_any", "proteins"): len(sc["prots"]["gain_any"])}
    rows = []
    for name, fv, cat, unit in FAN_COUNTS: rows.append([name, fv, ours[(cat, unit)], cat + " " + unit, "different data, sources and sequon definitions; ours is TCGA only"])
    rows.append(["High-confidence sequons abolished by somatic variants", 2898, len(sc["lostseq"]["loss_T12"]), "loss_T12 sequons (tier T1 or T2)", "Fan high-confidence = UniProt FT with evidence or manual assertion, HPRD, dbPTM, neXtProt, CDD"])
    rows.append(["  same, tier T1 only", "", len(sc["lostseq"]["loss_T1"]), "loss_T1 sequons", "T1 = experimental or combined evidence, or atlas"])
    rows.append(["High-confidence LOG variants in >= 3 cancer types", 13, len(rec["loss_T12"]), "loss_T12 variants in >= 3 TCGA projects", "Fan: DOID cancer types from many databases; ours: TCGA projects"])
    rows.append(["LOG variants in >= 3 cancers (signal peptide or Secreted/Membrane proteins, any NLG)", 24, len(rec["loss_any"]), "loss_any variants in >= 3 TCGA projects (no localisation filter)", "not like for like (no filter here)"])
    rows.append(["GOG variants in >= 3 cancers (signal peptide or Secreted/Membrane proteins)", 41, len(rec["gain_pu"]), "gain_pu variants in >= 3 TCGA projects", "Fan text says 41, printed Table 5 has 39 rows; our flag is Amendment 2 (plausible or unknown)"])
    rows.append(["  same, all gains", "", len(rec["gain_any"]), "gain_any variants in >= 3 TCGA projects", ""])
    rows.append(["Somatic nsSNVs collected", 1272878, sc["n_chosen"], "chosen-file missense rows (rows, not unique variants)", "Fan: five databases; ours: one TCGA file per case"])
    write_tsv(OUT["fan_counts"], ["quantity", "fan_2018", "ours", "our_definition", "note"], rows)
    log("=== Part 3: counts side by side (Fan 2018 | ours) ===")
    for r in rows: log(f"  {r[0]}: Fan {r[1] if r[1] != '' else '-'} | ours {r[2]} ({r[3]})")
    fk4 = {(r["acc"], r["pos"], r["ref"], r["alt"]) for r in fan if r["set"] == "LOG"}; fk5 = {(r["acc"], r["pos"], r["ref"], r["alt"]) for r in fan if r["set"] == "GOG"}
    orows = []
    for c in ("loss_T12", "loss_any", "gain_pu", "gain_any"):
        for v in sorted(rec[c], key=lambda x: (-len(vd[(c, x)][3]), -vd[(c, x)][1], x)):
            e = vd[(c, v)]
            orows.append([c, v[0], e[4], v[1], v[2], v[3], e[1], e[2], len(e[3]), ";".join(f"{p}:{n}" for p, n in sorted(e[3].items(), key=lambda x: (-x[1], x[0]))), "yes" if v in fk4 else "no", "yes" if v in fk5 else "no"])
    write_tsv(OUT["fan_recurrent_overlap"], ["our_category", "accession", "gene", "position", "ref", "alt", "n_cases", "n_cases_without_hypermutated", "n_tcga_projects", "projects", "in_fan_table4_LOG", "in_fan_table5_GOG"], orows)
    log(f"table {OUT['fan_recurrent_overlap']} ({len(orows)} rows) | md5 {H.md5_file(OUT['fan_recurrent_overlap'])}")
    for c, fk, nm in (("loss_T12", fk4, "Fan Table 4 (13 LOG)"), ("loss_any", fk4, "Fan Table 4 (13 LOG)"), ("gain_pu", fk5, "Fan Table 5 (39 GOG)"), ("gain_any", fk5, "Fan Table 5 (39 GOG)")):
        ov = sorted(rec[c] & fk); log(f"  our {c} variants in >= 3 TCGA projects: {len(rec[c])}; in {nm}: {len(ov)} " + ("(" + ", ".join(f"{vd[(c, v)][4]} {v[2]}{v[1]}{v[3]}" for v in ov) + ")" if ov else ""))

def part4(log, U, sc):
    vd = sc["vd"]; per = collections.defaultdict(collections.Counter); perseq = collections.defaultdict(collections.Counter)
    for (c, v) in vd: per[c][v[0]] += 1
    for c in ("loss_any", "loss_T12", "loss_T1", "gain_any", "gain_pu"):
        src = sc["lostseq"] if c.startswith("loss") else sc["gainseq"]
        for (a, n) in src[c]: perseq[c][a] += 1
    def rank(counter, acc):
        n = counter.get(acc, 0)
        return (1 + sum(1 for x in counter.values() if x > n)) if n > 0 else ""
    rows = []; log("=== Part 4: genes named in the Fan text (ours: distinct variants and sequons; rank 1 = most among all proteins) ===")
    for g, stmt in FAN_GENES:
        for a in sorted(a for a, u in U.items() if u["gene"] == g):
            ln = U[a]["length"]; nl = perseq["loss_any"].get(a, 0)
            rows.append([g, a, ln, stmt, per["loss_any"].get(a, 0), perseq["loss_any"].get(a, 0), perseq["loss_T12"].get(a, 0), perseq["loss_T1"].get(a, 0), rank(per["loss_any"], a),
                         per["gain_any"].get(a, 0), per["gain_pu"].get(a, 0), perseq["gain_any"].get(a, 0), rank(per["gain_any"], a), round(ln / nl, 1) if nl else ""])
            log(f"  {g} ({a}, {ln} aa): Fan: {stmt} | ours: loss variants {rows[-1][4]} (rank {rows[-1][8]}), lost sequons {rows[-1][5]} (T1 or T2 {rows[-1][6]}, T1 {rows[-1][7]}), gain variants {rows[-1][9]} (plausible or unknown {rows[-1][10]}, rank {rows[-1][12]})")
    write_tsv(OUT["fan_genes"], ["gene", "accession", "length", "fan_statement", "loss_variants", "lost_sequons_any_tier", "lost_sequons_T1_or_T2", "lost_sequons_T1", "rank_loss_variants",
                                 "gain_variants", "gain_variants_plausible_or_unknown", "gained_sequons", "rank_gain_variants", "residues_per_lost_sequon"], rows)
    log(f"table {OUT['fan_genes']} | md5 {H.md5_file(OUT['fan_genes'])}")
    for c in ("loss_any", "gain_any"):
        top = per[c].most_common(5); log(f"  top 5 proteins by distinct {c} variants: " + ", ".join(f"{U[a]['gene']} {n}" for a, n in top))

def section2(log):
    ok_tab = require_verdict(); sfx = "" if ok_tab else "_table_check_FAILED"
    log(f"Fan table checks C-F0 and C-F1 both passed: {'yes' if ok_tab else 'NO (labels carry _table_check_FAILED)'}")
    U = H.load_uniprot_main(); fan = fan_rows(U); stat = [fan_status(U, r) for r in fan]
    tier_of = {}
    with open(SITES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"): tier_of[(r["accession"], int(r["N_position"]))] = r["tier"]
    sc = scan(log, fan)
    log("=== Part 1: Fan variants looked up in the chosen-file mutation table (labelled cross-check, not blind) ===")
    part1(log, fan, stat, sc, tier_of, sfx); part2(log, sc, sfx); part3(log, fan, sc); part4(log, U, sc)
    log("Reminder: Fan et al. 2018 used other data and databases; cancer types are not comparable (DOID terms against TCGA projects); the position-3 excess and the top variants are not blind; "
        "a variant not seen in TCGA alone says nothing about Fan's other databases; Table 5 printed 39 rows against 41 in the text.")

SECTIONS = {"1": ("checks_controls", section1), "2": ("crosscheck_results", section2)}

def run_one(k):
    name, fn = SECTIONS[k]; ok = {"v": False}
    def wrapped(log):
        fn(log); ok["v"] = True
    H.run_section("P2_13", k, name, wrapped)
    return ok["v"]

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for k in (list(SECTIONS) if which == "all" else [which]):
        if not run_one(k): print(f"SECTION {k} DID NOT FINISH: the chain stops here (see the log)", flush=True); break
# ---- END OF SCRIPT ----
