# SqnPAN_06_helpers.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Reusable functions for SqnPAN scripts SqnPAN_06 onward: reuse them, do not rewrite them. Builds on SqnPAN_05_helpers.py (not modified).
# - classify_change: the loss/gain classifier (set difference of the sequons around the changed residue)
# - classify_literal: an independent version of analysis plan section 5, written rule by rule (used by control C1 as the cross-check)
# - build_tiers: tiers T1/T2/T3 (analysis plan decision 2A, atlas rule 3A); same logic as SqnPAN_05 section 3
# - load_topology, topology_category, plausibility: SqnPAN_05 decision 10 categories and the three-state flag of analysis plan Amendment 2
# - build_revmap: UniProt Ensembl field -> {ENST: [(accession, isoform_ok)]}, same logic as SqnPAN_05 section 4
import os, re, sys, csv, gzip, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import SqnPAN_05_helpers as H

AA = "ACDEFGHIKLMNPQRSTVWY"
T1CODES = {"ECO:0000269", "ECO:0007744"}
SELECTION = os.path.join(H.PAPER2, "results", "tables", "P2_05_case_file_selection.tsv")   # always read from the real folder
FILESUM = os.path.join(H.PROC, "P2_05_maf_file_summary.tsv")

# ---------------- classifier ----------------
def _is_sequon(get, n, p):
    """N-X-S/T with X not P; p = 1-based position of the N; get(i) = residue at 1-based i."""
    return p >= 1 and p + 2 <= n and get(p) == "N" and get(p + 1) != "P" and get(p + 2) in ("S", "T")

def classify_change(seq, pos, ref, alt):
    """(lost, gained): sorted lists of 1-based N positions of sequons lost / gained by the substitution ref->alt at pos.
    Set difference of the sequons that contain pos (N at pos-2, pos-1 or pos), before and after. S<->T changes nothing."""
    if not (1 <= pos <= len(seq)) or seq[pos - 1] != ref: raise ValueError(f"reference residue {ref} not at position {pos}")
    n = len(seq); base = lambda i: seq[i - 1]; mut = lambda i: alt if i == pos else seq[i - 1]
    cand = (pos - 2, pos - 1, pos)
    before = [p for p in cand if _is_sequon(base, n, p)]; after = [p for p in cand if _is_sequon(mut, n, p)]
    return sorted(p for p in before if p not in after), sorted(p for p in after if p not in before)

def classify_literal(seq, pos, ref, alt):
    """Independent version, one block per rule of analysis plan section 5 (loss: N changed; X changed to P; S/T changed to non-S/T.
    Gain: new N; new S/T at +2; P at position 2 changed to non-P)."""
    n = len(seq)
    def r(i): return seq[i - 1] if 1 <= i <= n else ""
    lost = []; gained = []
    if ref == "N" and alt != "N" and r(pos + 1) not in ("P", "") and r(pos + 2) in ("S", "T"): lost.append(pos)                       # loss 1
    if alt == "P" and ref != "P" and r(pos - 1) == "N" and r(pos + 1) in ("S", "T"): lost.append(pos - 1)                          # loss 2
    if ref in ("S", "T") and alt not in ("S", "T") and r(pos - 2) == "N" and r(pos - 1) not in ("P", ""): lost.append(pos - 2)     # loss 3
    if alt == "N" and ref != "N" and r(pos + 1) not in ("P", "") and r(pos + 2) in ("S", "T"): gained.append(pos)                   # gain 1
    if alt in ("S", "T") and ref not in ("S", "T") and r(pos - 2) == "N" and r(pos - 1) not in ("P", ""): gained.append(pos - 2)    # gain 2
    if ref == "P" and alt != "P" and r(pos - 1) == "N" and r(pos + 1) in ("S", "T"): gained.append(pos - 1)                        # gain 3
    return sorted(lost), sorted(gained)

# ---------------- tiers ----------------
def build_tiers(U):
    """Same logic as SqnPAN_05 section 3. Returns dict(tier_of {(acc, N position): 'T1'|'T2'|'T3'}, codes, annotated, atlas, counts)."""
    sites = {}
    for acc, u in U.items():
        for r in H.parse_glyco(u["glyco"]):
            if not r["note"].startswith("N-linked"): continue
            kind, pos = H.pos_kind(r["pos_raw"])
            if kind != "single": continue
            if (acc, pos) in sites: sites[(acc, pos)]["codes"] |= set(r["codes"]); continue
            sites[(acc, pos)] = {"codes": set(r["codes"]), "status": H.sequon_status(u["seq"], pos)}
    rows = H.read_xlsx(H.ATLAS, "Glycosite and glycoprotein DB")
    hi = next(i for i, r in enumerate(rows) if r and r[0] == "Database"); hdr = rows[hi]
    body = [r + [""] * (len(hdr) - len(r)) for r in rows[hi + 1:] if r and r[0] != ""]; ix = {h: i for i, h in enumerate(hdr)}
    uniq = set()
    for r in body:
        if r[ix["Database"]] != "Reviewed": continue
        acc = r[ix["Accession (UniProtKB)"]]
        try: pos = int(float(r[ix["Glycosylation location"]]))
        except ValueError: continue
        if acc in U: uniq.add((acc, pos))
    atlas_valid = {k for k in uniq if H.sequon_status(U[k[0]]["seq"], k[1]) == "N-X-S/T"}
    T1u = {k for k, s in sites.items() if s["status"] == "N-X-S/T" and s["codes"] & T1CODES}
    T2u = {k for k, s in sites.items() if s["status"] == "N-X-S/T" and not (s["codes"] & T1CODES)}
    atypical = {k for k, s in sites.items() if s["status"] != "N-X-S/T"}
    allseq = set()
    for acc, u in U.items(): allseq |= {(acc, i) for i in H.sequon_positions(u["seq"])[0]}
    T1 = T1u | atlas_valid; T2 = T2u - T1; T3 = allseq - T1 - T2
    tier_of = {}
    for t, S in (("T1", T1), ("T2", T2), ("T3", T3)):
        for k in S: tier_of[k] = t
    return {"tier_of": tier_of, "codes": {k: s["codes"] for k, s in sites.items()}, "annotated": set(sites), "atlas": atlas_valid,
            "counts": {"T1": len(T1), "T2": len(T2), "T3": len(T3), "all_sequons": len(allseq), "atlas_reviewed_sites": len(uniq),
                       "atlas_valid": len(atlas_valid), "uniprot_annotated_not_sequon": len(atypical), "tiers_disjoint_and_complete": len(T1) + len(T2) + len(T3) == len(allseq) == len(tier_of)}}

# ---------------- topology and gain plausibility ----------------
def load_topology(F):
    topo = {}
    for acc, r in F.items():
        def feats(col): return [f for f in H.parse_features(r[col]) if f["start"] is not None and f["end"] is not None]
        topo[acc] = {"topo": [(f["start"], f["end"], f["note"]) for f in feats("Topological domain")], "tm": [(f["start"], f["end"]) for f in feats("Transmembrane")],
                     "im": [(f["start"], f["end"]) for f in feats("Intramembrane")], "sig": [(f["start"], f["end"]) for f in feats("Signal peptide")]}
    return topo

def is_cyto(note): return note.startswith("Cytoplasmic")      # Amendment 2: non-cytoplasmic = a TOPO_DOM whose note does not start with "Cytoplasmic"

def _flags(d, p):
    cyto = any(a <= p <= b and is_cyto(n) for a, b, n in d["topo"]); noncyto = any(a <= p <= b and not is_cyto(n) for a, b, n in d["topo"])
    tm = any(a <= p <= b for a, b in d["tm"]); im = any(a <= p <= b for a, b in d["im"])
    return cyto, noncyto, tm, im

def topology_category(d, p):
    """Same categories and order as SqnPAN_05 section 3 (decision 10 table)."""
    cyto, noncyto, tm, im = _flags(d, p)
    if noncyto: return "TOPO_DOM non-cytoplasmic"
    if cyto: return "TOPO_DOM cytoplasmic"
    if tm: return "TRANSMEM"
    if im: return "INTRAMEM"
    if any(a <= p <= b for a, b in d["sig"]): return "SIGNAL"
    if d["topo"] or d["tm"]: return "no region at this position (protein has topology)"
    return "protein without topology annotation"

def plausibility(d, p, has_t1):
    """Analysis plan Amendment 2. Plausible = protein has a T1 site AND (N in a non-cytoplasmic TOPO_DOM, or protein has a signal peptide and N is not in a
    cytoplasmic TOPO_DOM). Implausible = not plausible and N in a cytoplasmic TOPO_DOM, TRANSMEM or INTRAMEM. Unknown = all other."""
    cyto, noncyto, tm, im = _flags(d, p)
    if has_t1 and (noncyto or (bool(d["sig"]) and not cyto)): return "plausible"
    if cyto or tm or im: return "implausible"
    return "unknown"

# ---------------- transcripts ----------------
def build_revmap(U):
    """{ENST: [(accession, isoform_ok)]}; isoform_ok = no isoform tag or tag -1 (rule 5A). Same logic as SqnPAN_05 section 4."""
    revmap = collections.defaultdict(list)
    for acc, u in U.items():
        for m in re.finditer(r"(ENST\d+)(?:\.\d+)?\s*(?:\[([^\]]*)\])?", u["ensembl"]):
            br = m.group(2); tag = None
            if br is not None: tag = re.match(r"^(.+?)(?:-(\d+))?$", br).group(2)
            revmap[m.group(1)].append((acc, tag in (None, "1")))
    return revmap
# ---- END OF HELPERS ----
