# SqnPAN_05_helpers.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Reusable functions for SqnPAN scripts SqnPAN_05 onward: reuse them, do not rewrite them.
# - paths and Log class (same log pattern as SqnPAN_04)
# - read_xlsx: reads one sheet of an .xlsx with the standard library only (no openpyxl or pandas needed). Cells are placed by their
#   column letter, so empty cells never shift the columns. Every value is returned as text exactly as stored in the file.
# - load_uniprot_main / load_uniprot_features: the two UniProt reviewed human 2026_03 tables (stream.gz and the features table)
# - parse_glyco, parse_features, note_class, pos_kind, sequon_status: parsers for the UniProt feature text and for sequons
import os, re, sys, csv, gzip, time, hashlib, zipfile, traceback, collections
import xml.etree.ElementTree as ET

PAPER2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(PAPER2, "data", "raw"); PROC = os.path.join(PAPER2, "data", "processed")
LOGS = os.environ.get("P2_LOGDIR", os.path.join(PAPER2, "results", "logs"))      # P2_LOGDIR is used only for test runs
TABLES = os.environ.get("P2_TABLEDIR", os.path.join(PAPER2, "results", "tables"))  # P2_TABLEDIR is used only for test runs
PROC_OUT = os.environ.get("P2_PROCDIR", PROC)                                     # P2_PROCDIR is used only for test runs
UNIPROT_MAIN = os.path.join(RAW, "Uniprot", "stream.gz")
UNIPROT_FEAT = os.path.join(RAW, "Uniprot", "uniprot_2026_03_features_stream.gz")
ATLAS = os.path.join(RAW, "nglycositeatlas", "HumanAll.xlsx")
CDR = os.path.join(RAW, "TCGA-CDR-SupplementalTableS1.xlsx")
BAILEY = os.path.join(RAW, "Bailey2018_TableS1_driver_genes.xlsx")
MANIFEST = os.path.join(RAW, "TCGA_masked_maf_manifest_2026-10-21.txt")
CASE_PROJECT = os.path.join(RAW, "GDC_TCGA_masked_maf_case_project_2026-10-22.tsv")
MAF_DIR = os.path.join(RAW, "maf_tcga")
USE_TSV = os.path.join(PROC, "P2_04_alphafold_use.tsv")
csv.field_size_limit(2 ** 31 - 1)

class Log:
    def __init__(self, prefix, k, name):
        os.makedirs(LOGS, exist_ok=True)
        self.path = os.path.join(LOGS, f"{prefix}_section{k}_{name}.txt"); self.k = k; self.prefix = prefix
        self.f = open(self.path, "a", encoding="utf-8")
        self("=" * 70); self(f"{prefix} section {k} ({name}) started {time.strftime('%Y-%m-%d %H:%M:%S')}")
    def __call__(self, *a):
        s = " ".join(str(x) for x in a); print(s, flush=True); self.f.write(s + "\n"); self.f.flush()
    def close(self):
        self(f"---- END OF SECTION {self.k} ---- {time.strftime('%Y-%m-%d %H:%M:%S')}"); self.f.close()

def md5_file(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""): h.update(b)
    return h.hexdigest()

# ---------------- xlsx reader (standard library only) ----------------
_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

def _col_index(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0): n = n * 26 + ord(ch) - 64
    return n - 1

def xlsx_sheet_names(path):
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
    return [s.get("name") for s in wb.iter(_NS + "sheet")]

def read_xlsx(path, sheet):
    """All rows of one sheet as a list of lists of text. List index = row number - 1 (blank rows are empty lists)."""
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(_NS + "si"):
                shared.append("".join(t.text or "" for t in si.iter(_NS + "t")))
        wb = ET.fromstring(z.read("xl/workbook.xml")); rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rid = {r.get("Id"): r.get("Target") for r in rels}; target = None
        for s in wb.iter(_NS + "sheet"):
            if s.get("name") == sheet: target = rid[s.get(_RNS + "id")]
        if target is None: raise KeyError(f"sheet {sheet!r} not found in {path}")
        member = target.lstrip("/") if target.startswith("/") else "xl/" + target
        rows = {}
        with z.open(member) as f:
            for _ev, el in ET.iterparse(f):
                if el.tag != _NS + "row": continue
                cells = {}
                for c in el.iter(_NS + "c"):
                    t = c.get("t"); v = c.find(_NS + "v")
                    if t == "s": val = shared[int(v.text)] if v is not None and v.text is not None else ""
                    elif t == "inlineStr": val = "".join(x.text or "" for x in c.iter(_NS + "t"))
                    else: val = v.text if v is not None and v.text is not None else ""
                    cells[_col_index(c.get("r"))] = val
                rn = int(el.get("r")); width = (max(cells) + 1) if cells else 0
                rows[rn] = [cells.get(i, "") for i in range(width)]
                el.clear()
    last = max(rows) if rows else 0
    return [rows.get(i, []) for i in range(1, last + 1)]

# ---------------- UniProt tables ----------------
def load_uniprot_main(path=UNIPROT_MAIN):
    """{entry: dict(name, gene, length, seq, glyco, ensembl, af)} from stream.gz (columns checked by the caller)."""
    out = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[r["Entry"]] = {"name": r["Entry Name"], "gene": r["Gene Names (primary)"], "length": int(r["Length"]), "seq": r["Sequence"],
                               "glyco": r["Glycosylation"], "ensembl": r["Ensembl"], "af": r["AlphaFoldDB"]}
    return out

def load_uniprot_features(path=UNIPROT_FEAT):
    """{entry: row dict} from the second table (feature columns as text)."""
    out = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"): out[r["Entry"]] = r
    return out

_REC = re.compile(r'CARBOHYD\s+([^;]+?);(.*?)(?=(?:;\s*)?CARBOHYD\s|\Z)', re.S)
def parse_glyco(text):
    """List of dict(pos_raw, note, codes) for every CARBOHYD record. codes = sorted set of ECO codes in the evidence field."""
    out = []
    for m in _REC.finditer(text or ""):
        rest = m.group(2); note = re.search(r'/note="([^"]*)"', rest); ev = re.search(r'/evidence="([^"]*)"', rest)
        out.append({"pos_raw": m.group(1).strip(), "note": note.group(1) if note else "",
                    "codes": sorted(set(re.findall(r'ECO:\d{7}', ev.group(1)))) if ev else []})
    return out

def note_class(note):
    if note.startswith("N-linked"): return "N-linked asparagine" if "asparagine" in note else "N-linked other residue"
    m = re.match(r"^(?:\([^)]*\)\s*)?([A-Za-z-]+)-linked", note)
    return (m.group(1) + "-linked") if m else "other"

def pos_kind(pos_raw):
    """('single', int) | ('isoform', None) | ('range', None) | ('other', None)"""
    if pos_raw.isdigit(): return "single", int(pos_raw)
    if re.match(r"^[A-Z0-9]+-\d+:\d+$", pos_raw): return "isoform", None
    if ".." in pos_raw: return "range", None
    return "other", None

_FSTART = re.compile(r'(?:^|(?<=; ))(?=[A-Z_]{3,10}\s+[<>?\d])')
_FHEAD = re.compile(r'^([A-Z_]{3,10})\s+([<>?]?\d+|\?)(?:\.\.([<>?]?\d+|\?))?(.*)$', re.S)
def parse_features(text):
    """List of dict(key, start, end, uncertain, note). start/end are ints or None (unknown); a single position has end == start."""
    text = text or ""; starts = [m.start() for m in _FSTART.finditer(text)]; out = []
    for n, st in enumerate(starts):
        rec = text[st:(starts[n + 1] if n + 1 < len(starts) else len(text))]
        m = _FHEAD.match(rec)
        if not m: continue
        a, b = m.group(2), m.group(3)
        def num(t):
            d = re.sub(r"\D", "", t) if t else ""
            return int(d) if d else None
        s = num(a); e = num(b) if b is not None else s
        note = re.search(r'/note="([^"]*)"', m.group(4))
        out.append({"key": m.group(1), "start": s, "end": e, "uncertain": any(c in (a or "") + (b or "") for c in "<>?"), "note": note.group(1) if note else ""})
    return out

def sequon_status(seq, i):
    """Status of the residue at 1-based position i as the N of a sequon."""
    if i < 1 or i > len(seq): return "position_outside_sequence"
    if seq[i - 1] != "N": return "residue_not_N"
    if i + 2 > len(seq): return "N_within_2_of_C_terminus"
    x, z = seq[i], seq[i + 1]
    if z in "ST": return "N-X-S/T" if x != "P" else "N-P-S/T"
    if z == "C" and x != "P": return "N-X-C"
    return "N_other_context"

def sequon_positions(seq):
    """1-based positions of N in N-X-S/T (X not P), the N-X-C positions and the N-P-S/T positions."""
    return ([m.start() + 1 for m in re.finditer(r"(?=N[^P][ST])", seq)], [m.start() + 1 for m in re.finditer(r"(?=N[^P]C)", seq)],
            [m.start() + 1 for m in re.finditer(r"(?=NP[ST])", seq)])

def run_section(prefix, k, name, fn):
    log = Log(prefix, k, name); t0 = time.time()
    try:
        fn(log)
    except Exception:
        log("ERROR:"); log(traceback.format_exc())
    finally:
        log(f"seconds: {round(time.time() - t0)}"); log.close()
# ---- END OF HELPERS ----
