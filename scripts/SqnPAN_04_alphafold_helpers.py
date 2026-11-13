# SqnPAN_04_alphafold_helpers.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Reusable functions for the AlphaFold human archive (UP000005640_9606_HUMAN_v6.tar). Later scripts import these: reuse them, do not rewrite them.
# The tar is never unpacked (rule 3A): models are read by seeking to their position in the tar, using the index
# data/processed/P2_04_alphafold_index.tsv made by SqnPAN_04_check_alphafold_archive.py section 1.
# Fragment rule (rule 2A): a residue covered by several fragments is taken from the fragment in which it lies
# furthest from an artificial cut end (the protein's own N and C termini are not cuts); ties go to the lower fragment number.
import os, re, gzip, csv

PAPER2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAR_PATH = os.path.join(PAPER2, "data", "raw", "alphafold", "UP000005640_9606_HUMAN_v6.tar")
INDEX_PATH = os.path.join(PAPER2, "data", "processed", "P2_04_alphafold_index.tsv")
NAME_RE = re.compile(r"^AF-([A-Z0-9]+(?:-\d+)?)-F(\d+)-model_v(\d+)\.(cif|pdb)\.gz$")

def parse_name(name):
    """'AF-P12345-F2-model_v6.cif.gz' -> (accession, fragment, version, kind) or None."""
    m = NAME_RE.match(os.path.basename(name))
    return None if m is None else (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))

def load_index(path=INDEX_PATH):
    """Returns {(accession, fragment, kind): (offset_data, size)}."""
    idx = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            idx[(r["accession"], int(r["fragment"]), r["kind"])] = (int(r["offset_data"]), int(r["size"]))
    return idx

def open_archive(path=TAR_PATH):
    return open(path, "rb")

def read_member(fh, offset, size):
    """Read one gzip member from the open tar at its data offset and return the decompressed text (gzip CRC is checked)."""
    fh.seek(offset)
    raw = fh.read(size)
    if len(raw) != size:
        raise IOError(f"short read at offset {offset}: {len(raw)} of {size} bytes")
    return gzip.decompress(raw).decode()

def read_model(fh, idx, accession, fragment=1, kind="cif"):
    off, size = idx[(accession, fragment, kind)]
    return read_member(fh, off, size)

def fragments_of(idx, accession, kind="cif"):
    return sorted(f for (a, f, k) in idx if a == accession and k == kind)

def _cif_value(lines, key):
    """Value of a non-loop mmCIF item; handles values on the same line and ;-delimited text blocks."""
    for i, l in enumerate(lines):
        if l.startswith(key + " ") or l == key:
            rest = l[len(key):].strip()
            if rest:
                return rest.strip('"')
            if i + 1 < len(lines) and lines[i + 1].startswith(";"):
                out = [lines[i + 1][1:]]
                for l2 in lines[i + 2:]:
                    if l2.startswith(";"):
                        break
                    out.append(l2)
                return "".join(x.strip() for x in out)
    raise KeyError(key)

def parse_cif(text):
    """Structure facts from one AlphaFold mmCIF: UniProt range, sequence, per-residue CA rows with pLDDT (B-factor column)."""
    lines = text.split("\n")
    d = {
        "accession": _cif_value(lines, "_ma_target_ref_db_details.db_accession"),
        "gene": _cif_value(lines, "_ma_target_ref_db_details.gene_name"),
        "align_begin": int(_cif_value(lines, "_ma_target_ref_db_details.seq_db_align_begin")),
        "align_end": int(_cif_value(lines, "_ma_target_ref_db_details.seq_db_align_end")),
        "seq_version_date": _cif_value(lines, "_ma_target_ref_db_details.seq_db_sequence_version_date"),
        "seq_checksum": _cif_value(lines, "_ma_target_ref_db_details.seq_db_sequence_checksum"),
        "sequence": _cif_value(lines, "_entity_poly.pdbx_seq_one_letter_code_can").replace(" ", ""),
    }
    cols, ca = [], []
    in_atoms = False
    for l in lines:
        if l.startswith("_atom_site."):
            cols.append(l.split()[0][len("_atom_site."):]); in_atoms = True; continue
        if in_atoms:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                v = l.split()
                if v[cols.index("label_atom_id")] == "CA":
                    ca.append((int(v[cols.index("label_seq_id")]), v[cols.index("label_comp_id")],
                               int(v[cols.index("pdbx_sifts_xref_db_num")]), v[cols.index("pdbx_sifts_xref_db_res")],
                               float(v[cols.index("B_iso_or_equiv")])))
            elif cols and (l.startswith("#") or l.startswith("loop_") or l.startswith("_")):
                break
    d["ca"] = ca  # list of (label_seq_id, residue name, UniProt number, UniProt residue, pLDDT)
    return d

def pdb_ca_count(text):
    return sum(1 for l in text.split("\n") if l.startswith("ATOM") and l[12:16].strip() == "CA")

def choose_fragment(frag_ranges, pos, protein_length):
    """frag_ranges: list of (fragment, begin, end) in UniProt numbering. Returns the fragment for UniProt position pos (rule 2A), or None."""
    best = None
    for f, b, e in frag_ranges:
        if not (b <= pos <= e):
            continue
        left = pos - b if b > 1 else float("inf")
        right = e - pos if e < protein_length else float("inf")
        dist = min(left, right)
        if best is None or dist > best[0] or (dist == best[0] and f < best[1]):
            best = (dist, f)
    return None if best is None else best[1]
# ---- END OF FILE ----
