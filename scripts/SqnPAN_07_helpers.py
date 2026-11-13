# SqnPAN_07_helpers.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Reusable functions for per-residue pLDDT and relative solvent exposure of the AlphaFold v6 human models: reuse them, do not rewrite them.
# Shared logic with script 06 of the companion LUAD study. NumPy only (no scipy, no Biopython, no freesasa).
# Method (pre-specification fixed 2026-10-24, header of SqnPAN_07_structure_features.py):
#   Shrake-Rupley accessible surface, probe 1.4 A, 401 Fibonacci points per atom (DSSP's point formula, 2N+1 with N = 200), DSSP atom radii,
#   relative exposure = residue area / Tien et al. 2013 theoretical maximum (Table 1, ALLOWED region).
# Constants and where they come from:
#   - DSSP radii: N 1.65, CA 1.87, C 1.76, O 1.40, any other (side-chain) atom 1.80, probe 1.40, surface sphere with 2*200+1 = 401 points. READ FROM THE SOURCE of
#     PDB-REDO/dssp (trunk, libdssp/src/dssp.cpp) on 2026-10-24 (constants kRadiusN, kRadiusCA, kRadiusC, kRadiusO, kRadiusSideAtom, kRadiusWater, kN = 200).
#     Tien et al. used the 2011 DSSP version, which was NOT checked. The Fibonacci point formula below is WRITTEN FROM MEMORY of DSSP's: any equal-area spiral is
#     valid, and controls C-A and C-E of SqnPAN_07 test the result.
#   - Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. Maximum allowed solvent accessibilites of residues in proteins. PLoS One 2013;8(11):e80635,
#     doi 10.1371/journal.pone.0080635, PMID 24278298. Values: Table 1, column Theoretical (ALLOWED region), copied from the PMC3836772 full text on 2026-10-24.
import os, sys, re, csv, gzip, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import SqnPAN_04_alphafold_helpers as AF
import SqnPAN_05_helpers as H

PROBE = 1.4
N_POINTS = 401
RADIUS_BY_ATOM = {"N": 1.65, "CA": 1.87, "C": 1.76, "O": 1.40}
RADIUS_SIDE = 1.80
TIEN_MAX = {"A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
            "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0}
THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L",
                "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "UNK": "X"}
assert len(TIEN_MAX) == 20 and sum(TIEN_MAX.values()) == 4040.0, "Tien table was edited"
ATOM_BLOCK = 1024     # atoms per neighbour-search block
SR_BLOCK = 128        # atoms per matrix block in the point test
TABLE_HEADER = ["accession", "position", "aa", "plddt", "asa", "rsa", "fragment", "n_fragments", "terminal", "cut_distance"]

def unit_points(n):
    """n = 2N+1 points on the unit sphere, DSSP's Fibonacci spiral (latitude asin(2i/n), longitude fmod(i, golden ratio) * 2 pi / golden ratio, i = -N..N)."""
    assert n % 2 == 1 and n >= 3
    N = (n - 1) // 2
    i = np.arange(-N, N + 1, dtype=float)
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    z = 2.0 * i / n
    lon = np.fmod(i, phi) * 2.0 * math.pi / phi
    rxy = np.sqrt(1.0 - z * z)
    return np.stack([rxy * np.cos(lon), rxy * np.sin(lon), z], axis=1)

# ---------------- neighbour search (cell grid, NumPy only) ----------------
def _grid(xyz, R):
    cs = 2.0 * float(R.max()) + 1e-6                      # cell size >= the largest overlap distance R_i + R_j
    cell = np.floor((xyz - xyz.min(axis=0)) / cs).astype(np.int64) + 1
    dims = cell.max(axis=0) + 2                           # one empty cell of padding on each side: no wrap-around in the keys
    key = (cell[:, 0] * dims[1] + cell[:, 1]) * dims[2] + cell[:, 2]
    order = np.argsort(key, kind="stable")
    offs = np.array([(a * dims[1] + b) * dims[2] + c for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)], dtype=np.int64)
    return key, order, key[order], offs

def _block_pairs(grid, xyz, R, ids):
    """Neighbour pairs (i, j) for the atoms in ids: spheres of radius R overlap (d < R_i + R_j), j != i. Returns i_local, j, q = xyz_j - xyz_i, d2."""
    key, order, skey, offs = grid
    il, jj = [], []
    for off in offs:
        nk = key[ids] + off
        lo = np.searchsorted(skey, nk, side="left"); hi = np.searchsorted(skey, nk, side="right")
        cnt = hi - lo; tot = int(cnt.sum())
        if tot == 0: continue
        il.append(np.repeat(np.arange(len(ids)), cnt))
        within = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        jj.append(order[np.repeat(lo, cnt) + within])
    if not il:
        return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros((0, 3)), np.zeros(0)
    il = np.concatenate(il); jj = np.concatenate(jj)
    gi = ids[il]
    q = xyz[jj] - xyz[gi]; d2 = (q * q).sum(axis=1)
    keep = (jj != gi) & (d2 < (R[gi] + R[jj]) ** 2)
    return il[keep], jj[keep], q[keep], d2[keep]

def neighbour_pair_set(xyz, R):
    """Set of (i, j) neighbour pairs found by the grid (used by control C-B only)."""
    grid = _grid(xyz, R); out = set()
    for b0 in range(0, len(xyz), ATOM_BLOCK):
        ids = np.arange(b0, min(b0 + ATOM_BLOCK, len(xyz)))
        il, jj, _, _ = _block_pairs(grid, xyz, R, ids)
        out.update(zip(ids[il].tolist(), jj.tolist()))
    return out

def exposed_points(xyz, R, U):
    """Number of exposed surface points for every atom. R = atom radius + probe. A point of atom i is exposed when it lies outside the expanded sphere of
    every other atom. With p = x_i + R_i u and q_j = x_j - x_i, 'inside sphere j' is |p - x_j|^2 < R_j^2, i.e. 2 R_i (u . q_j) > R_i^2 + |q_j|^2 - R_j^2."""
    n = len(xyz); npts = len(U)
    out = np.full(n, npts, dtype=np.int64)
    if n < 2: return out
    grid = _grid(xyz, R)
    for b0 in range(0, n, ATOM_BLOCK):
        ids = np.arange(b0, min(b0 + ATOM_BLOCK, n)); m = len(ids)
        il, jj, q, d2 = _block_pairs(grid, xyz, R, ids)
        if len(il) == 0: continue
        gi = ids[il]
        T = (R[gi] ** 2 + d2 - R[jj] ** 2) / (2.0 * R[gi])
        o = np.argsort(il, kind="stable"); il = il[o]; q = q[o]; T = T[o]
        cnt = np.bincount(il, minlength=m); rank = np.arange(len(il)) - (np.cumsum(cnt) - cnt)[il]
        K = int(cnt.max())
        Q = np.zeros((m, K, 3)); TT = np.full((m, K), np.inf)
        Q[il, rank] = q; TT[il, rank] = T
        for s in range(0, m, SR_BLOCK):
            e = min(s + SR_BLOCK, m)
            S = np.matmul(U, np.ascontiguousarray(np.transpose(Q[s:e], (0, 2, 1))))     # (atoms, points, neighbours)
            blocked = (S > TT[s:e][:, None, :]).any(axis=2)
            out[ids[s:e]] = npts - blocked.sum(axis=1)
    return out

def exposed_points_bruteforce(xyz, R, U):
    """Independent version for control C-B: direct distances against every other atom, no grid and no algebra."""
    n = len(xyz); out = np.zeros(n, dtype=np.int64); idx = np.arange(n)
    for i in range(n):
        P = xyz[i] + R[i] * U
        d = np.sqrt(((xyz - xyz[i]) ** 2).sum(axis=1))
        nb = np.where((d < R + R[i]) & (idx != i))[0]
        blocked = np.zeros(len(U), dtype=bool)
        for j in nb: blocked |= ((P - xyz[j]) ** 2).sum(axis=1) < R[j] ** 2
        out[i] = len(U) - int(blocked.sum())
    return out

def atom_areas(R, exposed, npts):
    return 4.0 * math.pi * R * R * exposed / npts

# ---------------- reading one AlphaFold model ----------------
def parse_atoms(lines):
    """Heavy atoms of one AlphaFold mmCIF. Returns coordinates, radii, residue ordinal per atom, residue names and pLDDT (B-factor of the CA atom) per residue."""
    cols = []; ix = None
    xyz, rad, seq = [], [], []; comp = {}; plddt = {}; skipped = {"H": 0, "HETATM": 0}
    for l in lines:
        if l.startswith("_atom_site."):
            cols.append(l.split()[0]); continue
        if l.startswith("ATOM") or l.startswith("HETATM"):
            if ix is None:
                ix = {k: cols.index("_atom_site." + k) for k in ("group_PDB", "type_symbol", "label_atom_id", "label_comp_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv")}
            v = l.split()
            if v[ix["group_PDB"]] != "ATOM": skipped["HETATM"] += 1; continue
            if v[ix["type_symbol"]] == "H": skipped["H"] += 1; continue
            name = v[ix["label_atom_id"]].strip('"'); s = int(v[ix["label_seq_id"]])
            xyz.append((float(v[ix["Cartn_x"]]), float(v[ix["Cartn_y"]]), float(v[ix["Cartn_z"]])))
            rad.append(RADIUS_BY_ATOM.get(name, RADIUS_SIDE)); seq.append(s)
            if s not in comp: comp[s] = v[ix["label_comp_id"]]
            if name == "CA": plddt[s] = float(v[ix["B_iso_or_equiv"]])
    seq = np.array(seq); s0 = int(seq.min()); n_res = int(seq.max()) - s0 + 1
    if len(comp) != n_res or len(plddt) != n_res: raise ValueError(f"residues not contiguous or CA missing: {len(comp)} residues, {len(plddt)} CA, span {n_res}")
    order = list(range(s0, s0 + n_res))
    return {"xyz": np.array(xyz), "radii": np.array(rad), "res": seq - s0, "n_res": n_res, "comp": [comp[s] for s in order],
            "plddt": np.array([plddt[s] for s in order]), "skipped": skipped}

def model_features(text, U):
    """Per-residue features of one model: UniProt range, one-letter sequence, pLDDT, residue area (A^2) and the number of atoms."""
    lines = text.split("\n")
    begin = int(AF._cif_value(lines, "_ma_target_ref_db_details.seq_db_align_begin")); end = int(AF._cif_value(lines, "_ma_target_ref_db_details.seq_db_align_end"))
    a = parse_atoms(lines)
    if begin + a["n_res"] - 1 != end: raise ValueError(f"UniProt range {begin}-{end} does not match {a['n_res']} residues")
    R = a["radii"] + PROBE
    cnt = exposed_points(a["xyz"], R, U)
    asa = np.bincount(a["res"], weights=atom_areas(R, cnt, len(U)), minlength=a["n_res"])
    return {"begin": begin, "end": end, "n_res": a["n_res"], "aa": "".join(THREE_TO_ONE.get(c, "X") for c in a["comp"]), "plddt": a["plddt"], "asa": asa,
            "n_atoms": len(R), "skipped": a["skipped"]}

# ---------------- context: archive, index, ranges, sequences ----------------
def load_use():
    """{accession: (gene, length, n_fragments)} for the proteins with status 'use' (rule 1A)."""
    out = {}
    with open(H.USE_TSV, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["status"] == "use": out[r["accession"]] = (r["gene"], int(r["uniprot_length"]), int(r["n_fragments"]))
    return out

class Ctx:
    """Open archive + index + fragment ranges + sequences for the usable proteins. One per process."""
    def __init__(self, n_points=N_POINTS):
        self.use = load_use(); self.U = unit_points(n_points)
        self.fh = AF.open_archive(); self.idx = AF.load_index()
        self.ranges = {}
        with open(os.path.join(H.PROC, "P2_04_alphafold_models.tsv"), newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["accession"] in self.use: self.ranges.setdefault(r["accession"], []).append((int(r["fragment"]), int(r["align_begin"]), int(r["align_end"])))
        self.seqs = {}
        with open(os.path.join(H.PROC, "P2_04_alphafold_sequences.tsv"), newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["accession"] in self.use: self.seqs[r["accession"]] = r["sequence"]
    def close(self): self.fh.close()

def protein_rows(acc, ctx, U=None):
    """All rows of one protein (one per UniProt position). Every fragment model is computed; each position takes the fragment chosen by rule 2A.
    Returns (list of tab-separated lines, stats)."""
    U = ctx.U if U is None else U
    gene, L, nfrag = ctx.use[acc]; seq = ctx.seqs[acc]; ranges = sorted(ctx.ranges[acc])
    if len(seq) != L or len(ranges) != nfrag: raise ValueError(f"{acc}: length or fragment count differs from the SqnPAN_04 tables")
    feats = {fr: model_features(AF.read_model(ctx.fh, ctx.idx, acc, fr, "cif"), U) for fr, b, e in ranges}
    lines = []; mism = 0; blank = 0; gt1 = 0
    for pos in range(1, L + 1):
        fr = 1 if nfrag == 1 else AF.choose_fragment(ranges, pos, L)
        m = feats[fr]; i = pos - m["begin"]
        aa = m["aa"][i]
        if aa != seq[pos - 1]: mism += 1
        mx = TIEN_MAX.get(aa)
        rsa = "" if mx is None else f"{m['asa'][i] / mx:.3f}"
        if mx is None: blank += 1
        elif m["asa"][i] / mx > 1.0: gt1 += 1
        term = 1 if (pos == 1 or pos == L) else 0
        if nfrag == 1: cut = ""
        else:
            d = min(pos - m["begin"] if m["begin"] > 1 else 10 ** 9, m["end"] - pos if m["end"] < L else 10 ** 9)
            cut = "" if d == 10 ** 9 else str(d)
        lines.append(f"{acc}\t{pos}\t{aa}\t{m['plddt'][i]:.2f}\t{m['asa'][i]:.1f}\t{rsa}\t{fr}\t{nfrag}\t{term}\t{cut}")
    return lines, {"rows": L, "mismatch": mism, "rsa_blank": blank, "rsa_gt1": gt1}

# ---------------- workers (module level so that Windows "spawn" can import them) ----------------
_W = {}
def worker_init():
    _W["ctx"] = Ctx()

def process_batch(task):
    """task = (batch id, [accessions], part directory). Writes part_<id>.tsv.gz and part_<id>.stats (the stats file is written last: a batch is complete only with both)."""
    bid, accs, partdir = task
    ctx = _W["ctx"]; t0 = time.time()
    path = os.path.join(partdir, f"part_{bid:05d}.tsv.gz"); tmp = path + ".tmp"; stats = os.path.join(partdir, f"part_{bid:05d}.stats")
    tot = {"rows": 0, "mismatch": 0, "rsa_blank": 0, "rsa_gt1": 0}
    with gzip.open(tmp, "wt", encoding="utf-8", newline="\n", compresslevel=6) as o:
        for acc in accs:
            try:
                lines, st = protein_rows(acc, ctx)
            except Exception as e:
                raise RuntimeError(f"{acc}: {type(e).__name__}: {e}")
            o.write("\n".join(lines)); o.write("\n")
            for k in tot: tot[k] += st[k]
    os.replace(tmp, path)
    secs = round(time.time() - t0, 1)
    with open(stats, "w") as o: o.write("\t".join(str(x) for x in (bid, len(accs), tot["rows"], tot["mismatch"], tot["rsa_blank"], tot["rsa_gt1"], secs)) + "\n")
    return bid, len(accs), tot, secs
# ---- END OF HELPERS ----
