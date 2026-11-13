# SqnPAN_14_context_helpers.py - shared code for A1 (context-aware Q2) and A2 (germline positive control), Amendment 4
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Shared code for SqnPAN_14_context_q2.py (A1) and SqnPAN_15_germline_q2.py (A2): reuse, do not rewrite. Importing it does nothing.
#  - genetic code, the 96 pyrimidine-centred SBS channels (COSMIC names, for example A[C>T]G), GENCODE v36 reader, translation
#  - transcript per UniProt protein: a GENCODE v36 transcript whose CDS translation equals the UniProt 2026_03 sequence exactly
#  - the SqnPAN_08 position list rebuilt with position numbers; codon opportunities per position (single-nucleotide changes that make a counted substitution)
#  - MAF reader for pool workers; mapping of a missense row to a UniProt position (same logic as SqnPAN_06 section 3, copied)
#  - Engine: the Q2 test of SqnPAN_08 / SqnPAN_09 (Mantel-Haenszel rate ratio over protein x residue strata, bootstrap over proteins, labels, BH families) with the
#    exposure of a position = its expected rate (mode ctx) or 1 (mode uniform = the SqnPAN_08 test); C2/C3 repeats, real runs, capability, final labels.
#    In mode uniform on the full SqnPAN_08 position list the Engine reproduces SqnPAN_08 line for line (control C-S of SqnPAN_14).
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, re, csv, gzip, math, time, collections, hashlib
import multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_06_helpers as P
import SqnPAN_08_selection_controls as S8

GENCODE = os.path.join(H.RAW, "gencode", "gencode.v36.pc_transcripts.fa.gz"); GENCODE_MD5 = "61eaea20b1275e2e19def3e5882ed7ae"
REAL_TABLES = os.path.join(H.PAPER2, "results", "tables")
NPZ8 = os.path.join(H.PROC, "P2_08_inputs.npz"); NPZ8_MD5 = "c663941bae9e31dae5f68590a7624a4d"
RUN1_P208 = os.path.join(REAL_TABLES, "P2_08_C2C3_run1.tsv"); RUN1_P208_MD5 = "c872da79c028e3010c46d255b7ff540f"
RES_P209 = os.path.join(REAL_TABLES, "P2_09_Q2_results.tsv"); RES_P209_MD5 = "2e25c9b14f3cc0984341648af07f8698"
NOT_TESTED = "Not tested (too few data)"; CAN = "Test can answer"; NOTSENS = "Inconclusive, test not sensitive"; FLAG_LOW = "low-sensitivity unit"
LEVEL1 = [("1.0", 1.0)]
GROUPS = S8.GROUPS; TIERS = S8.TIERS; fmt = S8.fmt
GORDER = {"N": 0, "ST": 1, "both": 2}
COLS_REAL = ["run", "config", "tier", "group", "scope", "family", "tested", "n_informative", "expected_glyco", "eligible_proteins", "num", "den", "log2rr", "lo95", "hi95",
             "lo90", "hi90", "p", "fdr", "n_inf", "n_dropped", "raw_label", "n_events_glyco", "n_events_bg", "n_positions_glyco", "n_positions_bg"]

# ---------------- genetic code and channels ----------------
_B = "TCAG"; _AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AAS[16 * i + 4 * j + k] for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
_COMP = str.maketrans("ACGT", "TGCA")
def comp(s): return s.translate(_COMP)
def revcomp(s): return s.translate(_COMP)[::-1]
CHANNELS = sorted(f"{a}[{r}>{x}]{b}" for r in "CT" for x in "ACGT" if x != r for a in "ACGT" for b in "ACGT")
assert len(CHANNELS) == 96
CH_IDX = {c: i for i, c in enumerate(CHANNELS)}
CH_TRI = [c[0] + c[2] + c[6] for c in CHANNELS]
CH_SUB = [c[2:5] for c in CHANNELS]
TRI32 = sorted(set(CH_TRI)); TRI_IDX = {t: i for i, t in enumerate(TRI32)}
CPG_CT = [i for i, c in enumerate(CHANNELS) if c[2:5] == "C>T" and c[6] == "G"]

def channel(five, ref, alt, three):
    """Channel index of a single-nucleotide change (any strand), or -1 if a base is not A/C/G/T or alt equals ref."""
    if alt == ref or any(x not in "ACGT" for x in (five, ref, alt, three)): return -1
    if ref in "AG": five, ref, alt, three = comp(three), comp(ref), comp(alt), comp(five)
    return CH_IDX[f"{five}[{ref}>{alt}]{three}"]

_CHT = {}
def chan_tx(tri, alt):
    k = (tri, alt); v = _CHT.get(k)
    if v is None:
        v = channel(tri[0], tri[1], alt, tri[2]) if len(tri) == 3 else -1; _CHT[k] = v
    return v

def tri_pyr(tri):
    if len(tri) != 3 or any(x not in "ACGT" for x in tri): return None
    return tri if tri[1] in "CT" else revcomp(tri)

def _chg(codon, cls):
    aa = CODON[codon]; out = []
    if aa == "*": return out
    for k in range(3):
        for alt in "ACGT":
            if alt == codon[k]: continue
            aa2 = CODON[codon[:k] + alt + codon[k + 1:]]
            if aa2 == "*" or aa2 == aa: continue
            if cls == "ST" and aa2 in "ST": continue
            out.append((k, alt, aa2))
    return out
# CHG[cls][codon] = single-nucleotide changes that make a counted substitution: N = any missense; ST = missense to a residue that is not S or T; any = any missense
CHG = {cls: {c: _chg(c, cls) for c in CODON} for cls in ("N", "ST", "any")}

def translate(cds):
    if len(cds) % 3: return None
    p = "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds), 3))
    if p.endswith("*"): p = p[:-1]
    return None if ("*" in p or "X" in p) else p

# ---------------- GENCODE and transcript map ----------------
def load_gencode(path=GENCODE):
    """({ENST without version: [dict(ver, gene, seq, cs, ce)]}, records, records without one CDS field). cs, ce = 0-based CDS start, end (exclusive).
    Versions per ID are ordered: not PAR_Y first, then higher version number first."""
    out = collections.defaultdict(list); cnt = {"n": 0, "bad": 0}
    def add(h, parts):
        cnt["n"] += 1; f = h.split("|"); tid, _, ver = f[0].partition(".")
        cds = [x for x in f if x.startswith("CDS:")]
        if len(cds) != 1: cnt["bad"] += 1; return
        a, b = cds[0][4:].split("-")
        out[tid].append(dict(ver=ver, gene=f[1].split(".")[0], seq="".join(parts).upper(), cs=int(a) - 1, ce=int(b)))
    h = None; parts = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        for line in fh:
            if line.startswith(">"):
                if h is not None: add(h, parts)
                h = line[1:].strip(); parts = []
            else: parts.append(line.strip())
        if h is not None: add(h, parts)
    def vkey(d):
        m = re.match(r"\d+", d["ver"]); return ("PAR" in d["ver"], -(int(m.group()) if m else 0))
    for v in out.values(): v.sort(key=vkey)
    return dict(out), cnt["n"], cnt["bad"]

def build_transcript_map(U, G):
    """{accession: (ENST, version index, number of candidates, listed by UniProt)}: the GENCODE transcript whose CDS translation equals the UniProt sequence;
    preference: listed in the UniProt Ensembl field without isoform tag or with tag -1, then lowest ENST number, then version order."""
    revmap = P.build_revmap(U); listed = collections.defaultdict(set)
    for enst, lst in revmap.items():
        for acc, ok in lst:
            if ok: listed[acc].add(enst)
    by_seq = collections.defaultdict(list)
    for tid, vl in G.items():
        for vi, d in enumerate(vl):
            p = translate(d["seq"][d["cs"]:d["ce"]])
            if p: by_seq[p].append((tid, vi))
    M = {}; c = collections.Counter()
    for acc in sorted(U):
        cand = by_seq.get(U[acc]["seq"])
        if not cand: c["proteins with no GENCODE v36 transcript translating to the UniProt sequence (left out of A1)"] += 1; continue
        cand = sorted(cand, key=lambda x: (x[0] not in listed[acc], int(x[0][4:]) if x[0][4:].isdigit() else 10 ** 12, x[1]))
        tid, vi = cand[0]; M[acc] = (tid, vi, len(cand), tid in listed[acc])
        c["proteins mapped, chosen transcript listed by UniProt" if tid in listed[acc] else "proteins mapped, chosen transcript not listed by UniProt"] += 1
        if len(cand) > 1: c["proteins mapped with more than one candidate transcript (rule applied)"] += 1
    return M, c

def cds_trinucleotides(G, M):
    """Counts of the 32 pyrimidine-centred trinucleotides over every CDS base (with both flanks present) of the mapped transcripts (each transcript once)."""
    cnt = collections.Counter()
    for tid, vi in sorted({(v[0], v[1]) for v in M.values()}):
        d = G[tid][vi]; T = d["seq"]; lo = max(d["cs"], 1); hi = min(d["ce"], len(T) - 1)
        cnt.update(T[i - 1:i + 2] for i in range(lo, hi))
    out = np.zeros(32)
    for tri, n in cnt.items():
        k = tri_pyr(tri)
        if k is not None: out[TRI_IDX[k]] += n
    return out

# ---------------- positions (SqnPAN_08 order) and opportunities ----------------
def rebuild_positions(U):
    """The SqnPAN_08 position list in the same order (copy of SqnPAN_08 section 1), with position numbers."""
    accs = sorted(U); by_acc = collections.defaultdict(dict)
    with open(S8.SITES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"): by_acc[r["accession"]][int(r["N_position"])] = int(r["tier"][1])
    P_prot, P_res, P_role, P_posn = [], [], [], []
    for i, a in enumerate(accs):
        s = U[a]["seq"]; sq = by_acc.get(a, {}); role = {}
        for p, t in sq.items(): role[p] = t; role[p + 2] = t
        xpos = {p + 1 for p in sq}
        for m in re.finditer("[NST]", s):
            j = m.start() + 1; res = 0 if m.group() == "N" else 1
            if j in role: code = role[j]
            elif j in xpos: continue
            else: code = 0
            P_prot.append(i); P_res.append(res); P_role.append(code); P_posn.append(j)
    return accs, np.array(P_prot, np.int32), np.array(P_res, np.int8), np.array(P_role, np.int8), np.array(P_posn, np.int32)

def opportunities(accs, P_prot, P_res, P_posn, U, G, M):
    """OPP (positions x 9, channel index or -1) and P_ctx (position has at least one usable change). Flanks come from the transcript sequence."""
    npos = len(P_prot); OPP = np.full((npos, 9), -1, np.int16); ctx = np.zeros(npos, bool); c = collections.Counter()
    pl = P_prot.tolist(); rl = P_res.tolist(); ql = P_posn.tolist()
    for k in range(npos):
        acc = accs[pl[k]]; mp_ = M.get(acc)
        if mp_ is None: c["positions in proteins without a mapped transcript"] += 1; continue
        d = G[mp_[0]][mp_[1]]; T = d["seq"]; i = d["cs"] + 3 * (ql[k] - 1); cod = T[i:i + 3]
        if CODON.get(cod) != U[acc]["seq"][ql[k] - 1]: c["codon does not code the UniProt residue (must be 0)"] += 1; continue
        j = 0
        for kk, alt, _a in CHG["N" if rl[k] == 0 else "ST"][cod]:
            idx = i + kk
            if idx < 1 or idx + 1 >= len(T): c["changes without both flanking bases (left out)"] += 1; continue
            ch = chan_tx(T[idx - 1:idx + 2], alt)
            if ch < 0: c["changes with a non-ACGT flank (left out)"] += 1; continue
            OPP[k, j] = ch; j += 1
        if j: ctx[k] = True; c["positions with codon opportunities"] += 1
        else: c["positions without a usable change"] += 1
    return OPP, ctx, c

# ---------------- MAF reading (pool worker) and SqnPAN_06 mapping ----------------
MAF_NEED = ["Hugo_Symbol", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", "Transcript_ID", "Protein_position", "Amino_acids",
            "Codons", "CDS_position", "CONTEXT", "SWISSPROT", "all_effects"]
MAF_KEEP = {"Missense_Mutation": "M", "Silent": "S", "Nonsense_Mutation": "X"}

def maf_worker(job):
    """Reads one MAF file. Returns tallies (CONTEXT checks over every SNV row), the #filedate line and the coding rows (missense, silent, nonsense):
    (class, Variant_Type, Transcript_ID, Protein_position, Amino_acids, Codons, CDS_position, CONTEXT, Reference_Allele, Tumor_Seq_Allele2, SWISSPROT, Hugo_Symbol,
    all_effects (missense rows of a chosen file only, otherwise empty))."""
    fid, path, chosen = job
    t = collections.Counter(); rows = []; hdr = None; ix = None; fdate = "(none)"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        for line in f:
            if line[0] == "#":
                if line.startswith("#filedate"): fdate = line.strip().split(" ", 1)[-1]
                continue
            if hdr is None:
                hdr = line.rstrip("\n").split("\t"); miss = [k for k in MAF_NEED if k not in hdr]
                if miss: return dict(fid=fid, error="missing columns: " + ", ".join(miss))
                ix = [hdr.index(k) for k in MAF_NEED]; continue
            v = line.rstrip("\n").split("\t"); t["data rows"] += 1
            vt = v[ix[2]]
            if vt == "SNP":
                t["SNV rows"] += 1; ctx = v[ix[10]]
                if len(ctx) == 11:
                    t["SNV rows with an 11-letter CONTEXT"] += 1
                    if ctx[5].upper() == v[ix[3]]: t["SNV rows whose CONTEXT centre base equals Reference_Allele"] += 1
            c = MAF_KEEP.get(v[ix[1]])
            if c is None: continue
            rows.append((c, vt, v[ix[5]], v[ix[6]], v[ix[7]], v[ix[8]], v[ix[9]], v[ix[10]], v[ix[3]], v[ix[4]], v[ix[11]], v[ix[0]],
                         v[ix[12]] if (chosen and c == "M") else ""))
    return dict(fid=fid, fdate=fdate, tally=t, rows=rows)

def map_missense(sw, ppos, aas, hugo, ae, U, revmap):
    """Same logic as SqnPAN_06 section 3 (rule 5A and Amendment 1): returns ((accession, position, ref, alt), 'direct' or 'rescued') or (None, reason)."""
    acc = sw.split(".")[0].split("&")[0]; reason = None; ref = alt = None
    if not acc: reason = "no SWISSPROT accession"
    elif acc not in U: reason = "accession not in UniProt reviewed 2026_03"
    else:
        try: pos = int(ppos.split("/")[0]); tot = int(ppos.split("/")[1]); ref, alt = aas.split("/")
        except Exception: reason = "position or amino acids unparsable"
        if reason is None:
            sq = U[acc]["seq"]; lenok = tot == len(sq); resok = 1 <= pos <= len(sq) and len(ref) == 1 and sq[pos - 1] == ref
            reason = None if (lenok and resok) else "length or reference residue differs"
    if reason is None:
        if len(alt) != 1 or alt not in P.AA or alt == ref: return None, "skipped"
        return (acc, pos, ref, alt), "direct"
    targets = set()
    for e in ae.split(";"):
        fl = e.split(",")
        if len(fl) < 4 or "missense_variant" not in fl[1]: continue
        m = re.match(r"^p\.([A-Z])(\d+)([A-Z])$", fl[2])
        if not m: continue
        rp = int(m.group(2))
        for acc2, okflag in revmap.get(fl[3], ()):
            s2 = U[acc2]["seq"]
            if okflag and 1 <= rp <= len(s2) and s2[rp - 1] == m.group(1): targets.add((acc2, rp, m.group(1), m.group(3)))
    if acc:
        own = {t for t in targets if t[0] == acc}
        if len(own) != 1: return None, "dropped"
        key = next(iter(own))
    else:
        if len(targets) != 1: return None, "dropped"
        key = next(iter(targets))
        if U[key[0]]["gene"] != hugo: return None, "dropped"
    if key[3] == key[2] or key[3] not in P.AA: return None, "skipped"
    return key, "rescued"

# ---------------- Engine ----------------
def load_npz(path):
    with np.load(path, allow_pickle=False) as z: return {k: z[k] for k in z.files}

EVMASK = {"main": lambda D: D["ev_main"], "nohyper": lambda D: D["ev_main"] & ~D["ev_hyper"], "union": lambda D: D["ev_union"],
          "norescue": lambda D: D["ev_norescue"], "germA": lambda D: D["ev_germA"], "germB": lambda D: D["ev_germB"]}

class Engine:
    """configs: list of dict(name, kind 'main'|'strict', mask (EVMASK key), mode 'ctx'|'uniform', spec (R_<spec> in D) or None, cancer (per cancer type scopes at T1)).
    restricted: positions without codon opportunities (P_ctx False) are left out of both classes."""
    def __init__(self, D, configs, restricted, tiers=(1, 2, 3), role=None):
        self.cfgs = [dict(c) for c in configs]; self.tiers = tuple(tiers)
        prot = D["P_prot"].astype(np.int64); res = D["P_res"].astype(np.int64)
        role = (D["P_role"] if role is None else role).astype(np.int8)
        if restricted: role = np.where(D["P_ctx"], role, -1).astype(np.int8)
        usable = D["P_usable"]; pb = D["P_pb"].astype(np.int64); rb = D["P_rb"].astype(np.int64)
        npos = len(prot); nprot = len(D["accs"]); self.projects = [str(x) for x in D["projects"]]; ncan = len(self.projects); self.ncan = ncan
        ev_pos = D["ev_pos"].astype(np.int64); ev_can = D["ev_cancer"].astype(np.int64); self.nev = len(ev_pos)
        self.R = {}; self.wpos = {}
        for c in self.cfgs:
            if c["mode"] == "ctx" and c["spec"] not in self.R:
                R = D["R_" + c["spec"]].astype(np.float64); self.R[c["spec"]] = R; O = D["OPP"].astype(np.int64)
                self.wpos[c["spec"]] = np.where(O >= 0, R[np.maximum(O, 0), 0], 0.0).sum(axis=1)
        need_ent = ncan > 0 and any(c["mode"] == "ctx" and c["cancer"] for c in self.cfgs)
        kinds = [k for k in ("main", "strict") if any(c["kind"] == k for c in self.cfgs)]
        self.K = {}
        for kind in kinds:
            valid = np.ones(npos, bool) if kind == "main" else (usable & (pb >= 0) & (rb >= 0))
            strat = prot * 2 + res if kind == "main" else ((prot * 2 + res) * 2 + np.maximum(pb, 0)) * 2 + np.maximum(rb, 0)
            for t in self.tiers:
                pool = np.flatnonzero(valid & ((role == t) | (role == 0)))
                uniq, st = np.unique(strat[pool], return_inverse=True); st = st.ravel()
                isg = role[pool] == t
                n1 = np.bincount(st, weights=isg.astype(np.float64), minlength=len(uniq)).round().astype(np.int64)
                n0 = np.bincount(st, minlength=len(uniq)).astype(np.int64) - n1
                if kind == "main": s_prot = uniq // 2; s_res = uniq % 2
                else: s_prot = uniq // 8; s_res = (uniq // 4) % 2
                loc = np.full(npos, -1, np.int64); loc[pool] = np.arange(len(pool))
                Kt = dict(pool=pool, st=st, n1=n1, n0=n0, tot=(n1 + n0).astype(np.float64), elig=(n1 > 0) & (n0 > 0), s_prot=s_prot, s_res=s_res, loc=loc,
                          npool=len(pool), isg_real=isg)
                if need_ent and kind == "main" and t == 1:
                    O = D["OPP"][pool]; ok = O >= 0
                    Kt["ent_p"] = np.nonzero(ok)[0].astype(np.int32); Kt["ent_ch"] = O[ok].astype(np.int8)
                self.K[(kind, t)] = Kt
        self.S = {}
        for c in self.cfgs:
            kind = c["kind"]; m_ev = EVMASK[c["mask"]](D)
            for t in self.tiers:
                Kt = self.K[(kind, t)]
                ei = np.flatnonzero(m_ev); j = Kt["loc"][ev_pos[ei]]; ok = j >= 0
                ei, j = ei[ok], j[ok]; s = Kt["st"][j]; ok = Kt["elig"][s]; ei, j, s = ei[ok], j[ok], s[ok]
                ns = self.nscope(c, t)
                W0, W1 = self.weights(c, t, Kt["isg_real"], ns)
                for g in GROUPS:
                    if g == "both": sel_s = Kt["elig"]; okg = np.ones(len(ei), bool)
                    else:
                        rc = 0 if g == "N" else 1; sel_s = Kt["elig"] & (Kt["s_res"] == rc); okg = Kt["s_res"][s] == rc
                    eprot = np.unique(Kt["s_prot"][sel_s]); m = len(eprot)
                    pmap = np.full(nprot, -1, np.int64); pmap[eprot] = np.arange(m)
                    gi, gj, gs = ei[okg], j[okg], s[okg]; lp = pmap[Kt["s_prot"][gs]]; can = ev_can[gi]
                    w1e = W1[gs, 0]; inf = [len(np.unique(lp))]; exp = [float(w1e.sum())]
                    if ns > 1:
                        w1c = W1[gs, 1 + can]
                        for cc in range(ncan):
                            mk = can == cc; inf.append(len(np.unique(lp[mk]))); exp.append(float(w1c[mk].sum()))
                    tested = [(a >= S8.MIN_INF and b >= S8.MIN_EXP) for a, b in zip(inf, exp)]
                    self.S[(c["name"], t, g)] = dict(ei=gi, j=gj, s=gs, lp=lp, can=can, m=m, ns=ns, inf=inf, exp=exp, tested=tested)

    def nscope(self, c, t): return 1 + self.ncan if (c["cancer"] and t == 1 and self.ncan > 0) else 1

    def weights(self, c, t, L, ns):
        """(w0, w1), arrays strata x scopes: E0/E and E1/E with E1, E0 = exposure of glycosite and background positions (0 where E = 0)."""
        Kt = self.K[(c["kind"], t)]; nst = len(Kt["n1"])
        if c["mode"] == "uniform":
            w0 = Kt["n0"] / Kt["tot"]; w1 = Kt["n1"] / Kt["tot"]     # P2_08: counts are kept by every shuffle, so the weights are label-independent
            return np.repeat(w0[:, None], ns, axis=1), np.repeat(w1[:, None], ns, axis=1)
        st = Kt["st"]; Lf = L.astype(np.float64); wp = self.wpos[c["spec"]][Kt["pool"]]
        E1 = [np.bincount(st, weights=wp * Lf, minlength=nst)]; E0 = [np.bincount(st, weights=wp * (1.0 - Lf), minlength=nst)]
        if ns > 1:
            ep = Kt["ent_p"]; key = st[ep] * 192 + L[ep].astype(np.int64) * 96 + Kt["ent_ch"].astype(np.int64)
            A = np.bincount(key, minlength=nst * 192).astype(np.float64).reshape(nst, 2, 96)
            Ec = A @ self.R[c["spec"]][:, 1:]
            E0.append(Ec[:, 0, :]); E1.append(Ec[:, 1, :])
        E0 = np.column_stack(E0); E1 = np.column_stack(E1); E = E0 + E1
        w0 = np.divide(E0, E, out=np.zeros_like(E), where=E > 0); w1 = np.divide(E1, E, out=np.zeros_like(E), where=E > 0)
        return w0, w1

    def real_labels(self): return {k: v["isg_real"] for k, v in self.K.items()}

    def shuffle(self, r_main, r_strict):
        """C2 shuffle, same steps and random draws as SqnPAN_08 run_repeat."""
        lab = {}
        for kind, rng in (("main", r_main), ("strict", r_strict)):
            for t in self.tiers:
                if (kind, t) not in self.K: continue
                Kt = self.K[(kind, t)]; st = Kt["st"]
                u = rng.random(Kt["npool"]); order = np.lexsort((u, st)); sts = st[order]
                start = np.searchsorted(sts, sts, side="left"); rank = np.arange(len(sts)) - start
                L = np.empty(len(sts), bool); L[order] = rank < Kt["n1"][sts]; lab[(kind, t)] = L
        return lab

    def evaluate(self, lab, r_boot, upl, levels, extras=False):
        """The labelling block of SqnPAN_08 run_repeat / SqnPAN_09 evaluate, with scope-specific weights. Returns records in SqnPAN_08 order."""
        B = S8.B; out = []
        for c in self.cfgs:
            cfg = c["name"]; kind = c["kind"]; recs = []
            for t in self.tiers:
                Kt = self.K[(kind, t)]; L = lab[(kind, t)]
                W0, W1 = self.weights(c, t, L, self.nscope(c, t))
                for g in GROUPS:
                    Sx = self.S[(cfg, t, g)]; m = Sx["m"]; ns = Sx["ns"]
                    isg = L[Sx["j"]]; w0p = W0[Sx["s"], 0]; w1p = W1[Sx["s"], 0]
                    if ns > 1: w0c = W0[Sx["s"], 1 + Sx["can"]]; w1c = W1[Sx["s"], 1 + Sx["can"]]
                    u = upl[t - 1][Sx["ei"]] if upl is not None else None
                    key0 = Sx["lp"] * ns; key1 = key0 + 1 + Sx["can"] if ns > 1 else None
                    def agg(wp_, wc_):
                        a = np.bincount(key0, weights=wp_, minlength=m * ns)
                        if ns > 1: a = a + np.bincount(key1, weights=wc_, minlength=m * ns)
                        return a.reshape(m, ns)
                    cols = []
                    for _name, rr in levels:
                        if rr == 1.0: cnt = np.ones(len(Sx["ei"]))
                        elif rr < 1.0: cnt = (u < rr).astype(np.float64)
                        else: cnt = 1.0 + (u < (rr - 1.0))
                        cols.append(agg(np.where(isg, cnt * w0p, 0.0), np.where(isg, cnt * w0c, 0.0) if ns > 1 else None))
                    cols.append(agg(np.where(isg, 0.0, w1p), np.where(isg, 0.0, w1c) if ns > 1 else None))
                    M = np.concatenate(cols, axis=1); BS = S8.boot_product(r_boot, m, M); obs = M.sum(axis=0)
                    if extras:
                        sel = Kt["elig"] if g == "both" else (Kt["elig"] & (Kt["s_res"] == (0 if g == "N" else 1)))
                        pos_g = int(Kt["n1"][sel].sum()); pos_b = int(Kt["n0"][sel].sum())
                    for li, (name, _rr) in enumerate(levels):
                        for sc in range(ns):
                            cn = li * ns + sc; cd = len(levels) * ns + sc
                            st = S8.unit_stats(obs[cn], obs[cd], BS[:, cn], BS[:, cd])
                            fam = "primary" if (t == 1 and g == "both" and sc == 0) else ("pooled_other" if sc == 0 else f"cancer_T1_{g}")
                            rec = dict(level=name, cfg=cfg, tier=t, group=g, scope="pooled" if sc == 0 else self.projects[sc - 1], family=fam,
                                       tested=Sx["tested"][sc], n_inf_prot=Sx["inf"][sc], expected=Sx["exp"][sc], m=m, **st)
                            if extras:
                                mk = np.ones(len(isg), bool) if sc == 0 else (Sx["can"] == sc - 1)
                                rec.update(ev_g=int(isg[mk].sum()), ev_b=int((~isg[mk]).sum()), pos_g=pos_g, pos_b=pos_b)
                            recs.append(rec)
            for name, _rr in levels:
                byfam = collections.defaultdict(list)
                for r in recs:
                    if r["level"] != name: continue
                    r["fdr"] = None
                    if r["tested"] and r["num"] > 0 and r["den"] > 0 and not math.isnan(r["p"]): byfam[r["family"]].append(r)
                for fam, rl in byfam.items():
                    for r, q in zip(rl, S8.bh([r["p"] for r in rl])): r["fdr"] = float(q)
            for r in recs:
                if not r["tested"]: lab_ = NOT_TESTED
                elif r["num"] == 0 or r["den"] == 0: lab_ = "Inconclusive (zero count)"
                elif r["n_dropped"] > S8.NAN_MAX * B: lab_ = "Inconclusive (unstable bootstrap)"
                elif r["fdr"] is not None and r["fdr"] < S8.FDR_CUT and abs(r["est"]) >= S8.CHG: lab_ = "Enriched" if r["est"] > 0 else "Depleted"
                elif r["lo90"] >= -S8.EQ and r["hi90"] <= S8.EQ: lab_ = "Equivalent"
                else: lab_ = "Inconclusive"
                r["label"] = lab_; out.append(r)
        return out

# ---------------- lines, controls, capability, real runs ----------------
def p208_line(rep, r):
    return ("\t".join(fmt(x) for x in [rep, r["level"], r["cfg"], f"T{r['tier']}", r["group"], r["scope"], r["family"], bool(r["tested"]), r["n_inf_prot"]])
            + "\t" + "\t".join([fmt(r["expected"]), fmt(r["m"]), fmt(r["num"]), fmt(r["den"]), fmt(r["est"]), fmt(r["lo95"]), fmt(r["hi95"]), fmt(r["lo90"]), fmt(r["hi90"]),
                                 fmt(r["p"], True), fmt(r["fdr"], True), fmt(r["n_inf"]), fmt(r["n_dropped"]), r["label"]]))

def real_line(run, r):
    return "\t".join([fmt(run), r["cfg"], f"T{r['tier']}", r["group"], r["scope"], r["family"], fmt(bool(r["tested"])), fmt(r["n_inf_prot"]), fmt(r["expected"]), fmt(r["m"]),
                      fmt(r["num"]), fmt(r["den"]), fmt(r["est"]), fmt(r["lo95"]), fmt(r["hi95"]), fmt(r["lo90"]), fmt(r["hi90"]), fmt(r["p"], True), fmt(r["fdr"], True),
                      fmt(r["n_inf"]), fmt(r["n_dropped"]), r["label"], fmt(r["ev_g"]), fmt(r["ev_b"]), fmt(r["pos_g"]), fmt(r["pos_b"])])

def control_repeat(eng, base, rep):
    """One C2/C3 repeat with the random streams of SqnPAN_08 (main shuffle, strict shuffle, planting, bootstrap)."""
    r_main, r_strict, r_plant, r_boot = [np.random.default_rng(x) for x in np.random.SeedSequence(base + rep).spawn(4)]
    lab = eng.shuffle(r_main, r_strict); upl = r_plant.random((len(S8.TIERS), eng.nev))
    return [p208_line(rep, r) for r in eng.evaluate(lab, r_boot, upl, S8.LEVELS)]

_W = {}
def _winit(npz_path, configs, restricted): _W["eng"] = Engine(load_npz(npz_path), configs, restricted)
def _wrep(job):
    run, base, rep, parts = job
    out = os.path.join(parts, f"run{run}_rep{rep:02d}.tsv")
    if os.path.exists(out): return rep, 0.0, True
    t0 = time.time(); lines = control_repeat(_W["eng"], base, rep); tmp = out + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as o: o.write("\n".join(lines) + "\n")
    os.replace(tmp, out); return rep, time.time() - t0, False

def check_seed_blocks(log, block, subs, n_rep=20):
    """Sub-bases inside one reserved block of 100,000 (Amendment 4 rule 5); repeat seeds base + 1 to base + n_rep never reach the next sub-base."""
    s = sorted(set(subs))
    for b in s:
        if not (block <= b and b + n_rep < block + 100000): raise RuntimeError(f"seed sub-base {b} leaves the reserved block {block}")
    for a, b in zip(s, s[1:]):
        if b - a <= n_rep: raise RuntimeError(f"seed sub-bases {a} and {b} overlap")
    log(f"seed sub-bases {s} inside the reserved block {block} to {block + 99999}: no overlap of repeat seeds (base + 1 to base + {n_rep})")

def run_controls(log, npz_path, configs, restricted, bases, parts, table_fn, workers):
    os.makedirs(parts, exist_ok=True); n_rep = S8.N_REP; cols = S8.COLS
    pool = mp.get_context("spawn").Pool(workers, initializer=_winit, initargs=(npz_path, configs, restricted)) if workers > 1 else None
    if pool is None: _winit(npz_path, configs, restricted)
    try:
        for run in (1, 2, 3):
            base = bases[run]; jobs = [(run, base, r, parts) for r in range(1, n_rep + 1)]; t0 = time.time(); done = 0
            log(f"--- run{run}: seed base {base}, {n_rep} repeats ---")
            it = pool.imap_unordered(_wrep, jobs) if pool else map(_wrep, jobs)
            for rep, sec, skipped in it:
                done += 1
                log(f"  run{run} repeat {rep}: {'already saved, skipped' if skipped else f'{sec:.1f} s'} | {done}/{n_rep} | elapsed {time.time() - t0:.0f} s")
            body = []
            for r in range(1, n_rep + 1):
                with open(os.path.join(parts, f"run{run}_rep{r:02d}.tsv"), encoding="utf-8") as f: body.append(f.read())
            with open(table_fn(run), "w", newline="", encoding="utf-8") as o: o.write("\t".join(cols) + "\n" + "".join(body))
            log(f"run{run} table {table_fn(run)} | md5 {H.md5_file(table_fn(run))} | seconds {time.time() - t0:.0f}")
            rows = list(csv.DictReader(open(table_fn(run), newline="", encoding="utf-8"), delimiter="\t"))
            for cfg in [c["name"] for c in configs]:
                sh = []
                for rep in range(1, n_rep + 1):
                    u = [r for r in rows if r["config"] == cfg and r["level"] == "1.0" and r["repeat"] == str(rep) and r["tested"] == "1"]
                    sh.append(f"{sum(1 for r in u if r['label'] in ('Enriched', 'Depleted'))}/{len(u)}")
                log(f"  run{run} C2 {cfg}: Changed / tested units per repeat: {' '.join(sh)}")
    finally:
        if pool: pool.close(); pool.join()

def capability(log, table_fn, cfg_names, cap_path):
    """SqnPAN_08 section 3 rules (C2 per set-up and run; capability and Equivalent-allowed per unit from runs 1 and 3). Returns (cap, c2ok, run1 equals run2)."""
    N_REP = S8.N_REP
    m1, m2, m3 = (H.md5_file(table_fn(k)) for k in (1, 2, 3)); same12 = m1 == m2
    log(f"run1 md5 {m1} | run2 md5 {m2} | identical: {'yes' if same12 else 'NO'} | run3 md5 {m3}")
    T = {}
    for k in (1, 3):
        T[k] = collections.defaultdict(list)
        with open(table_fn(k), newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"): T[k][(r["config"], r["tier"], r["group"], r["scope"])].append(r)
    c2ok = {}
    for k in (1, 3):
        for cfg in cfg_names:
            worst = 0.0; nfail = 0; fam = collections.Counter(); ntest = None
            for rep in range(1, N_REP + 1):
                u = [r for key, rl in T[k].items() if key[0] == cfg for r in rl if r["level"] == "1.0" and r["repeat"] == str(rep) and r["tested"] == "1"]
                ch = [r for r in u if r["label"] in ("Enriched", "Depleted")]
                for r in ch: fam[r["family"]] += 1
                share = len(ch) / len(u) if u else 0.0; worst = max(worst, share); ntest = len(u)
                if share > S8.C2_MAX: nfail += 1
            c2ok[(k, cfg)] = nfail == 0
            log(f"C2 run{k} {cfg}: tested units {ntest} | largest Changed share {worst:.3f} (limit {S8.C2_MAX}) | repeats over the limit {nfail} | Changed by family {dict(fam)} | {'PASSED' if nfail == 0 else 'FAILED'}")
    out = []; cap = {}; lab_count = collections.Counter()
    for key in sorted(T[1]):
        cfg = key[0]; rl1 = T[1][key]; rl3 = T[3][key]; tested = rl1[0]["tested"] == "1"
        if (rl3[0]["tested"] == "1") != tested: raise RuntimeError(f"tested status differs between run1 and run3 for {key}")
        rowv = list(key) + [rl1[0]["family"], int(tested), rl1[0]["n_informative"], rl1[0]["expected_glyco"]]; pw = []; eqa = []
        for rl in (rl1, rl3):
            def n_lab(level, labs): return sum(1 for r in rl if r["level"] == level and r["label"] in labs)
            right = {lv: n_lab(lv, (d,)) for lv, d in list(S8.STRONG.items()) + list(S8.WEAK.items())}
            eq1 = n_lab("1.0", ("Equivalent",)); eqs = [n_lab(lv, ("Equivalent",)) for lv in S8.STRONG]
            null = [float(r["log2rr"]) for r in rl if r["level"] == "1.0" and r["log2rr"] not in ("nan", "inf", "-inf")]
            cov90 = sum(1 for r in rl if r["level"] == "1.0" and r["lo90"] != "nan" and float(r["lo90"]) <= 0 <= float(r["hi90"]))
            cov95 = sum(1 for r in rl if r["level"] == "1.0" and r["lo95"] != "nan" and float(r["lo95"]) <= 0 <= float(r["hi95"]))
            med = float(np.median(null)) if null else float("nan")
            pw.append(all(right[lv] >= S8.PASS_N for lv in S8.STRONG)); eqa.append(eq1 >= S8.PASS_N and all(x <= S8.EQ_STRONG_MAX for x in eqs))
            rowv += [right["1.5"], right["0.667"], right["1.25"], right["0.8"], eq1, eqs[0], eqs[1], fmt(med), cov90, cov95]
        if not tested: label = NOT_TESTED; eq_ok = False
        elif not (c2ok[(1, cfg)] and c2ok[(3, cfg)]): label = "Inconclusive, test not sensitive (C2 failed)"; eq_ok = False
        elif all(pw): label = CAN; eq_ok = all(eqa)
        else: label = NOTSENS; eq_ok = False
        rowv += [label, "yes" if eq_ok else "no"]; out.append(rowv); cap[key] = (label, eq_ok)
        lab_count[(cfg, "pooled" if key[3] == "pooled" else "cancer type", label, eq_ok)] += 1
    hdr = ["config", "tier", "group", "scope", "family", "tested", "n_informative", "expected_glyco"]
    for k in (1, 3): hdr += [f"run{k}_{x}" for x in ("changed_right_1.5", "changed_right_0.667", "changed_right_1.25", "changed_right_0.8", "equivalent_1.0",
                                                    "equivalent_1.5", "equivalent_0.667", "median_null_log2rr", "null_cover90", "null_cover95")]
    hdr += ["capability", "equivalent_allowed"]
    with open(cap_path, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(hdr); w.writerows(out)
    log(f"capability table {cap_path}: {len(out)} units | md5 {H.md5_file(cap_path)}")
    for k in sorted(lab_count): log(f"  {k}: {lab_count[k]}")
    log("--- pooled units (C3 counts of 20: right-direction Changed at 1.5, 0.667, 1.25, 0.8; Equivalent at 1.0, 1.5, 0.667; run1 then run3) ---")
    for rowv in out:
        if rowv[3] == "pooled": log("  " + " | ".join(str(x) for x in rowv[:4]) + " | run1 " + " ".join(str(x) for x in rowv[8:15]) + " | run3 " + " ".join(str(x) for x in rowv[18:25]) + f" | {rowv[-2]} | Eq {rowv[-1]}")
    return cap, c2ok, same12

def read_cap(path):
    with open(path, newline="", encoding="utf-8") as f:
        return {(r["config"], r["tier"], r["group"], r["scope"]): (r["capability"], r["equivalent_allowed"] == "yes") for r in csv.DictReader(f, delimiter="\t")}

def final_label(raw, cap, eq_ok):
    """SqnPAN_09 decision 1A. Returns (final label, flag)."""
    if raw == NOT_TESTED: return NOT_TESTED, ""
    sens = cap == CAN
    if raw in ("Enriched", "Depleted"): return raw, ("" if sens else FLAG_LOW)
    if raw == "Equivalent": return ("Equivalent", "") if eq_ok else ("Inconclusive, Equivalent not allowed", "")
    return (raw, "") if sens else (NOTSENS, "")

def real_runs(log, eng, bases, path_fn):
    lab = eng.real_labels()
    for run in (1, 2, 3):
        base = bases[run]; t0 = time.time()
        recs = eng.evaluate(lab, np.random.default_rng(np.random.SeedSequence(base)), None, LEVEL1, extras=True)
        with open(path_fn(run), "w", newline="", encoding="utf-8") as o: o.write("\t".join(COLS_REAL) + "\n" + "".join(real_line(run, r) + "\n" for r in recs))
        log(f"real run{run} (seed base {base}): {len(recs)} units | {time.time() - t0:.0f} s | md5 {H.md5_file(path_fn(run))} | content md5 without the run column {body_md5(path_fn(run))}")

def body_md5(path):
    h = hashlib.md5()
    with open(path, encoding="utf-8") as f:
        for line in f: h.update(line.split("\t", 1)[1].encode("utf-8"))
    return h.hexdigest()

def read_run(path):
    with open(path, newline="", encoding="utf-8") as f: return {(r["config"], r["tier"], r["group"], r["scope"]): r for r in csv.DictReader(f, delimiter="\t")}

def ff(x):
    try: v = float(x)
    except Exception: return str(x)
    return "nan" if math.isnan(v) else ("inf" if v == float("inf") else ("-inf" if v == float("-inf") else f"{v:.3f}"))

def assemble(log, path_fn, cap, cfg_order, out_path, suffix="", extra_hdr=(), extra_fn=None):
    """Final labels (SqnPAN_09 rule 1A) from run 1, run 3 beside it; writes the results table. Returns (runs, fin)."""
    runs = {k: read_run(path_fn(k)) for k in (1, 2, 3)}; m = {k: body_md5(path_fn(k)) for k in (1, 2, 3)}
    log(f"real-run content md5 without the run column: run1 {m[1]} | run2 {m[2]} | identical: {'yes' if m[1] == m[2] else 'NO'} | run3 {m[3]}")
    if set(cap) != set(runs[1]) or set(cap) != set(runs[3]): raise RuntimeError("unit lists of the real runs and of the capability table differ")
    def ukey(k): return (cfg_order.index(k[0]), k[1], GORDER[k[2]], 0 if k[3] == "pooled" else 1, k[3])
    keys = sorted(cap, key=ukey); fin = {1: {}, 3: {}}
    for key in keys:
        for k in (1, 3):
            lab, flag = final_label(runs[k][key]["raw_label"], cap[key][0], cap[key][1])
            if suffix and lab != NOT_TESTED: lab = lab + suffix
            fin[k][key] = (lab, flag)
    hdr = ["config", "tier", "group", "scope", "family", "tested", "n_informative", "expected_glyco", "n_events_glyco", "n_events_bg", "n_positions_glyco", "n_positions_bg",
           "log2rr", "rate_ratio", "lo95", "hi95", "lo90", "hi90", "p", "fdr", "raw_label", "capability", "equivalent_allowed", "final_label", "flag",
           "run3_log2rr", "run3_lo95", "run3_hi95", "run3_fdr", "run3_raw_label", "run3_final_label"] + list(extra_hdr)
    out = []
    for key in keys:
        r = runs[1][key]; r3 = runs[3][key]
        rr = fmt(2.0 ** float(r["log2rr"])) if r["log2rr"] not in ("nan",) else "nan"
        out.append([key[0], key[1], key[2], key[3], r["family"], r["tested"], r["n_informative"], r["expected_glyco"], r["n_events_glyco"], r["n_events_bg"], r["n_positions_glyco"],
                    r["n_positions_bg"], r["log2rr"], rr, r["lo95"], r["hi95"], r["lo90"], r["hi90"], r["p"], r["fdr"], r["raw_label"], cap[key][0],
                    "yes" if cap[key][1] else "no", fin[1][key][0], fin[1][key][1], r3["log2rr"], r3["lo95"], r3["hi95"], r3["fdr"], r3["raw_label"], fin[3][key][0]]
                   + (list(extra_fn(key)) if extra_fn else []))
    with open(out_path, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(hdr); w.writerows(out)
    log(f"results table {out_path} | {len(out)} units | md5 {H.md5_file(out_path)}")
    nd = sum(1 for key in keys if runs[1][key]["raw_label"] != runs[3][key]["raw_label"]); nf = sum(1 for key in keys if fin[1][key][0] != fin[3][key][0])
    log(f"seed comparison (run 3 against run 1): raw label differs in {nd} units, final label differs in {nf} units")
    return runs, fin

def unit_line(runs, fin, cap, key):
    r = runs[1][key]; fl, flag = fin[1][key]
    return (f"log2RR {ff(r['log2rr'])} (95% CI {ff(r['lo95'])} to {ff(r['hi95'])}; 90% CI {ff(r['lo90'])} to {ff(r['hi90'])}) FDR {ff(r['fdr']) if r['fdr'] else '-'}"
            f" | events glyco {r['n_events_glyco']} bg {r['n_events_bg']} | raw {r['raw_label']} | capability {cap[key][0]} | FINAL {fl}" + (f" [{flag}]" if flag else ""))

def run_chain(prefix, sections, start):
    keys = list(sections)
    if start not in keys: raise SystemExit(f"unknown section {start}; use all or one of {keys}")
    for k in keys[keys.index(start):]:
        name, fn = sections[k]; ok = {"v": False}
        def wrapped(log, fn=fn): fn(log); ok["v"] = True
        H.run_section(prefix, k, name, wrapped)
        if not ok["v"]: print(f"SECTION {k} DID NOT FINISH: the chain stops here (see the log)", flush=True); return False
    return True
# ---- END OF HELPERS ----
