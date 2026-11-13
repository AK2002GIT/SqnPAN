# SqnPAN_20_hotspots.py - B4 single-site hotspots (Amendment 4; added after the primary analysis)
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# PRE-SPECIFICATION (dated 2026-11-07). Post hoc: no label of SqnPAN_09, SqnPAN_14 or SqnPAN_19 changes.
# Question: are individual T1 glycosite residues (the N, and the S/T at +2 of a T1 sequon) hit by somatic substitutions in more cases than their trinucleotide
#   context predicts, given the mutation burden of their own protein in each cancer type (single-site hotspots)?
# Unit: one T1 glycosite position of the SqnPAN_14 position list (role 1, with codon opportunities). Recurrence k = number of cases with a counted substitution there
#   (the SqnPAN_08 / Engine events: one per case and position; at S/T only substitutions to a non-S/T residue count).
# Expected: pool = positions of the Engine T1 pool (role 0 or 1, codon opportunities). Stratum = protein x residue type (N, S/T). For cancer type c the n[s,c] events of
#   stratum s are shared over its positions by their context rates w[j,c] (sum of the A1 rates R_main (main) or R_nohyper (nohyper) over the position's counted
#   single-nucleotide changes): pi[j,c] = w[j,c] / W[s,c]. Under the null X_j = sum over c of Binomial(n[s,c], pi[j,c]); p = exact P(X_j >= k). The protein's own burden
#   per cancer type is kept (conditioning on n[s,c]).
# Families and labels: T1 family = T1 positions in strata with at least one event of the set-up; BH over the family. Hotspot = BH < 0.05 AND k >= 3; otherwise Not hotspot.
#   Set-ups: main (primary; all cases) and nohyper (cases with more than 1,000 missense rows left out; own family). A main Hotspot that is not a nohyper Hotspot is
#   flagged "hypermutator-driven". Background comparison (descriptive): the same test at background positions (role 0) of the same strata, own BH family; shares of
#   Hotspot and of (p < 0.001 and k >= 3) at T1 against background, odds ratio with a Woolf 95% CI (0.5 added to each cell).
# Controls (section 2, judged before the real result is computed in section 3):
#   Calibration: background positions of the main family outside the 299 Bailey driver genes; real count with p < 0.001 against the mean count in the first 50 C2 main
#     repeats. PASSED if the ratio is <= 2.0 (p < 0.01 reported). If FAILED, every T1 label gets the suffix _calibration_FAILED.
#   C2 null: 200 repeats; in each, every event of cancer c in stratum s is moved to a position of s drawn with probabilities pi[.,c]; then the full T1 test. PASSED for a
#     set-up if the share of repeats with at least one T1 Hotspot is <= 0.07 (main and nohyper each). If FAILED, the set-up's labels read "Inconclusive (null control failed)".
#   C3 planted (main): 20 repeats x 100 random positions of the T1 family; each gets d extra cases (d = 3, 5, 8; the cancer type of each extra case drawn from the real events
#     of its stratum), p recomputed in the affected strata, BH over the family. PASSED if recovery (planted position labelled Hotspot) is >= 0.80 at d = 8 (3 and 5
#     reported). If FAILED, "Not hotspot" reads "Not hotspot, test not sensitive".
#   Positive control (main, real data): 8 known driver hotspot residues, all background positions: CTNNB1 S33, S37, T41, S45; U2AF1 S34; FGFR3 S249; ERBB2 S310;
#     PIK3CA N345 (residues checked against UniProt 2026_03; TP53 S241 left out because it lies in a sequon). PASSED if at least 6 of 8 have k >= 3 and
#     p < 0.05 / 8. If FAILED, every T1 label gets the suffix _positive_control_FAILED. Not blind: these hotspots are well known.
# Bailey driver overlap of the hotspot genes: descriptive (alias step of SqnPAN_12).
# Section 4 writes the inputs of the independent ActiveDriver test (SqnPAN_21_activedriver.R, pre-specification in its header).
# Not blind: the pooled Q2 and A1 results (Equivalent) were read before; no single-site result was read before.
# Seeds (block 33300000 to 33399999 reserved for B4): C2 33300000 (repeat r: base + r), C3 33320000, ActiveDriver null proteins and sites 33340000 (+ 1 to 3),
#   ActiveDriver runs 33360000 and 33380000 (used in SqnPAN_21).
# Usage (from the repository folder):  python scripts/SqnPAN_20_hotspots.py all
#   "all" runs sections 1 to 4 and stops at the first section that does not finish; "2", "3" or "4" restart from that section.
# Logs: results\logs\P2_20_section<k>_<name>.txt. Outputs listed in each section.
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import sys, re, csv, gzip, math, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import SqnPAN_05_helpers as H
import SqnPAN_08_selection_controls as S8
import SqnPAN_14_context_helpers as C

BLOCK = 33300000
SEEDS = dict(c2=33300000, c3=33320000, adnull=33340000, adrun=33360000, adpos=33380000)
ENG14 = os.path.join(H.PROC_OUT, "P2_14_engine.npz"); ENG14_MD5 = "11da2bc41ad5271d27d099b1dace5f3e"
POS14 = os.path.join(H.PROC_OUT, "P2_14_positions.npz"); POS14_MD5 = "85e2842262142b016493b7a3aa09e04a"
OK14_2 = os.path.join(H.PROC_OUT, "P2_14_section2_ok.txt")
MD5_IN = [(ENG14, ENG14_MD5), (POS14, POS14_MD5), (H.UNIPROT_MAIN, S8.MD5[H.UNIPROT_MAIN]), (H.BAILEY, "144a9ccbbf038236decc1226df24620b")]
MD5_AD = [(S8.SITES, S8.MD5[S8.SITES]), (S8.PROTS, S8.MD5[S8.PROTS]), (S8.MUTS, S8.MD5[S8.MUTS]), (S8.RESID, S8.MD5[S8.RESID])]
OK1 = os.path.join(H.PROC_OUT, "P2_20_section1_ok.txt"); OK2 = os.path.join(H.PROC_OUT, "P2_20_section2_ok.txt"); OK3 = os.path.join(H.PROC_OUT, "P2_20_section3_ok.txt")
VERDICT = os.path.join(H.PROC_OUT, "P2_20_controls_verdict.txt")
T_C2 = os.path.join(H.TABLES, "P2_20_C2_null.tsv"); T_C3 = os.path.join(H.TABLES, "P2_20_C3_planted.tsv"); T_POS = os.path.join(H.TABLES, "P2_20_positive_control.tsv")
T_HOT = os.path.join(H.TABLES, "P2_20_T1_hotspots.tsv"); T_SUM = os.path.join(H.TABLES, "P2_20_summary.tsv"); T_BAIL = os.path.join(H.TABLES, "P2_20_bailey_overlap.tsv")
T_BGH = os.path.join(H.TABLES, "P2_20_background_hotspots.tsv")
AD = {k: os.path.join(H.PROC_OUT, f"P2_20_ad_{k}.tsv") for k in ("sequences", "disorder", "sites_real", "sites_null", "sites_pos", "proteins")}
AD["mutations"] = os.path.join(H.PROC_OUT, "P2_20_ad_mutations.tsv.gz")
N_C2, N_CAL, N_C3, N_PLANT, DS = 200, 50, 20, 100, (3, 5, 8)
FDR_CUT, KMIN, C2_MAX, CAL_MAX, C3_MIN = 0.05, 3, 0.07, 2.0, 0.80
POSCTRL = [("CTNNB1", "P35222", 33, "S"), ("CTNNB1", "P35222", 37, "S"), ("CTNNB1", "P35222", 41, "T"), ("CTNNB1", "P35222", 45, "S"), ("U2AF1", "Q01081", 34, "S"),
           ("FGFR3", "P22607", 249, "S"), ("ERBB2", "P04626", 310, "S"), ("PIK3CA", "P42336", 345, "N")]
ALIAS = [("FAM46D", "TENT5D", "Q8NEK8"), ("H3F3A", "H3-3A", "P84243"), ("H3F3C", "H3-5", "Q6NXT2"), ("HIST1H1C", "H1-2", "P16403"), ("HIST1H1E", "H1-4", "P10412"),
         ("RQCD1", "CNOT9", "Q92600"), ("TCEB1", "ELOC", "Q15369"), ("WHSC1", "NSD2", "O96028")]   # the P2_12 alias step

def need_ok(path, what):
    if not os.path.exists(path) or not open(path).read().startswith("OK"): raise RuntimeError(f"{what} has not finished ({os.path.basename(path)} missing)")

def check_md5(log, pairs):
    for path, rec in pairs:
        got = H.md5_file(path); log(f"md5 {os.path.basename(path)}: {got} ({'same as recorded' if got == rec else 'DIFFERENT, recorded ' + rec})")
        if got != rec: raise RuntimeError("an input md5 differs from the recorded value: stop and decide")

def bailey_genes(log, U):
    rows = H.read_xlsx(H.BAILEY, "Table S1"); hdr = rows[3]
    gi = hdr.index("Gene"); genes = [r[gi] for r in rows[4:] if len(r) > gi and r[gi]]
    log(f"Bailey Table S1 rows with a gene: {len(genes)} (expected 739); unique genes {len(set(genes))} (expected 299)")
    if len(genes) != 739 or len(set(genes)) != 299: raise RuntimeError("Bailey Table S1 counts differ from the recorded counts")
    drv = set(genes)
    for old, new, acc in ALIAS:
        ok = acc in U and U[acc]["gene"] == new
        if ok and old in drv: drv.discard(old); drv.add(new)
        log(f"  alias {old} -> {new} ({acc}): guard {'ok' if ok else 'FAILED (kept as is)'}")
    return drv

# ---------------- the model ----------------
class Model:
    """Pool, strata, per-cancer context weights and event counts for one set-up (main or nohyper)."""
    def __init__(self, D, Pz, setup):
        role = D["P_role"].astype(np.int8); ctx = Pz["P_ctx"]; prot = D["P_prot"].astype(np.int64); res = D["P_res"].astype(np.int64)
        inpool = ctx & ((role == 0) | (role == 1)); strat_all = prot * 2 + res
        idx = np.flatnonzero(inpool); o = np.argsort(strat_all[idx], kind="stable"); self.pool = idx[o]      # pool sorted by stratum
        st_raw = strat_all[self.pool]; self.ukeys, self.st = np.unique(st_raw, return_inverse=True); self.st = self.st.ravel()
        self.nst = len(self.ukeys); self.npool = len(self.pool)
        self.first = np.searchsorted(self.st, np.arange(self.nst), "left"); self.last = np.searchsorted(self.st, np.arange(self.nst), "right")
        self.role = role[self.pool]; self.prot = prot[self.pool]; self.res = res[self.pool]; self.posn = Pz["P_posn"][self.pool].astype(np.int64)
        self.loc = np.full(len(role), -1, np.int64); self.loc[self.pool] = np.arange(self.npool)
        self.projects = [str(x) for x in D["projects"]]; nc = len(self.projects); self.nc = nc
        R = D["R_main" if setup == "main" else "R_nohyper"].astype(np.float64); O = D["OPP"][self.pool].astype(np.int64); ok = O >= 0; Oc = np.maximum(O, 0)
        self.w = np.zeros((self.npool, nc), np.float32)
        for c in range(nc): self.w[:, c] = np.where(ok, R[Oc, 1 + c], 0.0).sum(axis=1)
        self.W = np.zeros((self.nst, nc))
        for c in range(nc): self.W[:, c] = np.bincount(self.st, weights=self.w[:, c].astype(np.float64), minlength=self.nst)
        m = D["ev_main"] & (~D["ev_hyper"] if setup == "nohyper" else True)
        ei = np.flatnonzero(m); j = self.loc[D["ev_pos"][ei].astype(np.int64)]; keep = j >= 0
        self.ev_j = j[keep]; self.ev_c = D["ev_cancer"][ei[keep]].astype(np.int64); self.ev_s = self.st[self.ev_j]
        self.n_out = int((~keep).sum())
        self.n = np.bincount(self.ev_s * nc + self.ev_c, minlength=self.nst * nc).reshape(self.nst, nc).astype(np.int64)
        self.ns_ev = self.n.sum(axis=1)
        self.k_real = np.bincount(self.ev_j, minlength=self.npool).astype(np.int64)
        self.cases_by_c = None
        self.prep_tail(self.n)

    def prep_tail(self, n):
        """L0[j] = log P(X_j = 0) and lam[j] = E[X_j] for the counts n."""
        self.L0 = np.zeros(self.npool); self.lam = np.zeros(self.npool)
        for c in range(self.nc):
            nj = n[self.st, c].astype(np.float64)
            if not nj.any(): continue
            pi = np.divide(self.w[:, c].astype(np.float64), self.W[self.st, c], out=np.zeros(self.npool), where=self.W[self.st, c] > 0)
            pi = np.minimum(pi, 1.0)
            with np.errstate(divide="ignore"): self.L0 += np.where(nj > 0, nj * np.log1p(-pi), 0.0)
            self.lam += nj * pi

    def tail(self, j, k, n):
        """Exact P(X_j >= k), X_j = sum over cancer types of Binomial(n[s,c], pi[j,c])."""
        if k <= 0: return 1.0
        if k == 1: return float(-np.expm1(self.L0[j])) if np.isfinite(self.L0[j]) else 1.0
        s = self.st[j]; cs = np.flatnonzero(n[s] > 0); pmf = np.zeros(k); pmf[0] = 1.0
        for c in cs:
            nn = int(n[s, c]); W = self.W[s, c]; p = float(self.w[j, c]) / W if W > 0 else 0.0
            if p <= 0: continue
            if p >= 1.0 - 1e-12:
                new = np.zeros(k)
                if nn < k: new[nn:] = pmf[:k - nn]
                pmf = new; continue
            m = min(k - 1, nn); b = np.empty(m + 1); b[0] = math.exp(nn * math.log1p(-p)); r = p / (1.0 - p)
            for x in range(m): b[x + 1] = b[x] * (nn - x) / (x + 1) * r
            pmf = np.convolve(pmf, b)[:k]
        return max(0.0, 1.0 - float(pmf.sum()))

    def pvals(self, k, sel, n=None):
        """p for the positions sel (indices into the pool) with counts k[sel]; p = 1 where k = 0."""
        n = self.n if n is None else n; out = np.ones(len(sel)); kk = k[sel]
        one = kk == 1; out[one] = np.where(np.isfinite(self.L0[sel[one]]), -np.expm1(self.L0[sel[one]]), 1.0)
        for i in np.flatnonzero(kk >= 2): out[i] = self.tail(int(sel[i]), int(kk[i]), n)
        return out

    def family(self, roleval):
        return np.flatnonzero((self.role == roleval) & (self.ns_ev[self.st] > 0))

    def simulate(self, rng):
        """C2: every event moved to a position of its stratum drawn with probabilities pi[., c]."""
        jj = np.empty(len(self.ev_j), np.int64)
        for c in range(self.nc):
            mk = np.flatnonzero(self.ev_c == c)
            if not len(mk): continue
            cum = np.cumsum(self.w[:, c].astype(np.float64)); start = np.concatenate([[0.0], cum])[self.first]
            s = self.ev_s[mk]; u = start[s] + rng.random(len(mk)) * self.W[s, c]
            j = np.searchsorted(cum, u, side="right"); j = np.minimum(np.maximum(j, self.first[s]), self.last[s] - 1); jj[mk] = j
        return np.bincount(jj, minlength=self.npool).astype(np.int64)

def call(p, k, fam):
    q = S8.bh(p) if len(p) else p
    return q, (q < FDR_CUT) & (k[fam] >= KMIN)

def load(log):
    D = C.load_npz(ENG14); Pz = C.load_npz(POS14)
    if [str(x) for x in D["accs"]] != [str(x) for x in Pz["accs"]]: raise RuntimeError("accession order differs between SqnPAN_14 engine and positions")
    return D, Pz

# ======================= Section 1: structure check =======================
def section1(log):
    need_ok(OK14_2, "P2_14 section 2")
    log("Pre-specification: header of this script (fixed 2026-11-07). Post hoc, added after the primary analysis.")
    log(f"Python {sys.version.split()[0]} | NumPy {np.__version__}")
    check_md5(log, MD5_IN); C.check_seed_blocks(log, BLOCK, SEEDS.values(), n_rep=N_C2)
    D, Pz = load(log); U = H.load_uniprot_main(); accs = [str(x) for x in D["accs"]]
    for setup in ("main", "nohyper"):
        t0 = time.time(); M = Model(D, Pz, setup)
        fam1 = M.family(1); fam0 = M.family(0)
        log(f"[{setup}] pool positions {M.npool:,} in {M.nst:,} strata | events in the pool {len(M.ev_j):,} (outside the pool {M.n_out:,}) | "
            f"strata with events {(M.ns_ev > 0).sum():,} | T1 family positions {len(fam1):,} | background family positions {len(fam0):,} | {time.time() - t0:.0f} s")
        if setup == "main":
            bad = 0
            for j in range(0, M.npool, 997):
                a = accs[M.prot[j]]; aa = U[a]["seq"][M.posn[j] - 1]
                if aa not in ("N" if M.res[j] == 0 else "ST"): bad += 1
            log(f"[{setup}] residue check on every 997th pool position: {bad} wrong residues (must be 0)")
            if bad: raise RuntimeError("position numbers do not match UniProt")
            if abs(float(M.lam.sum()) - len(M.ev_j)) > 1e-6 * len(M.ev_j) + 1: raise RuntimeError("expected counts do not sum to the event count")
            log(f"[{setup}] sum of expected counts {M.lam.sum():.1f} equals the event count (check of the weights)")
            a2i = {a: i for i, a in enumerate(accs)}
            for g, acc, pos, aa in POSCTRL:
                i = a2i[acc]; s = U[acc]["seq"]; hit = np.flatnonzero((M.prot == i) & (M.posn == pos))
                log(f"  positive control {g} {aa}{pos} ({acc}): UniProt residue {s[pos - 1]}, in pool {len(hit) == 1}, role {int(M.role[hit[0]]) if len(hit) else '-'}")
                if s[pos - 1] != aa or len(hit) != 1 or M.role[hit[0]] != 0: raise RuntimeError("a positive-control position is not a background pool position")
    with open(OK1, "w") as o: o.write("OK\n")

# ======================= Section 2: controls =======================
def section2(log):
    need_ok(OK1, "section 1"); check_md5(log, MD5_IN)
    D, Pz = load(log); U = H.load_uniprot_main(); accs = [str(x) for x in D["accs"]]; drv = bailey_genes(log, U)
    verdict = {}; c2rows = []
    for setup in ("main", "nohyper"):
        M = Model(D, Pz, setup); fam1 = M.family(1); fam0 = M.family(0)
        nondrv = np.array([U[accs[M.prot[j]]]["gene"] not in drv for j in fam0], bool); fam0n = fam0[nondrv]
        t0 = time.time(); anyhot = 0; cal_null = []
        for r in range(1, N_C2 + 1):
            rng = np.random.default_rng(SEEDS["c2"] + r + (0 if setup == "main" else 1000))
            if setup == "nohyper" and r == 1: C.check_seed_blocks(log, BLOCK, [SEEDS["c2"] + 1000], n_rep=N_C2)
            k = M.simulate(rng); p = M.pvals(k, fam1); q, hot = call(p, k, fam1); nh = int(hot.sum()); anyhot += nh > 0
            row = dict(setup=setup, repeat=r, T1_hotspots=nh, T1_min_fdr=float(q.min()) if len(q) else 1.0)
            if setup == "main" and r <= N_CAL:
                pb = M.pvals(k, fam0n); row["bg_nondriver_p001"] = int((pb < 0.001).sum()); row["bg_nondriver_p01"] = int((pb < 0.01).sum())
                cal_null.append((row["bg_nondriver_p001"], row["bg_nondriver_p01"]))
            c2rows.append(row)
            if r in (1, 10, 50, N_C2): log(f"[{setup}] C2 repeat {r}: T1 hotspots {nh} | elapsed {time.time() - t0:.0f} s")
        share = anyhot / N_C2; ok = share <= C2_MAX; verdict[f"C2_{setup}"] = ok
        log(f"[{setup}] C2 null: share of repeats with >= 1 T1 Hotspot {share:.3f} (rule <= {C2_MAX}) -> {'PASSED' if ok else 'FAILED'}")
        if setup != "main": continue
        # calibration (background outside Bailey drivers, real against the null mean)
        k = M.k_real; pb = M.pvals(k, fam0n); obs1 = int((pb < 0.001).sum()); obs2 = int((pb < 0.01).sum())
        m1 = float(np.mean([a for a, _ in cal_null])); m2 = float(np.mean([b for _, b in cal_null]))
        ratio = obs1 / m1 if m1 > 0 else float("inf"); ok = ratio <= CAL_MAX; verdict["calibration"] = ok
        log(f"calibration: background non-driver positions {len(fam0n):,}; p < 0.001 real {obs1} against null mean {m1:.1f} (ratio {ratio:.2f}, rule <= {CAL_MAX}) -> "
            f"{'PASSED' if ok else 'FAILED'}; p < 0.01 real {obs2} against {m2:.1f} (ratio {obs2 / m2 if m2 > 0 else float('nan'):.2f}, reported)")
        # positive control
        a2i = {a: i for i, a in enumerate(accs)}; prow = []; npass = 0; sel = []
        for g, acc, pos, aa in POSCTRL: sel.append(int(np.flatnonzero((M.prot == a2i[acc]) & (M.posn == pos))[0]))
        sel = np.array(sel); pp = M.pvals(k, sel)
        for (g, acc, pos, aa), j, p in zip(POSCTRL, sel, pp):
            hit = (k[j] >= KMIN) and (p < 0.05 / len(POSCTRL)); npass += hit
            prow.append(dict(gene=g, accession=acc, residue=f"{aa}{pos}", k=int(k[j]), expected=float(M.lam[j]), p=float(p), called=bool(hit)))
            log(f"  positive control {g} {aa}{pos}: k {k[j]}, expected {M.lam[j]:.3f}, p {p:.3g} -> {'called' if hit else 'not called'}")
        ok = npass >= 6; verdict["positive"] = ok; log(f"positive control: {npass} of {len(POSCTRL)} called (rule >= 6) -> {'PASSED' if ok else 'FAILED'}")
        write_tsv(T_POS, prow, log)
        # C3 planted (main, real data)
        p_real = M.pvals(k, fam1); fam_set_strata = M.st[fam1]; c3rows = []; t0 = time.time()
        ev_by_s = collections.defaultdict(list)
        for s, c in zip(M.ev_s.tolist(), M.ev_c.tolist()): ev_by_s[s].append(c)
        for r in range(1, N_C3 + 1):
            rng = np.random.default_rng(SEEDS["c3"] + r); pick = rng.choice(len(fam1), N_PLANT, replace=False); pj = fam1[pick]
            draws = [rng.integers(0, len(ev_by_s[int(M.st[j])]), max(DS)) for j in pj]
            for d in DS:
                k2 = k.copy(); n2 = M.n.copy(); aff = set()
                for j, dr in zip(pj, draws):
                    s = int(M.st[j]); k2[j] += d; aff.add(s)
                    for x in dr[:d]: n2[s, ev_by_s[s][x]] += 1
                M.prep_tail(n2)
                p2 = p_real.copy(); redo = np.flatnonzero(np.isin(fam_set_strata, list(aff)))
                p2[redo] = M.pvals(k2, fam1[redo], n2)
                q, hot = call(p2, k2, fam1); rec = float(hot[pick].mean())
                c3rows.append(dict(repeat=r, d=d, planted=N_PLANT, recovered=int(hot[pick].sum()), recovery=rec, other_hotspots=int(hot.sum() - hot[pick].sum())))
            M.prep_tail(M.n)
            if r in (1, N_C3): log(f"C3 repeat {r} done | elapsed {time.time() - t0:.0f} s")
        for d in DS:
            rr = [x["recovery"] for x in c3rows if x["d"] == d]; log(f"C3 planted d = {d}: mean recovery {np.mean(rr):.3f} (min {min(rr):.2f}, max {max(rr):.2f})")
        rec8 = float(np.mean([x["recovery"] for x in c3rows if x["d"] == 8])); ok = rec8 >= C3_MIN; verdict["C3"] = ok
        log(f"C3 planted: recovery at d = 8 {rec8:.3f} (rule >= {C3_MIN}) -> {'PASSED' if ok else 'FAILED'}")
        write_tsv(T_C3, c3rows, log)
    write_tsv(T_C2, c2rows, log)
    with open(VERDICT, "w") as o:
        for kk, v in verdict.items(): o.write(f"{kk}\t{'PASSED' if v else 'FAILED'}\n")
    log("verdicts: " + ", ".join(f"{kk} {'PASSED' if v else 'FAILED'}" for kk, v in verdict.items()) + f" | {VERDICT} md5 {H.md5_file(VERDICT)}")
    with open(OK2, "w") as o: o.write("OK\n")

def write_tsv(path, rows, log):
    if not rows: return
    cols = list(rows[0].keys())
    for r in rows:
        for c in r:
            if c not in cols: cols.append(c)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(cols)
        for r in rows: w.writerow([S8.fmt(r.get(c, "")) if isinstance(r.get(c, ""), float) else r.get(c, "") for c in cols])
    log(f"wrote {path} ({len(rows)} rows) md5 {H.md5_file(path)}")

def read_verdict():
    return {l.split("\t")[0]: l.split("\t")[1].strip() == "PASSED" for l in open(VERDICT) if "\t" in l}

def woolf(a, b, c, d):
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5; l = math.log(a * d / (b * c)); se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return math.exp(l), math.exp(l - 1.96 * se), math.exp(l + 1.96 * se)

# ======================= Section 3: real result =======================
def section3(log):
    need_ok(OK2, "section 2"); check_md5(log, MD5_IN); V = read_verdict()
    log("control verdicts read: " + ", ".join(f"{k} {'PASSED' if v else 'FAILED'}" for k, v in V.items()))
    D, Pz = load(log); U = H.load_uniprot_main(); accs = [str(x) for x in D["accs"]]; drv = bailey_genes(log, U)
    res = {}; summ = []
    for setup in ("main", "nohyper"):
        M = Model(D, Pz, setup); k = M.k_real; fam1 = M.family(1); fam0 = M.family(0)
        p1 = M.pvals(k, fam1); q1, h1 = call(p1, k, fam1); p0 = M.pvals(k, fam0); q0, h0 = call(p0, k, fam0)
        bycase = collections.defaultdict(collections.Counter)
        for j, c in zip(M.ev_j.tolist(), M.ev_c.tolist()):
            if M.role[j] == 1: bycase[j][M.projects[c]] += 1
        res[setup] = dict(M=M, fam1=fam1, p1=p1, q1=q1, h1=h1, fam0=fam0, p0=p0, q0=q0, h0=h0, bycase=bycase)
        n1, n0 = len(fam1), len(fam0); a = int(h1.sum()); b = int(h0.sum())
        a2 = int(((p1 < 0.001) & (k[fam1] >= KMIN)).sum()); b2 = int(((p0 < 0.001) & (k[fam0] >= KMIN)).sum())
        o1 = woolf(a, n1 - a, b, n0 - b); o2 = woolf(a2, n1 - a2, b2, n0 - b2)
        log(f"[{setup}] T1 family {n1:,} positions: Hotspot {a} | background family {n0:,}: Hotspot {b} | odds ratio T1 vs background {o1[0]:.2f} (95% CI {o1[1]:.2f} to {o1[2]:.2f})")
        log(f"[{setup}] nominal p < 0.001 and k >= 3: T1 {a2} ({a2 / n1:.5f}), background {b2} ({b2 / n0:.5f}); odds ratio {o2[0]:.2f} (95% CI {o2[1]:.2f} to {o2[2]:.2f})")
        log(f"[{setup}] T1 positions with k >= 1: {(k[fam1] >= 1).sum():,}, k >= 3: {(k[fam1] >= 3).sum():,}; observed / expected over the T1 family {k[fam1].sum() / M.lam[fam1].sum():.3f}")
        summ.append(dict(setup=setup, T1_positions=n1, T1_hotspots=a, bg_positions=n0, bg_hotspots=b, OR_hotspot=o1[0], OR_lo=o1[1], OR_hi=o1[2],
                         T1_nominal=a2, bg_nominal=b2, OR_nominal=o2[0], ORn_lo=o2[1], ORn_hi=o2[2], T1_obs=int(k[fam1].sum()), T1_exp=float(M.lam[fam1].sum())))
    # labels (main primary, nohyper beside)
    Mm = res["main"]["M"]; nh = res["nohyper"]; nh_hot = set(nh["fam1"][nh["h1"]].tolist())
    pos_nh = {int(j): i for i, j in enumerate(nh["fam1"])}
    def lab(setup, hot):
        if not V.get(f"C2_{setup}", False): return "Inconclusive (null control failed)"
        s = "Hotspot" if hot else ("Not hotspot" if V.get("C3", False) else "Not hotspot, test not sensitive")
        if not V.get("calibration", False): s += "_calibration_FAILED"
        if not V.get("positive", False): s += "_positive_control_FAILED"
        return s
    rows = []; r = res["main"]
    for i, j in enumerate(r["fam1"]):
        kj = int(Mm.k_real[j])
        if kj == 0: continue
        acc = accs[Mm.prot[j]]; g = U[acc]["gene"]; inh = pos_nh.get(int(j))
        hot = bool(r["h1"][i]); hot_nh = int(j) in nh_hot
        rows.append(dict(accession=acc, gene=g, position=int(Mm.posn[j]), residue=U[acc]["seq"][Mm.posn[j] - 1], k_main=kj, expected_main=float(Mm.lam[j]),
                         obs_exp=kj / Mm.lam[j] if Mm.lam[j] > 0 else float("inf"), p_main=float(r["p1"][i]), fdr_main=float(r["q1"][i]), label_main=lab("main", hot),
                         k_nohyper=int(nh["M"].k_real[j]), p_nohyper=float(nh["p1"][inh]) if inh is not None else 1.0,
                         fdr_nohyper=float(nh["q1"][inh]) if inh is not None else 1.0, label_nohyper=lab("nohyper", hot_nh),
                         hypermutator_driven=bool(hot and not hot_nh), bailey_driver=g in drv,
                         cases_by_cancer=",".join(f"{c.replace('TCGA-', '')}:{n}" for c, n in sorted(r["bycase"][int(j)].items()))))
    rows.sort(key=lambda x: (x["p_main"], -x["k_main"], x["gene"], x["position"]))
    write_tsv(T_HOT, rows, log)
    hot = [x for x in rows if x["label_main"].startswith("Hotspot")]
    log(f"T1 hotspots (main): {len(hot)}; also nohyper: {sum(1 for x in hot if not x['hypermutator_driven'])}; hypermutator-driven: {sum(1 for x in hot if x['hypermutator_driven'])}")
    for x in hot[:40]: log(f"  {x['gene']} {x['residue']}{x['position']} k {x['k_main']} exp {x['expected_main']:.3f} FDR {x['fdr_main']:.3g} | nohyper {x['label_nohyper']} | Bailey {x['bailey_driver']} | {x['cases_by_cancer']}")
    hn = [x for x in rows if x["label_nohyper"].startswith("Hotspot") and not x["label_main"].startswith("Hotspot")]
    log(f"nohyper Hotspots that are not main Hotspots: {len(hn)}" + "".join(f"\n  {x['gene']} {x['residue']}{x['position']} k_nohyper {x['k_nohyper']} FDR {x['fdr_nohyper']:.3g}" for x in hn[:20]))
    bg = []; r = res["main"]
    for i in np.flatnonzero(r["h0"]):
        j = int(r["fam0"][i]); acc = accs[Mm.prot[j]]; g = U[acc]["gene"]
        bg.append(dict(accession=acc, gene=g, position=int(Mm.posn[j]), residue=U[acc]["seq"][Mm.posn[j] - 1], k=int(Mm.k_real[j]), expected=float(Mm.lam[j]),
                       p=float(r["p0"][i]), fdr=float(r["q0"][i]), bailey_driver=g in drv))
    bg.sort(key=lambda x: (x["p"], -x["k"])); write_tsv(T_BGH, bg, log)
    log(f"background hotspots (main, descriptive): {len(bg)}, in Bailey driver genes {sum(x['bailey_driver'] for x in bg)}")
    hg = sorted({x["gene"] for x in hot}); bail = [dict(gene=g, bailey_driver=g in drv, n_hotspots=sum(1 for x in hot if x["gene"] == g)) for g in hg]
    write_tsv(T_BAIL, bail, log) if bail else log("no T1 hotspot gene, no Bailey table written")
    log(f"Bailey overlap (descriptive): {sum(x['bailey_driver'] for x in bail)} of {len(bail)} T1 hotspot genes are Bailey drivers")
    write_tsv(T_SUM, summ, log)
    with open(OK3, "w") as o: o.write("OK\n")

# ======================= Section 4: ActiveDriver inputs =======================
def section4(log):
    need_ok(OK1, "section 1"); check_md5(log, MD5_AD + MD5_IN[2:3])
    U = H.load_uniprot_main(); t0 = time.time()
    prot = {}
    with open(S8.PROTS, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"): prot[r["accession"]] = r
    t1 = collections.defaultdict(list); anyseq = collections.defaultdict(set)
    with open(S8.SITES, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = int(r["N_position"]); anyseq[r["accession"]].update((p, p + 1, p + 2))
            if r["tier"] == "T1": t1[r["accession"]].append(p)
    muts = collections.defaultdict(list)
    with gzip.open(S8.MUTS, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["in_chosen_file"] == "1": muts[r["accession"]].append((r["case_id"], int(r["position"]), r["ref"], r["alt"]))
    pos_acc = sorted({a for _, a, _, _ in POSCTRL})
    real = sorted(a for a in t1 if prot.get(a, {}).get("alphafold_status") == "use" and prot[a]["n_fragments"] == "1" and muts.get(a))
    need = set(real) | set(pos_acc)
    log(f"ActiveDriver proteins: T1 sequon {len(t1):,}; AlphaFold use with one fragment and >= 1 chosen-file missense {len(real):,}; positive-control proteins {pos_acc}")
    plddt = collections.defaultdict(dict); badaa = 0
    with gzip.open(S8.RESID, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            a = r["accession"]
            if a in need and r["fragment"] == "1" and r["plddt"] not in ("", "NA"):
                plddt[a][int(r["position"])] = float(r["plddt"])
                if U[a]["seq"][int(r["position"]) - 1] != r["aa"]: badaa += 1
    log(f"pLDDT rows read for {len(plddt):,} proteins; residue differs from UniProt: {badaa} (must be 0)")
    if badaa: raise RuntimeError("pLDDT residues do not match UniProt")
    drop = [a for a in real if len(plddt[a]) != len(U[a]["seq"])]
    real = [a for a in real if a not in drop]; log(f"proteins without a full pLDDT row set (left out): {len(drop)}; kept {len(real):,}")
    seqf = open(AD["sequences"], "w", newline="", encoding="utf-8"); disf = open(AD["disorder"], "w", newline="", encoding="utf-8")
    seqf.write("gene\tsequence\n"); disf.write("gene\tdisorder\n"); nd = 0
    for a in sorted(need):
        s = U[a]["seq"]
        if len(plddt.get(a, {})) == len(s): d = "".join("1" if plddt[a][i] < 50 else "0" for i in range(1, len(s) + 1))
        else: d = "0" * len(s); nd += 1; log(f"  {a}: no full pLDDT set, disorder written as all 0 (positive control only)")
        seqf.write(f"{a}\t{s}\n"); disf.write(f"{a}\t{d}\n")
    seqf.close(); disf.close()
    nm = 0
    with gzip.open(AD["mutations"], "wt", newline="", encoding="utf-8") as f:
        f.write("gene\tsample_id\tposition\twt_residue\tmut_residue\n")
        for a in sorted(need):
            for cid, p, ref, alt in muts.get(a, []): f.write(f"{a}\t{cid}\t{p}\t{ref}\t{alt}\n"); nm += 1
    log(f"mutation rows written {nm:,}")
    with open(AD["sites_real"], "w", newline="", encoding="utf-8") as f:
        f.write("gene\tposition\tresidue\tkinase\n")
        for a in real:
            for p in sorted(set(t1[a])): f.write(f"{a}\t{p}\tN\t\n")
    with open(AD["sites_pos"], "w", newline="", encoding="utf-8") as f:
        f.write("gene\tposition\tresidue\tkinase\n")
        for g, a, p, aa in POSCTRL: f.write(f"{a}\t{p}\t{aa}\t\n")
    rng = np.random.default_rng(SEEDS["adnull"]); sub = sorted(rng.choice(real, min(1000, len(real)), replace=False).tolist()); short = 0
    with open(AD["sites_null"], "w", newline="", encoding="utf-8") as f:
        f.write("rep\tgene\tposition\tresidue\tkinase\n")
        for rep in (1, 2, 3):
            rr = np.random.default_rng(SEEDS["adnull"] + rep)
            for a in sub:
                s = U[a]["seq"]; cand = [i + 1 for i, x in enumerate(s) if x == "N" and (i + 1) not in anyseq[a]]; m = len(set(t1[a]))
                if len(cand) < m: short += 1
                for p in sorted(rr.choice(cand, min(m, len(cand)), replace=False).tolist()) if cand else []: f.write(f"{rep}\t{a}\t{p}\tN\t\n")
    log(f"null subset {len(sub)} proteins x 3 repeats; protein-repeats with fewer non-sequon N than T1 sites (all taken): {short}")
    with open(AD["proteins"], "w", newline="", encoding="utf-8") as f:
        f.write("gene\tsymbol\tlength\tn_T1_sites\tn_mutations\tset\n")
        for a in sorted(need):
            st = ("real" if a in real else "") + (",null" if a in sub else "") + (",positive" if a in pos_acc else "")
            f.write(f"{a}\t{U[a]['gene']}\t{len(U[a]['seq'])}\t{len(set(t1.get(a, [])))}\t{len(muts.get(a, []))}\t{st.strip(',')}\n")
    for k_, pth in AD.items(): log(f"ActiveDriver input {k_}: {pth} md5 {H.md5_file(pth)}")
    log(f"section 4 {time.time() - t0:.0f} s")

SECTIONS = {"1": ("structure_check", section1), "2": ("controls", section2), "3": ("real_result", section3), "4": ("activedriver_inputs", section4)}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    C.run_chain("P2_20", SECTIONS, "1" if which == "all" else which)
# ---- END OF SCRIPT ----
